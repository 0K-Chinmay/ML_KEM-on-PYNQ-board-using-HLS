
#include "unified.h"
#include <stdio.h>

const lane_t RC[24] = {
    0x0000000000000001ULL,0x0000000000008082ULL,0x800000000000808AULL,0x8000000080008000ULL,
    0x000000000000808BULL,0x0000000080000001ULL,0x8000000080008081ULL,0x8000000000008009ULL,
    0x000000000000008AULL,0x0000000000000088ULL,0x0000000080008009ULL,0x000000008000000AULL,
    0x000000008000808BULL,0x800000000000008BULL,0x8000000000008089ULL,0x8000000000008003ULL,
    0x8000000000008002ULL,0x8000000000000080ULL,0x000000000000800AULL,0x800000008000000AULL,
    0x8000000080008081ULL,0x8000000000008080ULL,0x0000000080000001ULL,0x8000000080008008ULL
};

    // Rotation offsets for rho step
    const int rho_offsets[25] = {
        0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56, 14
    };

void keccak_f1600(lane_t state[25]) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=state complete
    
    lane_t C[5], D[5], B[25];
#pragma HLS ARRAY_PARTITION variable=C complete
#pragma HLS ARRAY_PARTITION variable=D complete
#pragma HLS ARRAY_PARTITION variable=B complete
    
    // Round constants for Keccak-f[1600]

    
    for (int round = 0; round < 24; round++) {
#pragma HLS LOOP_TRIPCOUNT min=24 max=24
#pragma HLS PIPELINE II=1
        
        // Theta step
        for (int x = 0; x < 5; x++) {
#pragma HLS UNROLL
            C[x] = state[x] ^ state[x + 5] ^ state[x + 10] ^ state[x + 15] ^ state[x + 20];
        }

        for (int x = 0; x < 5; x++) {
#pragma HLS UNROLL
            D[x] = C[(x + 4) % 5] ^ ((C[(x + 1) % 5] << 1) | (C[(x + 1) % 5] >> 63));
        }

        for (int x = 0; x < 5; x++) {
#pragma HLS UNROLL
            for (int y = 0; y < 5; y++) {
#pragma HLS UNROLL
                state[5 * y + x] ^= D[x];
            }
        }



        // Rho and Pi steps
        for (int i = 0; i < 25; i++) {
        #pragma HLS UNROLL
            int x = i % 5;
            int y = i / 5;
            
            // Pi mapping
            int new_x = y;
            int new_y = (2 * x + 3 * y) % 5;
            int pi_index = new_x + 5 * new_y;

            // Rho offset
            int rho_offset = rho_offsets[i];

            if (rho_offset == 0) {
                B[pi_index] = state[i];
            } else {
                B[pi_index] = (state[i] << rho_offset) | (state[i] >> (64 - rho_offset));
            }
        }


        // Chi step
        for (int y = 0; y < 5; y++) {
#pragma HLS UNROLL
            for (int x = 0; x < 5; x++) {
#pragma HLS UNROLL
                state[5 * y + x] = B[5 * y + x] ^ ((~B[5 * y + (x + 1) % 5]) & B[5 * y + (x + 2) % 5]);
            }
        }


        
        // Iota step
        state[0] ^= RC[round];

    }
}
// SHA3-256 hash function
void shake256(const byte_t* input, int input_len, byte_t* output, int output_len) {
#pragma HLS INTERFACE m_axi port=input offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=output offset=slave bundle=gmem1
#pragma HLS INTERFACE s_axilite port=input_len bundle=control
#pragma HLS INTERFACE s_axilite port=output_len bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    lane_t state[25];
#pragma HLS ARRAY_PARTITION variable=state complete
    
    // Initialize state
    for (int i = 0; i < 25; i++) {
#pragma HLS UNROLL
        state[i] = 0;
    }
    
    // Absorbing phase
    int pos = 0;
    
    // Process full blocks
    while (pos + SHAKE256_RATE <= input_len) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=10
        
        // XOR input block into state
        for (int i = 0; i < SHAKE256_RATE; i += 8) {
#pragma HLS PIPELINE II=1
            lane_t block_word = 0;
            for (int j = 0; j < 8 && (i + j) < SHAKE256_RATE; j++) {
#pragma HLS UNROLL
                block_word |= (lane_t)input[pos + i + j] << (8 * j);
            }
            state[i / 8] ^= block_word;
        }
        
        keccak_f1600(state);
        pos += SHAKE256_RATE;
    }
    
    // Process remaining bytes and padding
    byte_t last_block[SHAKE256_RATE];
#pragma HLS ARRAY_PARTITION variable=last_block complete
    
    // Copy remaining input bytes
    int remaining = input_len - pos;
    for (int i = 0; i < remaining; i++) {
#pragma HLS UNROLL
        last_block[i] = input[pos + i];
    }
    
    // Add padding (0x1f for SHAKE256)
    last_block[remaining] = 0x1f;
    for (int i = remaining + 1; i < SHAKE256_RATE - 1; i++) {
#pragma HLS UNROLL
        last_block[i] = 0;
    }
    last_block[SHAKE256_RATE - 1] = 0x80;
    
    // XOR last block into state
    for (int i = 0; i < SHAKE256_RATE; i += 8) {
#pragma HLS PIPELINE II=1
        lane_t block_word = 0;
        for (int j = 0; j < 8; j++) {
#pragma HLS UNROLL
            block_word |= (lane_t)last_block[i + j] << (8 * j);
        }
        state[i / 8] ^= block_word;
    }
    
    keccak_f1600(state);
    
    // Squeezing phase
    int output_pos = 0;
    
    while (output_pos < output_len) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=10
        
        // Extract bytes from current state
        int extract_len = (output_len - output_pos < SHAKE256_RATE) ? 
                         (output_len - output_pos) : SHAKE256_RATE;
        
        for (int i = 0; i < extract_len; i += 8) {
#pragma HLS PIPELINE II=1
            lane_t lane_data = state[i / 8];
            for (int j = 0; j < 8 && (i + j) < extract_len; j++) {
#pragma HLS UNROLL
                output[output_pos + i + j] = (byte_t)(lane_data >> (8 * j));
            }
        }
        
        output_pos += extract_len;
        
        // If more output needed, permute state
        if (output_pos < output_len) {
            keccak_f1600(state);
        }
    }
}

void sha3_256(const byte_t* input, int input_len, byte_t output[32]) {
#pragma HLS INTERFACE m_axi port=input offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=output offset=slave bundle=gmem1
#pragma HLS INTERFACE s_axilite port=input_len bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    lane_t state[25];
#pragma HLS ARRAY_PARTITION variable=state complete
    
    // Initialize state
    for (int i = 0; i < 25; i++) {
#pragma HLS UNROLL
        state[i] = 0;
    }
    
    // Absorbing phase (similar to SHAKE256 but with different rate)
    int rate = 136; // SHA3-256 rate
    int pos = 0;
    
    // Process full blocks
    while (pos + rate <= input_len) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=10
        
        // XOR input block into state
        for (int i = 0; i < rate; i += 8) {
#pragma HLS PIPELINE II=1
            lane_t block_word = 0;
            for (int j = 0; j < 8 && (i + j) < rate; j++) {
#pragma HLS UNROLL
                block_word |= (lane_t)input[pos + i + j] << (8 * j);
            }
            state[i / 8] ^= block_word;
        }
        
        keccak_f1600(state);
        pos += rate;
    }
    
    // Process remaining bytes and padding
    byte_t last_block[136];
#pragma HLS ARRAY_PARTITION variable=last_block complete
    
    // Copy remaining input bytes
    int remaining = input_len - pos;
    for (int i = 0; i < remaining; i++) {
#pragma HLS UNROLL
        last_block[i] = input[pos + i];
    }
    
    // Add padding (0x06 for SHA3-256)
    last_block[remaining] = 0x06;
    for (int i = remaining + 1; i < rate - 1; i++) {
#pragma HLS UNROLL
        last_block[i] = 0;
    }
    last_block[rate - 1] = 0x80;
    
    // XOR last block into state
    for (int i = 0; i < rate; i += 8) {
#pragma HLS PIPELINE II=1
        lane_t block_word = 0;
        for (int j = 0; j < 8; j++) {
#pragma HLS UNROLL
            block_word |= (lane_t)last_block[i + j] << (8 * j);
        }
        state[i / 8] ^= block_word;
    }
    
    keccak_f1600(state);
    
    // Extract 256 bits (32 bytes)
    for (int i = 0; i < 32; i += 8) {
#pragma HLS PIPELINE II=1
        lane_t lane_data = state[i / 8];
        for (int j = 0; j < 8; j++) {
#pragma HLS UNROLL
            output[i + j] = (byte_t)(lane_data >> (8 * j));
        }
    }
}

void sha3_512(const byte_t* input, int input_len, byte_t output[64]) {
#pragma HLS INTERFACE m_axi port=input offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=output offset=slave bundle=gmem1
#pragma HLS INTERFACE s_axilite port=input_len bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    lane_t state[25];
#pragma HLS ARRAY_PARTITION variable=state complete

    // Initialize state to zero
    for (int i = 0; i < 25; i++) {
#pragma HLS UNROLL
        state[i] = 0;
    }

    const int rate = 72; // 72 bytes = 576 bits for SHA3-512
    int pos = 0;

    // Final block and padding
    byte_t last_block[rate];
#pragma HLS ARRAY_PARTITION variable=last_block complete

    // 1. Zero the block
    for (int i = 0; i < rate; i++) {
#pragma HLS UNROLL
        last_block[i] = 0;
    }

    // 2. Copy remaining bytes
    int remaining = input_len - pos;
    for (int i = 0; i < remaining; i++) {
#pragma HLS UNROLL
        last_block[i] = input[pos + i];
    }
 

    // 3. Padding: 0x06 and 0x80
    last_block[remaining] ^= 0x06;
    last_block[rate - 1] ^= 0x80;



    // 4. Absorb final padded block
    for (int i = 0; i < rate; i += 8) {
#pragma HLS PIPELINE II=1
        ap_uint<64> block_word = 0;
        for (int j = 0; j < 8; j++) {
#pragma HLS UNROLL
            block_word |= (ap_uint<64>)last_block[i + j] << (8 * j);
        }
        
    // Compute original linear index
    int idx = i / 8;
    state[idx] ^= block_word;
         //printf("\n");
    }
          
    keccak_f1600(state);

    // Extract 64 bytes (SHA3-512)
    for (int i = 0; i < 64; i++) {
#pragma HLS PIPELINE II=1
        int lane_index = i / 8;
        int byte_offset = i % 8;
        output[i] = (byte_t)(state[lane_index] >> (8 * byte_offset));
    }
}



// SHAKE128 implementation
void shake128(const byte_t* input, int input_len, byte_t* output, int output_len) {
#pragma HLS INTERFACE m_axi port=input offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=output offset=slave bundle=gmem1
#pragma HLS INTERFACE s_axilite port=input_len bundle=control
#pragma HLS INTERFACE s_axilite port=output_len bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    lane_t state[25];
#pragma HLS ARRAY_PARTITION variable=state complete
    
    // Initialize state
    for (int i = 0; i < 25; i++) {
#pragma HLS UNROLL
        state[i] = 0;
    }
    
    // Absorbing phase
    int pos = 0;
    
    // Process full blocks
    while (pos + SHAKE128_RATE <= input_len) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=10
        
        // XOR input block into state
        for (int i = 0; i < SHAKE128_RATE; i += 8) {
#pragma HLS PIPELINE II=1
            lane_t block_word = 0;
            for (int j = 0; j < 8 && (i + j) < SHAKE128_RATE; j++) {
#pragma HLS UNROLL
                block_word |= (lane_t)input[pos + i + j] << (8 * j);
            }
            state[i / 8] ^= block_word;
        }
        
        keccak_f1600(state);
        pos += SHAKE128_RATE;
    }
    
    // Process remaining bytes and padding
    byte_t last_block[SHAKE128_RATE];
#pragma HLS ARRAY_PARTITION variable=last_block complete
    
    // Copy remaining input bytes
    int remaining = input_len - pos;
    for (int i = 0; i < remaining; i++) {
#pragma HLS UNROLL
        last_block[i] = input[pos + i];
    }
    
    // Add padding (0x1f for SHAKE128)
    last_block[remaining] = 0x1f;
    for (int i = remaining + 1; i < SHAKE128_RATE - 1; i++) {
#pragma HLS UNROLL
        last_block[i] = 0;
    }
    last_block[SHAKE128_RATE - 1] = 0x80;
    
    // XOR last block into state
    for (int i = 0; i < SHAKE128_RATE; i += 8) {
#pragma HLS PIPELINE II=1
        lane_t block_word = 0;
        for (int j = 0; j < 8; j++) {
#pragma HLS UNROLL
            block_word |= (lane_t)last_block[i + j] << (8 * j);
        }
        state[i / 8] ^= block_word;
    }
    
    keccak_f1600(state);
    
    // Squeezing phase
    int output_pos = 0;
    
    while (output_pos < output_len) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=10
        
        // Extract bytes from current state
        int extract_len = (output_len - output_pos < SHAKE128_RATE) ? 
                         (output_len - output_pos) : SHAKE128_RATE;
        
        for (int i = 0; i < extract_len; i += 8) {
#pragma HLS PIPELINE II=1
            lane_t lane_data = state[i / 8];
            for (int j = 0; j < 8 && (i + j) < extract_len; j++) {
#pragma HLS UNROLL
                output[output_pos + i + j] = (byte_t)(lane_data >> (8 * j));
            }
        }
        
        output_pos += extract_len;
        
        // If more output needed, permute state
        if (output_pos < output_len) {
            keccak_f1600(state);
        }
    }
}

void G(const byte_t* input, int input_len, byte_t output[64]) {
#pragma HLS INLINE off
    sha3_512(input, input_len, output);
}

// H function: SHA3-256(input)
void H(const byte_t* input, int input_len, byte_t output[32]) {
#pragma HLS INLINE off
    sha3_256(input, input_len, output);
}

// XOF function: SHAKE128(input, output_len)
void XOF(const byte_t* input, int input_len, byte_t* output, int output_len) {
#pragma HLS INLINE off
    shake128(input, input_len, output, output_len);
}