#include "unified.h"

// Keccak-f[1600] implementation (from PRF module)


// PRF function implementation
void prf_eta(int eta, const byte_t s[32], byte_t b, byte_t* output) {
#pragma HLS INLINE off
    
    // Input buffer: s (32 bytes) || b (1 byte)
    byte_t input[33];
#pragma HLS ARRAY_PARTITION variable=input complete
    
    // Copy s into input buffer
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        input[i] = s[i];
    }
    
    // Append b
    input[32] = b;
    
    // Output length is 64 * eta bytes
    int output_len = 64 * eta;
    //print_hex(input, 33, "");
    // Call SHAKE256
    shake256(input, 33, output, output_len);
}

