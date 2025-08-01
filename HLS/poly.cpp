#include "unified.h"
// Montgomery reduction
coeff_t montgomery_reduce(int32_t a) {
#pragma HLS INLINE
    int32_t t = (int32_t)((ap_uint<16>)a * QINV);
    t = (a - t * MLKEM_Q) >> 16;
    return (coeff_t)t;
}

// Barrett reduction
coeff_t barrett_reduce(coeff_t a) {
#pragma HLS INLINE
    int32_t t = ((int32_t)a * 20159) >> 26;
    return a - t * MLKEM_Q;
}

// Conditional subtraction
coeff_t csubq(coeff_t a) {
#pragma HLS INLINE
    if(a >= MLKEM_Q)  
        return (a - MLKEM_Q); 
    else
        return a;
}

// NTT forward transform
ap_uint<8> reverse_bits(ap_uint<8> x, int bits) {
    ap_uint<8> y = 0;
    for (int i = 0; i < bits; i++) {
        y = (y << 1) | (x & 1);
        x >>= 1;
    }
    return y;
}
const int N=256;
const int MOD = 3329;
// NTT forward transform

void ntt_forward(poly_t* r) {
   int k = 1;
    
    for (int l = 128; l >= 2; l = l >> 1) {        
        for (int start = 0; start < 256; start = start + 2 * l) {
            #pragma HLS UNROLL
            ap_uint<16> zeta = ntt_zetas[k];
            k = k + 1;
            
            for (int j = start; j < start + l; j++) {
                #pragma HLS UNROLL
                ap_uint<32> t = (ap_uint<32>)zeta * r->coeffs[j + l];
                ap_uint<16> t_mod = t % MOD;
                
                ap_uint<32> diff = (ap_uint<32>)r->coeffs[j] + MOD - t_mod;
                r->coeffs[j + l] = diff % MOD;
                
                ap_uint<32> sum = (ap_uint<32>)r->coeffs[j] + t_mod;
                r->coeffs[j] = sum % MOD;
            }
        }
    }
    
    // Final modular reduction (optional since we're doing it in the loop)
    for (int j = 0; j < 256; j++) {
        r->coeffs[j] = barrett_reduce( r->coeffs[j]) ;
    }

}


// Base multiplication function for 2 coefficients and zeta
void ntt_base_multiplication(int16_t *r0, int16_t *r1,
                                           int16_t a0, int16_t a1,
                                           int16_t b0, int16_t b1,
                                           int16_t zeta) {
    int32_t t0 = (int32_t)a0 * b0;
    int32_t t1 = (int32_t)zeta * a1 % MLKEM_Q;
    t1 = t1 * b1 % MLKEM_Q;
    *r0 = (t0 + t1) % MLKEM_Q;

    int32_t t2 = (int32_t)a1 * b0 + (int32_t)a0 * b1;
    *r1 = t2 % MLKEM_Q;
}

// Perform base-wise multiplication using NTT with zeta coefficients
void poly_basemul_montgomery(poly_t *r, const poly_t *a, const poly_t *b) {
    for (int i = 0; i < 64; i++) {
        int16_t a0 = a->coeffs[4 * i + 0];
        int16_t a1 = a->coeffs[4 * i + 1];
        int16_t a2 = a->coeffs[4 * i + 2];
        int16_t a3 = a->coeffs[4 * i + 3];

        int16_t b0 = b->coeffs[4 * i + 0];
        int16_t b1 = b->coeffs[4 * i + 1];
        int16_t b2 = b->coeffs[4 * i + 2];
        int16_t b3 = b->coeffs[4 * i + 3];

        int16_t r0, r1, r2, r3;

        // first pair
        ntt_base_multiplication(&r0, &r1, a0, a1, b0, b1, ntt_zetas[64 + i]);

        // second pair with -zeta
        ntt_base_multiplication(&r2, &r3, a2, a3, b2, b3, (MLKEM_Q - ntt_zetas[64 + i]) % MLKEM_Q);

        r->coeffs[4 * i + 0] = r0;
        r->coeffs[4 * i + 1] = r1;
        r->coeffs[4 * i + 2] = r2;
        r->coeffs[4 * i + 3] = r3;
    }
}



// Polynomial addition
void poly_add(poly_t* r, const poly_t* a, const poly_t* b) {
#pragma HLS INLINE 
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete
#pragma HLS ARRAY_PARTITION variable=a->coeffs complete
#pragma HLS ARRAY_PARTITION variable=b->coeffs complete

    for (int i = 0; i < MLKEM_N; i++) {        
//#pragma HLS PIPELINE 
#pragma HLS UNROLL
        r->coeffs[i] = barrett_reduce(a->coeffs[i] + b->coeffs[i]);
    }
}

// Polynomial subtraction
void poly_sub(poly_t* r, const poly_t* a, const poly_t* b) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete
#pragma HLS ARRAY_PARTITION variable=a->coeffs complete
#pragma HLS ARRAY_PARTITION variable=b->coeffs complete

    for (int i = 0; i < MLKEM_N; i++) {

#pragma HLS PIPELINE II=1
        r->coeffs[i] = a->coeffs[i] - b->coeffs[i] + MLKEM_Q;
    }
}

// Reduce polynomial coefficients modulo q
void poly_reduce(poly_t* r) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete

    for (int i = 0; i < MLKEM_N; i++) {
#pragma HLS PIPELINE II=1
        r->coeffs[i] = barrett_reduce(r->coeffs[i]);
    }
}

// Centered binomial distribution sampling

void poly_cbd_eta1(poly_t* r, const byte_t* buf) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete
#pragma HLS ARRAY_PARTITION variable=buf complete
//print_hex(buf, 64 * MLKEM_ETA1 , "");

    for (int i = 0; i < MLKEM_N / 4; i++) {
#pragma HLS PIPELINE II=1
        // Combine 3 bytes into 24 bits
        uint32_t t = buf[3*i];
        t |= ((uint32_t)buf[3*i + 1]) << 8;
        t |= ((uint32_t)buf[3*i + 2]) << 16;

        for (int j = 0; j < 4; j++) {
#pragma HLS UNROLL
            uint32_t d = t & 0x3F;  // extract 6 bits for one coefficient

            // Count a = number of 1s in first 3 bits
            uint32_t a = (d >> 0 & 0x1) + (d >> 1 & 0x1) + (d >> 2 & 0x1);
            // Count b = number of 1s in next 3 bits
            uint32_t b = (d >> 3 & 0x1) + (d >> 4 & 0x1) + (d >> 5 & 0x1);

            r->coeffs[4*i + j] = (coeff_t)(a - b + MLKEM_Q); // mod wraparound-safe

            t >>= 6; // move to next 6 bits
        }
    }
}


void poly_cbd_eta2(poly_t* r, const byte_t* buf) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete
#pragma HLS ARRAY_PARTITION variable=buf complete

    for (int i = 0; i < MLKEM_N; i++) {
#pragma HLS PIPELINE II=1
        uint32_t t = 0;
        int byte_pos = i / 4;  // Each coefficient uses 4 bits from eta=2
        int bit_pos = (i % 4) * 2;
        
        // Extract 4 bits for CBD(eta=2)
        t = (buf[byte_pos] >> bit_pos) & 0x0F;
        
        // Compute a - b where a and b are 2-bit values
        uint32_t a = (t & 0x03) + ((t >> 2) & 0x03);
        uint32_t b = ((t >> 1) & 0x01) + ((t >> 3) & 0x01);
        r->coeffs[i] = (coeff_t)(a - b + MLKEM_Q);
    }
}

// Serialize polynomial to bytes
void poly_tobytes(byte_t* r, const poly_t* a) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r complete
#pragma HLS ARRAY_PARTITION variable=a->coeffs complete

    for (int i = 0; i < MLKEM_N; i += 2) {
#pragma HLS PIPELINE II=1
        coeff_t t0 = csubq(a->coeffs[i]);
        coeff_t t1 = csubq(a->coeffs[i + 1]);
        
        r[3 * i / 2] = (byte_t)t0;
        r[3 * i / 2 + 1] = (byte_t)((t0 >> 8) | (t1 << 4));
        r[3 * i / 2 + 2] = (byte_t)(t1 >> 4);
    }
}

// Deserialize polynomial from bytes
void poly_frombytes(poly_t* r, const byte_t* a) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete
#pragma HLS ARRAY_PARTITION variable=a complete



    for (int i = 0; i < MLKEM_N; i += 2) {
#pragma HLS PIPELINE II=1
        r->coeffs[i] = ((coeff_t)a[3 * i / 2] | ((coeff_t)a[3 * i / 2 + 1] << 8)) & 0xFFF;
        r->coeffs[i + 1] = ((coeff_t)a[3 * i / 2 + 1] >> 4 | ((coeff_t)a[3 * i / 2 + 2] << 4)) & 0xFFF;
    }
}

// Uniform sampling from XOF
void poly_uniform(poly_t* r, const byte_t* seed, byte_t nonce) {
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=r->coeffs complete

    byte_t input[34];
#pragma HLS ARRAY_PARTITION variable=input complete
    
    // Copy seed
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        input[i] = seed[i];
    }
    input[33] = (nonce == 2) | (nonce == 3);
    input[32] = nonce % 2;
   // printf("\n%x\n",nonce);
    //print_poly(r[0]);
    // Generate random bytes using SHAKE128
    byte_t buf[REJ_UNIFORM_BUFLEN];
    
#pragma HLS ARRAY_PARTITION variable=buf complete

shake128(input, 34, buf, REJ_UNIFORM_BUFLEN);

    // Rejection sampling
    int ctr = 0;
    for (int i = 0; i < REJ_UNIFORM_BUFLEN && ctr < MLKEM_N; i += 3) {
#pragma HLS PIPELINE II=1
        coeff_t val1 = ((coeff_t)buf[i] | ((coeff_t)buf[i + 1] << 8)) & 0xFFF;
        coeff_t val2 = ((coeff_t)buf[i + 1] >> 4 | ((coeff_t)buf[i + 2] << 4)) & 0xFFF;
        
        if (val1 < MLKEM_Q && ctr < MLKEM_N) {
            r->coeffs[ctr++] = val1;
        }
        if (val2 < MLKEM_Q && ctr < MLKEM_N) {
            r->coeffs[ctr++] = val2;
        }
    }
}

