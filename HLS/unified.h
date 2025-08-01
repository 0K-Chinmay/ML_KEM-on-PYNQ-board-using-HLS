#ifndef MLKEM_UNIFIED_H
#define MLKEM_UNIFIED_H

#include "ap_int.h"
#include <stdint.h>
#include <cstring>
#include <iostream>
#include <iomanip>
// ============================================================================
// ML-KEM-512 PARAMETERS AND CONSTANTS
// ============================================================================

// ML-KEM-512 parameters
const int MLKEM_N = 256;           // Polynomial degree
const int MLKEM_Q = 3329;          // Modulus
const int MLKEM_K = 2;             // Matrix dimension
const int MLKEM_ETA1 = 3;          // Noise parameter 1
const int MLKEM_ETA2 = 2;          // Noise parameter 2
const int MLKEM_DU = 10;           // Compression parameter u
const int MLKEM_DV = 4;            // Compression parameter v

// Key sizes
const int MLKEM_SYMBYTES = 32;     // Size of symmetric key
const int MLKEM_SSBYTES = 32;      // Size of shared secret
const int MLKEM_POLYBYTES = 384;   // Size of polynomial in bytes
const int MLKEM_POLYCOMPRESSEDBYTES_DU = 320;  // Compressed polynomial size for du=10
const int MLKEM_POLYCOMPRESSEDBYTES_DV = 128;  // Compressed polynomial size for dv=4
// Basic types
typedef ap_uint<8> byte_t;
typedef ap_uint<16> coeff_t;       // Coefficient type (can hold values up to q-1)
typedef ap_uint<64> lane_t;
    // For Keccak permutation

// Public key size: k * polybytes + 32
const int MLKEM_PUBLICKEYBYTES = MLKEM_K * MLKEM_POLYBYTES + MLKEM_SYMBYTES;  // 800 bytes
// Private key size: k * polybytes + publickey + 32 + 32
const int MLKEM_SECRETKEYBYTES = MLKEM_K * MLKEM_POLYBYTES + MLKEM_PUBLICKEYBYTES + MLKEM_SYMBYTES + MLKEM_SYMBYTES;  // 1632 bytes

// Constants for SHAKE128
const int SHAKE128_RATE = 168;       // 1344 bits / 8 = 168 bytes
const int SHAKE128_CAPACITY = 32;    // 256 bits / 8 = 32 bytes
const int SHAKE256_RATE = 136; // 1088 bits / 8 = 136 bytes
const int SHAKE256_CAPACITY = 64;
// Rejection sampling bounds
const int REJ_UNIFORM_BUFLEN = 504;  // Buffer length for uniform sampling
const int REJ_UNIFORM_ETA_BUFLEN = 256;  // Buffer length for CBD sampling

// NTT constants
const int NTT_ZETAS_SIZE = 128;
const coeff_t ntt_zetas[NTT_ZETAS_SIZE] = {1, 1729, 2580, 3289, 2642, 630, 1897, 848, 1062, 1919, 193, 797, 2786, 3260, 569, 1746, 296, 2447, 1339, 
1476, 3046, 56, 2240, 1333, 1426, 2094, 535, 2882, 2393, 2879, 1974, 821, 289, 331, 3253, 1756, 1197, 2304, 2277, 2055, 650, 1977, 2513, 632, 2865,
 33, 1320, 1915, 2319, 1435, 807, 452, 1438, 2868, 1534, 2402, 2647, 2617, 1481, 648, 2474, 3110, 1227, 910, 17, 2761, 583, 2649, 1637, 723, 2288, 
 1100, 1409, 2662, 3281, 233, 756, 2156, 3015, 3050, 1703, 1651, 2789, 1789, 1847, 952, 1461, 2687, 939, 2308, 2437, 2388, 733, 2337, 268, 641, 1584,
  2298, 2037, 3220, 375, 2549, 2090, 1645, 1063, 319, 2773, 757, 2099, 561, 2466, 2594, 2804, 1092, 403, 1026, 1143, 2150, 2775, 886, 1722, 1212, 1874,
   1029, 2110, 2935, 885, 2154};


// Montgomery reduction constants
const uint16_t QINV = 62209;  // q^(-1) mod 2^16
const uint16_t MONT = 2285;   // 2^16 mod q

// ============================================================================
// TYPE DEFINITIONS
// ============================================================================

void print_hex(const byte_t* data, int len, const std::string& label) ;
// Polynomial type
struct poly_t {
    coeff_t coeffs[MLKEM_N];
};

// Vector of polynomials
struct polyvec_t {
    poly_t vec[MLKEM_K];
};

// Matrix type
struct matrix_t {
    polyvec_t rows[MLKEM_K];
};
 void print_poly(const poly_t pv);

// ============================================================================
// POLYNOMIAL OPERATIONS
// ============================================================================

// Reduction functions
coeff_t montgomery_reduce(int32_t a);
coeff_t barrett_reduce(coeff_t a);
coeff_t csubq(coeff_t a);

// NTT operations
void ntt_forward(poly_t* r);
void ntt_inverse(poly_t* r);

// Polynomial arithmetic
void poly_basemul_montgomery(poly_t* r, const poly_t* a, const poly_t* b);
void poly_add(poly_t* r, const poly_t* a, const poly_t* b);
void poly_sub(poly_t* r, const poly_t* a, const poly_t* b);
void poly_reduce(poly_t* r);

// Polynomial sampling
void poly_cbd_eta1(poly_t* r, const byte_t* buf);
void poly_cbd_eta2(poly_t* r, const byte_t* buf);
void poly_uniform(poly_t* r, const byte_t* seed, byte_t nonce);

// Polynomial serialization
void poly_tobytes(byte_t* r, const poly_t* a);
void poly_frombytes(poly_t* r, const byte_t* a);

// ============================================================================
// POLYNOMIAL VECTOR OPERATIONS
// ============================================================================

// NTT operations on vectors
void polyvec_ntt(polyvec_t* r);


// Vector arithmetic
void polyvec_add(polyvec_t* r, const polyvec_t* a, const polyvec_t* b);
void polyvec_sub(polyvec_t* r, const polyvec_t* a, const polyvec_t* b);
void polyvec_reduce(polyvec_t* r);
void polyvec_pointwise_acc_montgomery(poly_t* r, const polyvec_t* a, const polyvec_t* b);
void ntt_base_multiplication(int16_t *r0, int16_t *r1,
                                           int16_t a0, int16_t a1,
                                           int16_t b0, int16_t b1,
                                           int16_t zeta);
// Matrix-vector operations
void matrix_vector_mul(polyvec_t* r, const matrix_t* A, const polyvec_t* s);
void matrix_transpose_vector_mul(polyvec_t* r, const matrix_t* A, const polyvec_t* s);

// Matrix manipulation
void matrix_expand(matrix_t* A, const byte_t* rho);
void matrix_transpose(matrix_t* At, const matrix_t* A);

// Vector sampling
void polyvec_cbd_eta1(polyvec_t* r, const byte_t* buf);
void polyvec_cbd_eta2(polyvec_t* r, const byte_t* buf);

// Vector serialization
void polyvec_tobytes(byte_t* r, const polyvec_t* a);
void polyvec_frombytes(polyvec_t* r, const byte_t* a);

// ============================================================================
// CRYPTOGRAPHIC HASH FUNCTIONS
// ============================================================================

// Keccak permutation
void keccak_f1600(lane_t state[25]);

// SHA3-256 hash
void sha3_256(const byte_t* input, int input_len, byte_t output[32]);

void sha3_512(const byte_t* input, int input_len, byte_t output[64]);

// SHAKE128 XOF
void shake128(const byte_t* input, int input_len, byte_t* output, int output_len);

// SHAKE256 XOF
void shake256(const byte_t* input, int input_len, byte_t* output, int output_len);

// G function = SHA3-256(input)
void G(const byte_t* input, int input_len, byte_t output[64]);

// H function = SHA3-256(input)
void H(const byte_t* input, int input_len, byte_t output[32]);

// XOF function = SHAKE128(input, output_len)
void XOF(const byte_t* input, int input_len, byte_t* output, int output_len);

// Pseudo-random function with domain separation
void prf_eta(int eta, const byte_t s[32], byte_t b, byte_t* output);

// ============================================================================
// KEY GENERATION FUNCTIONS
// ============================================================================



// Top-level key generation function for HLS
void mlkem512_keygen_top(const byte_t seed[32],const byte_t z[32], byte_t pk[MLKEM_PUBLICKEYBYTES], byte_t sk[MLKEM_SECRETKEYBYTES]);







// ============================================================================
// INLINE HELPER FUNCTIONS
// ============================================================================

// Conditional subtraction of q
inline coeff_t csubq_inline(coeff_t a) {
    if(a >= MLKEM_Q) 
        return (a - MLKEM_Q) ;
    else 
        return a;
}

// Modular reduction
inline coeff_t mod_q(int32_t a) {
    return (coeff_t)(a % MLKEM_Q);
}

// Freeze to [0, q-1]
inline coeff_t freeze(coeff_t a) {
    return csubq_inline(a);
}

#endif // MLKEM_UNIFIED_H