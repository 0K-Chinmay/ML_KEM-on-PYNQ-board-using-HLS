#include "unified.h"
// Matrix type for A matrix


// Vector NTT forward transform
void polyvec_ntt(polyvec_t* r) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        ntt_forward(&r->vec[i]);
    }
}

// Vector NTT i

// Vector addition
void polyvec_add(polyvec_t* r, const polyvec_t* a, const polyvec_t* b) {
#pragma HLS INLINE 
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL 
        poly_add(&r->vec[i], &a->vec[i], &b->vec[i]);
    }
}

// Vector subtraction
void polyvec_sub(polyvec_t* r, const polyvec_t* a, const polyvec_t* b) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_sub(&r->vec[i], &a->vec[i], &b->vec[i]);
    }
}

// Vector reduction
void polyvec_reduce(polyvec_t* r) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_reduce(&r->vec[i]);
    }
}

// Point-wise multiplication and accumulation
void polyvec_pointwise_acc_montgomery(poly_t* r, const polyvec_t* a, const polyvec_t* b) {
#pragma HLS INLINE off
    poly_t temp;
#pragma HLS ARRAY_PARTITION variable=temp.coeffs complete
    
    // Initialize result with first multiplication
    poly_basemul_montgomery(r, &a->vec[0], &b->vec[0]);
    
    // Accumulate remaining multiplications
    for (int i = 1; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_basemul_montgomery(&temp, &a->vec[i], &b->vec[i]);
        poly_add(r, r, &temp);
    }
}

// Matrix-vector multiplication: r = A * s
void matrix_vector_mul(polyvec_t* r, const matrix_t* A, const polyvec_t* s) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        polyvec_pointwise_acc_montgomery(&r->vec[i], &A->rows[i], s);
    }
}
// Field arithmetic

// Vector sampling with CBD eta1
void polyvec_cbd_eta1(polyvec_t* r, const byte_t* buf) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_cbd_eta1(&r->vec[i], buf + i * MLKEM_ETA1 * MLKEM_N / 4);
    }
}

// Vector sampling with CBD eta2
void polyvec_cbd_eta2(polyvec_t* r, const byte_t* buf) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_cbd_eta2(&r->vec[i], buf + i * MLKEM_ETA2 * MLKEM_N / 4);
    }
}

// Vector serialization
void polyvec_tobytes(byte_t* r, const polyvec_t* a) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_tobytes(r + i * MLKEM_POLYBYTES, &a->vec[i]);
    }
}

// Vector deserialization
void polyvec_frombytes(polyvec_t* r, const byte_t* a) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        poly_frombytes(&r->vec[i], a + i * MLKEM_POLYBYTES);
    }
}

// Matrix generation from seed
void matrix_expand(matrix_t* A, const byte_t* rho) {
#pragma HLS INLINE off
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        for (int j = 0; j < MLKEM_K; j++) {
#pragma HLS UNROLL
            poly_uniform(&A->rows[i].vec[j], rho, (byte_t)(i * MLKEM_K + j));
        }
    }
}

