#include "unified.h"
#include <iomanip>
#include <iostream>
// Key generation function
/*void print_hex(const byte_t* data, int len, const std::string& label) {
    std::cout << label << ": ";
    for (int i = 0; i < len; i++) {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)data[i];
    
    }
    std::cout << std::dec << std::endl;
}

void print_poly(const poly_t pv) {

        for (int j = 0; j < MLKEM_N; j++) {
            printf("%d ", (pv.coeffs[j].to_uint()) % 3329);
            if ((j + 1) % 16 == 0) printf("\n"); // Format nicely in rows
        }

}

void print_polyvec(const polyvec_t pv) {
    for (int i = 0; i < MLKEM_K; i++) {
        printf("Polynomial %d:\n", i);
        for (int j = 0; j < MLKEM_N; j++) {
            printf("%u ", (pv.vec[i].coeffs[j].to_uint()) % 3329);
            if ((j + 1) % 16 == 0) printf("\n"); // Format nicely in rows
        }
        printf("\n");
    }
}

void print_matrix(const matrix_t A) {
    for(int k = 0 ; k <MLKEM_K;k++)
    for (int i = 0; i < MLKEM_K; i++) {
        printf("Polynomial %d:\n", i);
        for (int j = 0; j < MLKEM_N; j++) {
            printf("%u ", (A.rows[k].vec[i].coeffs[j].to_uint()) % 3329);
            if ((j + 1) % 16 == 0) printf("\n"); // Format nicely in rows
        }
        printf("\n");
    }
}*/

void mlkem512_keygen_top(const byte_t d[32],const byte_t z[32], byte_t pk[MLKEM_PUBLICKEYBYTES], byte_t sk[MLKEM_SECRETKEYBYTES]) {
#pragma HLS INTERFACE m_axi port=z offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=d offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=pk offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=sk offset=slave bundle=gmem2
#pragma HLS INTERFACE s_axilite port=return bundle=control

    // Local variables
    byte_t buf[64];
    byte_t rho[32], sigma[32];
    polyvec_t s_hat, e_hat, pkpv;
    poly_t s_hato;
    matrix_t A;
    
#pragma HLS ARRAY_PARTITION variable=buf complete
#pragma HLS ARRAY_PARTITION variable=rho complete
#pragma HLS ARRAY_PARTITION variable=sigma complete

    // Step 1: (rho, sigma) := G(d)
    G(d, 32, buf);

    // Split the 64-byte output
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        rho[i] = buf[i];
        sigma[i] = buf[i + 32];
    }

    // Step 2: Generate matrix A from rho
    matrix_expand(&A, rho);



    //print_matrix(A);    // Step 3: Generate secret vector s from sigma
    byte_t prf_buf_s[MLKEM_K * 64 * MLKEM_ETA1];
    byte_t prf_buf_s1[64 * MLKEM_ETA1];
    byte_t prf_buf_s2[64 * MLKEM_ETA1];
#pragma HLS ARRAY_PARTITION variable=prf_buf_s complete
    
    // Generate PRF output for each polynomial in s
    
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        prf_eta(MLKEM_ETA1, sigma, (byte_t)i, prf_buf_s + i * 64 * MLKEM_ETA1);
    }
    prf_eta(MLKEM_ETA1, sigma, (byte_t)0, prf_buf_s1 );
    //prf_eta(MLKEM_ETA1, sigma, (byte_t)1, prf_buf_s2 );

    polyvec_cbd_eta1(&s_hat, prf_buf_s);
    
    // Step 4: Generate error vector e from sigma
    byte_t prf_buf_e[MLKEM_K * 64 * MLKEM_ETA1];
#pragma HLS ARRAY_PARTITION variable=prf_buf_e complete
    
    // Generate PRF output for each polynomial in e
    for (int i = 0; i < MLKEM_K; i++) {
#pragma HLS UNROLL
        prf_eta(MLKEM_ETA1, sigma, (byte_t)(i + MLKEM_K), prf_buf_e + i * 64 * MLKEM_ETA1);
    }
    
    // Sample e using CBD
    polyvec_cbd_eta1(&e_hat, prf_buf_e);
    
    // Step 5: Transform s and e to NTT domain
    
    polyvec_ntt(&s_hat);

    polyvec_ntt(&e_hat);

    // Step 6: Compute t = A * s + e
    matrix_vector_mul(&pkpv, &A, &s_hat);

    polyvec_add(&pkpv, &pkpv, &e_hat);
    
    // Step 7: Reduce t modulo q
    polyvec_reduce(&pkpv);
    

    polyvec_tobytes(pk, &pkpv);
    //print_polyvec(pkpv);
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        pk[MLKEM_K * MLKEM_POLYBYTES + i] = rho[i];
    }
    
 
    polyvec_tobytes(sk, &s_hat);
    
    // Copy public key to secret key
    for (int i = 0; i < MLKEM_PUBLICKEYBYTES; i++) {
#pragma HLS PIPELINE II=1
        sk[MLKEM_K * MLKEM_POLYBYTES + i] = pk[i];
    }
    
    // Compute H(pk) and store in secret key
    byte_t pk_hash[32];
#pragma HLS ARRAY_PARTITION variable=pk_hash complete
    H(pk, MLKEM_PUBLICKEYBYTES, pk_hash);
    
    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        sk[MLKEM_K * MLKEM_POLYBYTES + MLKEM_PUBLICKEYBYTES + i] = pk_hash[i];
    }
    

    for (int i = 0; i < 32; i++) {
#pragma HLS UNROLL
        sk[MLKEM_K * MLKEM_POLYBYTES + MLKEM_PUBLICKEYBYTES + 32 + i] = z[i];
    }
}


