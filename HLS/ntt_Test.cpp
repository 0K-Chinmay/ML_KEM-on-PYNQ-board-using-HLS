#include <iostream>
#include <stdint.h>
#include "ap_int.h"

#include "unified.h"

// Main test
int main2() {
    poly_t p;

    // Initialize polynomial with simple data
    for (int i = 0; i < MLKEM_N; i++) {
        p.coeffs[i] = i+1 % 3329;
    }

    std::cout << "Input Polynomial:\n";
    //print_poly(p);

    // Run NTT
    ntt_forward(&p);

    std::cout << "Output Polynomial:\n";
    //print_poly(p);

    return 0;
}
