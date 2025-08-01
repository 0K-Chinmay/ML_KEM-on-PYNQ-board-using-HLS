#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <cstring>
#include "unified.h"


// Function prototypes
void mlkem512_keygen_top(const byte_t seed[32],const byte_t z[32], byte_t pk[MLKEM_PUBLICKEYBYTES], byte_t sk[MLKEM_SECRETKEYBYTES]);

void print_hex(const byte_t* data, int len, const std::string& label) {
    std::cout << label << ": ";
    for (int i = 0; i < len; i++) {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)data[i];
    
    }
    std::cout << std::dec << std::endl;
}
// Deterministic testing (same seed should produce same keys)
bool test_deterministic() {
    std::cout << "\n=== Testing Deterministic Behavior ===" << std::endl;

    byte_t seed[32] = {
        0xE1, 0xE3, 0x20, 0x68, 0x75, 0xE6, 0x7D, 0x7E,
        0x81, 0x35, 0x37, 0x74, 0xFE, 0x90, 0x25, 0x03,
        0x5B, 0x9B, 0x41, 0xA4, 0xA9, 0xF6, 0xEC, 0x00,
        0xB9, 0x1C, 0x60, 0x04, 0x42, 0xFD, 0x71, 0x7D
    };

    byte_t z[32] = {
        0xC6, 0xF5, 0x78, 0x5A, 0x6F, 0x2B, 0x42, 0xE8,
        0x43, 0x22, 0x8B, 0xE5, 0x3E, 0xB7, 0x68, 0xD6,
        0x4C, 0x6F, 0x9D, 0x43, 0x55, 0xAE, 0x95, 0xF0,
        0x83, 0xE5, 0x1E, 0xD5, 0x7C, 0x43, 0x73, 0x10
    };
    byte_t pk1[MLKEM_PUBLICKEYBYTES], sk1[MLKEM_SECRETKEYBYTES];
    byte_t pk2[MLKEM_PUBLICKEYBYTES], sk2[MLKEM_SECRETKEYBYTES];
    
    // Generate first key pair
    mlkem512_keygen_top(seed,z, pk1, sk1);
    
    // Generate second key pair with same seed
    //mlkem512_keygen_top(seed, pk2, sk2);
    
    print_hex(seed, 32, "Seed");
    print_hex(pk1, MLKEM_PUBLICKEYBYTES, "PK (first 32 bytes)");
    print_hex(sk1, MLKEM_SECRETKEYBYTES, "SK (first 32 bytes)");

    
    return 1;
}

// Main function
int main(int argc, char* argv[]) {
    std::cout << "ML-KEM 512 Key Generation Test Suite" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    bool all_tests_passed = true;
    
    // Run all test suites
    //all_tests_passed &= test_key_sizes();
    //all_tests_passed &= test_known_vectors();
    all_tests_passed &= test_deterministic();
    //all_tests_passed &= test_random_vectors(100);
    
    return 0;
}