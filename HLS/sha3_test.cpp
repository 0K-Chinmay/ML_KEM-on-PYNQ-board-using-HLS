#include <stdio.h>
#include <stdint.h>
#include "unified.h"
// Alias for byte
//typedef uint8_t byte_t;
//typedef uint64_t lane_t;

// Declare SHA3-512 function (defined elsewhere, like in your HLS project)
//void sha3_512(const byte_t* input, int input_len, byte_t output[64]);

// Helper to convert hex string to byte array
void hexstr_to_bytes(const char* hex, byte_t* bytes, int len) {
    for (int i = 0; i < len; i++) {
        sscanf(hex + 2*i, "%2hhx", &bytes[i]);
    }
}

void print_bytes(const byte_t* bytes, int len) {
    for (int i = 0; i < len; i++) {
        printf("%02x", bytes[i]);
    }
    printf("\n");
}

int main12() {
    // Input as hex string
    const char* hex_input = "e1e3206875e67d7e81353774fe9025035b9b41a4a9f6ec00b91c600442fd717d";
    const int input_len = 32;

    // Convert to byte array
    byte_t input_bytes[32];
    hexstr_to_bytes(hex_input, input_bytes, input_len);

    // Output buffer for SHA3-512 result
    byte_t output[64];

    // Call SHA3-512
    sha3_512(input_bytes, input_len, output);

    // Print the 512-bit (64-byte) hash output
    printf("SHA3-512 Output:\n");
    print_bytes(output, 64);

    return 0;
}
