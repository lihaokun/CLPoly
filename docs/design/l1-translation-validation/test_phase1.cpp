#include <cstdint>
#include <vector>
#include <stdexcept>

// Test 1: while loop
uint64_t find_first_nonzero(const uint64_t* arr, int n) {
    int i = 0;
    while (i < n) {
        if (arr[i] != 0) break;
        i++;
    }
    return i;
}

// Test 2: throw
uint64_t safe_div(uint64_t a, uint64_t b) {
    if (b == 0) throw std::runtime_error("div by zero");
    return a / b;
}

// Test 3: for + continue
uint64_t sum_odd(const uint64_t* arr, int n) {
    uint64_t sum = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] % 2 == 0) continue;
        sum += arr[i];
    }
    return sum;
}
