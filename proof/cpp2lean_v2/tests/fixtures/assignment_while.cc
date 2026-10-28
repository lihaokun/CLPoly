#include <cstdint>

uint64_t assignment_while(uint64_t a, uint64_t b) {
    uint64_t c;
    while ((c = a % b)) {
        a = b;
        b = c;
    }
    return b;
}
