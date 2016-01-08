
#include <stdlib.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    uint64_t value = atol(argv[1]), b;
    int64_t res;
    //return __builtin_cpu_supports_popcount();
    return (int)__builtin_clzll(value);
    //__asm__ volatile ("lzcnt %1, %0" : "=r"(res) : "r"(value));
    //__asm__ volatile ("popcnt %1, %0" : "=r"(b) : "r"(a));
    //return res;
}
