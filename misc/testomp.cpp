#include <omp.h>
#include <stdio.h>

int main(const int argc, const char *const *const argv) {
    int result = 0;
    #pragma omp parallel for
    for (int i = 0; i < 1000*1000*1000; i++) {
        if (i%(1000*1000) == 0) { printf("%d\n", i/(1000*1000)); }
        for (int j = 0; j < 1000; j++) {
            result = result * 3 + i * j;
        }
    }
    printf("Dummyresult=%d\n", result);
}

