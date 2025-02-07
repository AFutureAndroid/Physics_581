#include <stdio.h>

void compute_weights(int n, double x[], double w[]) {
    for (int j = 0; j < n; j++) {
        w[j] = 1.0;
        for (int k = 0; k < n; k++) {
            if (k != j) {
                w[j] *= 1.0 / (x[j] - x[k]);
            }
        }
    }
}