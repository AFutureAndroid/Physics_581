/* Physics 581 Assignment 3
*/

// Question 2 Part a

/*
Letting j = 1, k = 2

then H_j = 2x and H_k = -2+4x^2
*/

#include <stdio.h>
#include <math.h>

#define EPSILON 1e-16
#define STEP 0.01

double H1(double x) {
    return 2*x;
}

double H2(double x) {
    return -2+4*x*x;
}

int main() {
    double start = 0.0001;
    double integrand = 1.0;
    while (fabs(integrand) > EPSILON) {
        integrand = exp(-start*start)*H1(start)*H2(start);
        printf("Value chack: %.18f",integrand);
        start += 0.00001;
    }
    printf("x_max = %.4f\n", start);
    return 0;
}

