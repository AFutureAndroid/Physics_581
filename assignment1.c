// Question 3 
//We start by storing the n=10 Chebyshev polynomial

/** 
#include <stdio.h>

int main(void){
    int x = 1;
    int T = -1 + 50*x*x - 400*x*x*x*x + 1120*x*x*x*x*x*x - 1280*x*x*x*x*x*x*x*x + 512*x*x*x*x*x*x*x*x*x*x;
    int xvals[] = {-1, -0.5, 0, 0.5, 1};

    printf("%d\n", xvals[2]);
    return 0;
}
**/


#include <stdio.h>
#include <math.h>
#include <string.h>

void chebyshev_ten(int m, double *grid, double *values, double *first, double *second, double step){
    for (int i=0; i<m; i++) {
        grid[i] = -1.0 + i*step;
        double j = grid[i];
        values[i] = cos(10*acos(grid[i]));
        first[i] = (100*j) - (1600*j*j*j) + (6720*pow(j,5)) - (10240*pow(j,7)) + (5120*pow(j,9));
        second[i] = 100 - 4800*j*j + 33600*pow(j,4) - 71680*pow(j,6) + 46080*pow(j,8);
    }
}

void fd_method(int m, double *grid, double *values, double *first_fdm, double *second_fdm, double step){
    // First derivative
    first_fdm[0] = (values[1] - values[0]) / step;
    first_fdm[m-1] = (values[m - 1] - values[m - 2]) / step;
    for (int i=1; i < m - 1; i++) {
        first_fdm[i] = (values[i + 1] - values[i - 1]) / (2 * step);
    }

    // Second derivative
    second_fdm[0] = (values[2] + values[0] - 2*values[1]) / (step*step);
    second_fdm[m-1] = (values[m-1] + values[m-3] - 2*values[m-2]) / (step*step);
    for(int i=1; i < m-1; i++) {
        second_fdm[i] = (values[i+1] + values[i-1] - 2*values[i]) / (step*step);
    }
}


void errors(int m, double *first, double *second, double *first_fdm, double *second_fdm, double first_err, double second_err){
    first_err = 0.0;
    second_err = 0.0;

    for (int i = 0; i < m; i++) {
        first_err += fabs(first[i] - first_fdm[i]);
        printf("%f|", first[i]);
        printf("%f|", first_fdm[i]);
        printf("\n");
        second_err += fabs(second[i] - second_fdm[i]);
    }

    first_err = first_err / m;
    second_err = second_err / m;
}


int main() {
    int m = 10; // Number of points
    double grid[m], values[m], first[m], second[m];
    double step = 2.0 / (m-1);

    // Compute the Chebyshev polynomial values on the grid
    chebyshev_ten(m, grid, values, first, second, step);

    // Print the grid points and corresponding T10(x) values
    printf("x\t\t\tT10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\n", grid[i], values[i]);
    }

    printf("x\t\t\tT'10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\n", grid[i], first[i]);
    }

    printf("x\t\t\tT''10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\n", grid[i], second[i]);
    }

    // Approximate first and second derivatives from listed data
    double *first_fdm = first;
    double *second_fdm = second;
    fd_method(m, grid, values, first_fdm, second_fdm, step);
    printf("x\t\t\tT10(x)\t\tT'10(x)\t\tT''10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\t%f\t%f\n", grid[i], values[i], first_fdm[i], second_fdm[i]);
    }

    // Error calculation
    double first_err, second_err;
    errors(m, first, second, first_fdm, second_fdm, first_err, second_err);
    
    printf("\nTotal Error in T'10(x): %f\n", first_err);
    printf("Total Error in T''10(x): %f\n", second_err);

    return 0;
}



