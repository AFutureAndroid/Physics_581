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

void chebyshev_ten(int m, double *grid, double *values, double *first, double *second){
    double step = 2.0 / (m-1);
    for (int i=0; i<m; i++) {
        grid[i] = -1.0 + i*step;
        values[i] = cos(10*acos(grid[i]));
        first[i] = (100*i) - (1600*i*i*i) + (6720*pow(i,5)) - (10240*pow(i,7)) + (5120*pow(i,9));
        second[i] = 100 - 4800*i*i + 33600*pow(i,4) - 71680*pow(i,6) + 46080*pow(i,8);
    }
}


int main() {
    int m = 10; // Number of points
    double grid[m], values[m], first[m], second[m];

    // Compute the Chebyshev polynomial values on the grid
    chebyshev_ten(m, grid, values, first, second);

    // Print the grid points and corresponding T10(x) values
    printf("x\t\tT10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\n", grid[i], values[i]);
    }

    printf("x\t\tT'10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\n", grid[i], first[i]);
    }

    printf("x\t\tT''10(x)\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\n", grid[i], second[i]);
    }
    return 0;
}


