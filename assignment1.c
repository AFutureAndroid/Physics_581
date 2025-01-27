// Question 3 
//We start by storing the n=10 Chebyshev polynomial

#include <stdio.h>

int main(void){
    int x = 1;
    int T = -1 + 50*x*x - 400*x*x*x*x + 1120*x*x*x*x*x*x - 1280*x*x*x*x*x*x*x*x + 512*x*x*x*x*x*x*x*x*x*x;
    int xvals[] = {-1, -0.5, 0, 0.5, 1};

    printf("%d\n", xvals[2]);
    return 0;

}
