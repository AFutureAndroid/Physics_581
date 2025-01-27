// tutorial for c from tutorialspoint.com

#include <stdio.h>
#define PI 3.14159
#define AREA(r) (PI*r*r)

int main() {
    int radius = 5;
    float area = AREA(radius);

    printf("Area: %f", area);
    return 0;
}




