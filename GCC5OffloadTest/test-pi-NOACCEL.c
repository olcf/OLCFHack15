#include <stdio.h>

#define N 2000000000
int main(void)
{

    double pi = 0.0;
    for (long long ii = 0.0; ii < n; ii = ii + 1) {
	double t = (ii + 0.5) / n;
	double s = 4.0 / (1.0 + t * t);
	pi = pi + s;
    }
    printf("pi=%11.10f\n", pi / N);
    return 0;
}
