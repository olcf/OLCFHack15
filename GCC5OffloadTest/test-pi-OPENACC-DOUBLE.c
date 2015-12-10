#include <stdio.h>

#define N 2000000000
/* Maximum parallelism allowed by Nvidia*/
#define vl 1024
double calcpi(long long n);

int main(void)
{

    double pi = 0.0;
    pi = calcpi(N);

    printf("pi=%11.10f\n", pi / N);
    return 0;

}

double calcpi(long long n)
{
    double pi = 0.0;
#pragma acc parallel vector_length(vl)
#pragma acc loop reduction(+:pi)
    for (long long ii = 0.0; ii < n; ii = ii + 1) {
	float t = (ii + 0.5f) / n;
	float s = 4.0f / (1.0f + t * t);
	pi = pi + s;
    }
    return pi;
}				//end calcpi
