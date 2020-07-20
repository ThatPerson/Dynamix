#include <stdio.h>
#include <math.h>

#define double long double

int main(void) {
	double a;
	a  =  0.000000002;
	a  =  a * a;
	printf("%Le\n",a);
	return 1;
}
