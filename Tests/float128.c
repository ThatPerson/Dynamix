#include <stdio.h>
#include <math.h>
#include <quadmath.h>
int main(void) {
	__float128 a;
	char buf[128];
	a  =  2.0q;
	a  =  sqrtq(a);
	int n  =  quadmath_snprintf(buf, sizeof(buf), "%.36Qe", a);
	printf("%s\n", buf);
	a  =  a * 0.000000001q;
	n  =  quadmath_snprintf(buf, sizeof(buf), "%.36Qe", a);
	printf("%s\n", buf);
	return 1;
}
