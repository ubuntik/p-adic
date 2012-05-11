#include "p-adic.h"
#define G_MAX 1
#define G_MIN (-2)


int main()
{
	pa_num *pa;
	pa_num *x, *n;
	int gamma = -1;
	int ind;
	complex res;

	printf("Test#3: Wavelet functions\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	pa = init_pa_num(G_MIN, G_MAX);
	pa->x[0] = 1;
	pa->x[1] = 2;
	pa->x[2] = 1;
	pa->x[3] = 1;

	print_pa_num(pa);
	printf("Number: %g\n", from_canonic_to_float(pa));

	res = character(pa);
	printf("Character: %g + i%g\n", creal(res), cimag(res));

	x = init_pa_num(-1,0);
	x->x[0] = 1;
	x->x[1] = 1;
	print_pa_num(x);
	printf("x >> %g\n", from_canonic_to_float(x));

	n = init_pa_num(-1, -1);
	n->x[0] = 1;
	print_pa_num(n);
	printf("n >> %g\n", from_canonic_to_float(n));

	gamma = 0;

	ind = indicator(x, n, gamma);
	printf("gamma = %d, indicator = %d\n\n", gamma, ind);

	res = wavelet(x, n, gamma, 1);
	printf(">>%g + i%g\n", creal(res), cimag(res));

	return 0;
}

