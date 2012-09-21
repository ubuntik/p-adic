#include "../src/p-adic.h"
#define G_MAX 1
#define G_MIN (-2)


int main()
{
	pa_num *pa;
	pa_num *x, *n;
	int gamma = -1;
	int ind;
	int i, j;
	complex res;

	printf("Test#3: Wavelet functions\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	pa = init_pa_num(G_MIN, G_MAX);
	set_x_by_gamma(pa, G_MIN, 1);
	set_x_by_gamma(pa, G_MIN + 1, 2);
	set_x_by_gamma(pa, G_MAX - 1, 1);
	set_x_by_gamma(pa, G_MAX, 1);

	print_pa_num(pa);
	printf("Number: %g\n", from_canonic_to_float(pa));

	res = character(pa);
	printf("Character: %g + i%g\n", creal(res), cimag(res));
	free_pa_num(pa);

	for (i = 0; i < P; i++) {
		x = init_pa_num(-1,0);
		set_x_by_gamma(x, -1, 1);
		set_x_by_gamma(x, 0, i);
		print_pa_num(x);
		printf("x >> %g\n", from_canonic_to_float(x));

		n = init_pa_num(-1, -1);
		set_x_by_gamma(n, -1, 1);
		print_pa_num(n);
		printf("n >> %g\n", from_canonic_to_float(n));

		gamma = 0;

		ind = indicator(x, n, gamma);
		printf("gamma = %d, indicator = %d\n", gamma, ind);

		for (j = 1; j < P; j++) {
			res = wavelet(x, n, gamma, j);
			printf("(j = %d) >> %g + i%g\n", j, creal(res), cimag(res));
		}
		printf("\n");
		free_pa_num(x);
		free_pa_num(n);
	}

	return 0;
}

