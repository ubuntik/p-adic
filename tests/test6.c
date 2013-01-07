#include "../src/p-adic.h"

double function(pa_num* pa)
{
	double ret;
	ret = 1.0 / (p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= (1.0 / P)) ? 0.0 : ret;
}

int main()
{
	double (*pfunc)(pa_num* pnum);
	int gmax, gmin, gamma;
	pa_num *n = NULL;
	complex ret;
	PADIC_ERR err = ESUCCESS;

	n = (pa_num *)malloc(sizeof(pa_num));
	if (n == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}

	printf("Test#6: Wavelet integrals\n");

	pfunc = function;
	gmin = 0;
	gmax = 3;
	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);
	gamma = 0;
	printf("Current gamma: %d\n", gamma);
	err = init_pa_num(n, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init of number\n");
		exit(err);
	}
	ret = wavelet_integral(pfunc, n, gamma, 1, gmin, gmax);
	printf(">>>>> Result: %f + i * %f\n\n", crealf(ret), cimagf(ret));
	free_pa_num(n);
	return 0;
}

