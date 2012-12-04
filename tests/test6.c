#include "../src/p-adic.h"

double p_norm_wrapped(pa_num* pa) { return (double) p_norm(pa); }

double function(pa_num* pa)
{
	double ret;
	ret = p_norm_wrapped(pa) * p_norm_wrapped(pa);
	return 1 / ret;
}

int main()
{
	double (*pfunc)(pa_num* pnum);
	int gmax, gmin, gamma;
	pa_num *n;
	complex ret;

	printf("Test#6: Wavelet integrals\n");

	pfunc = function;
	gmin = -4;
	gmax = 0;
	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);
	gamma = 0;
	printf("Current gamma: %d\n", gamma);
	n = init_pa_num(-4, 0);
	ret = wavelet_integral(pfunc, n, gamma, 0, gmin, gmax);
	printf(">>>>> Result: %f + i * %f\n\n", crealf(ret), cimagf(ret));
	free_pa_num(n);
	return 0;
}

