#include "../src/p-adic.h"

double function(pa_num* pa)
{
	double ret;
	printf("norma pa = %.16f\n", p_norm(pa));
	printf("1 / p = %.16f\n", (1.0 / P));
	printf("diff: %.16f\n", (p_norm(pa) - (1.0 / P)));
	ret = 1.0 / (p_norm(pa) * p_norm(pa));
//	return (p_norm(pa) <= (1.0 / P)) ? 0.0 : ret;
	return ((1.0 / P) >= p_norm(pa)) ? 0.0 : ret;
}

int main()
{
	double (*pfunc)(pa_num* pnum);
	int gmax, gmin, gamma;
	pa_num *n;
	complex ret;

	printf("Test#6: Wavelet integrals\n");

	pfunc = function;
	gmin = -5;
	gmax = 0;
	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);
	gamma = 0;
	printf("Current gamma: %d\n", gamma);
	n = init_pa_num(gmin, gmax);
	ret = wavelet_integral(pfunc, n, gamma, 1, gmin, gmax);
	printf(">>>>> Result: %f + i * %f\n\n", crealf(ret), cimagf(ret));
	free_pa_num(n);
	return 0;
}

