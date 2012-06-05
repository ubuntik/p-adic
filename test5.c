#include "p-adic.h"

float wrapped_indicator(pa_num *pa)
{
	pa_num *n;
	int gamma = -1;
	float ret;

	n = init_pa_num(-2, -2);
	set_x_by_gamma(n, -2, 1);
	ret = (float)indicator(pa, n, gamma);
	free_pa_num(n);

	return ret;
}

float p_norma_wrapped(pa_num* pa) { return (float) p_norma(pa); }

float function(pa_num* pa)
{
	float ret;
	//ret = wrapped_indicator(pa) * p_norma_wrapped(pa);
	//ret = wrapped_indicator(pa);
	ret = p_norma_wrapped(pa) * p_norma_wrapped(pa);
	return 1 / ret;
}

int main()
{
	float res;
	float res1, res2, res3;
	float (*pfunc)(pa_num* pnum);
	int i, gmax, gmin, gamma, fs_sz;
	pa_num **fs;
	complex ret;

	printf("Test#5: Integrals\n");

#if 0
	pfunc = function;
	gmin = -1;
	gmax = 0;
	res1 = integral(pfunc, gmin, gmax);
	gmin = -1;
	gmax = 1;
	res2 = integral(pfunc, gmin, gmax);
	gmin = -1;
	gmax = 2;
	res3 = integral(pfunc, gmin, gmax);

	printf("Zp: %f\nB1: %f\nB2: %f\n", res1, res2, res3);

	printf("\n>>> Indicator function <<<\n");
	pfunc = wrapped_indicator;

	for (gmax = 0; gmax < 4; gmax ++){
		for (gmin = 0; gmin > -4; gmin--){
			res = integral(pfunc, gmin, gmax);
			printf("gmin = %d gmax = %d result = %f\n", \
							gmin, gmax, res);
		}
	}

	printf("\n>>> p-norma function <<<\n");
	pfunc = p_norma_wrapped;

	gmax = 1;
	for (gmin = 0; gmin > -10; gmin--){
		res = integral(pfunc, gmin, gmax);
		printf("gmin = %d gmax = %d result = %f\n", gmin, gmax, res);
	}
	gmin = -1;
	for (gmax = 0; gmax < 10; gmax ++){
		res = integral(pfunc, gmin, gmax);
		printf("gmin = %d gmax = %d result = %f\n", gmin, gmax, res);
	}
#endif

	printf("\n>>> wavelet basis <<<\n");
	gmin = 0;
	gmax = 4;
	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);
	//for (gamma = gmin; gamma <= gmax; gamma++) {
		gamma = 0;
		printf("Current gamma: %d\n", gamma);
		fs_sz = fspace_sz(gamma, gmax);
		fs = gen_factor_space(gamma, gmax);
		printf("The power of factor space: %d\n", fs_sz);
		for (i = 0; i < fs_sz; i++) {
			//ret = wavelet_integral(fs[i], gamma, 1, gamma - 1, gmax);
			ret = wavelet_integral(fs[i], gamma, 1, gmin, gmax);
			printf(">>> Result: %f + i * %f\n\n", crealf(ret), cimagf(ret));
			free_pa_num(fs[i]);
		}
		free(fs);
	//}
	return 0;
}

