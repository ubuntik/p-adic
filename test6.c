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
	ret = p_norma_wrapped(pa) * p_norma_wrapped(pa);
	return 1 / ret;
}

int main()
{
	float (*pfunc)(pa_num* pnum);
	int i, gmax, gmin, gamma, fs_sz;
	pa_num **fs;
	complex ret;

	printf("Test#6: Wavelet integrals\n");

	pfunc = function;
	gmin = -2;
	gmax = 2;
	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);
	//for (gamma = gmin; gamma <= gmax; gamma++) {
		gamma = 0;
		printf("Current gamma: %d\n", gamma);
		fs_sz = fspace_sz(gamma, gmax);
		fs = gen_factor_space(gamma, gmax);
		printf("The power of factor space: %d\n", fs_sz);
		for (i = 0; i < fs_sz; i++) {
			//ret = wavelet_integral(fs[i], gamma, 1, gamma - 1, gmax);
			ret = wavelet_integral(pfunc, fs[i], gamma, 1, gmin, gmax);
			printf(">>> Result: %f + i * %f\n\n", crealf(ret), cimagf(ret));
			free_pa_num(fs[i]);
		}
		free(fs);
	//}
	return 0;
}

