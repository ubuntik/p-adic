#include "p-adic.h"

float wrapped_indicator(pa_num *pa)
{
	pa_num *n;
	int gamma = 0;
	float ret;

	n = init_pa_num(-1, -1);
	set_x_by_gamma(n, -1, 1);
	ret = (float)indicator(pa, n, gamma);
	free_pa_num(n);

	return ret;
}

float p_norma_wrapped(pa_num* pa) { return (float) p_norma(pa); }

int main()
{
	float res;
	float (*pfunc)(pa_num* pnum);
	int i, gmax, gmin, gamma, fs_sz;
	pa_num **fs;
	complex ret;

	printf("Test#5: Integrals\n");

#if 0
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
	for (gmin = 0; gmin > -15; gmin--){
		res = integral(pfunc, gmin, gmax);
		printf("gmin = %d gmax = %d result = %f\n", gmin, gmax, res);
	}
	gmin = -1;
	for (gmax = 0; gmax < 15; gmax ++){
		res = integral(pfunc, gmin, gmax);
		printf("gmin = %d gmax = %d result = %f\n", gmin, gmax, res);
	}
#endif


	printf("\n>>> wavelet basis <<<\n");
	gamma = 0;
	gmax = 1;
	fs_sz = fspace_sz(gamma, gmax);
	fs = gen_factor_space(gamma, gmax);
	printf("fs sz = %d\n", fs_sz);
	for (i = 0; i < fs_sz; i++) {
		ret = wavelet_integral(fs[i], fs[i], gamma, 1, gamma - 1, gmax);
		printf(">>>> %f + i * %f\n", crealf(ret), cimagf(ret));

	}

	return 0;
}

