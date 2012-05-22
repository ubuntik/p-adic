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
	int gmax, gmin;

	printf("Test#5: Integrals\n");

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

	return 0;
}

