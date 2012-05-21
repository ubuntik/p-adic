#include "p-adic.h"

float wrapped_indicator(pa_num *pa)
{
	pa_num *n;
	int gamma = 0;
	float ret;

	n = init_pa_num(0, 0);
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

//	pfunc = wrapped_indicator;
	pfunc = p_norma_wrapped;

	for (gmax = 0; gmax < 2; gmax ++){
		for (gmin = 0; gmin > -1; gmin--){
			res = integral(pfunc, gmin, gmax);
			printf(">>> gmin = %d; gmax = %d; result = %f\n\n", \
							gmin, gmax, res);
		}
	}

	return 0;
}

