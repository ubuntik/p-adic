#include "../src/p-adic.h"

int gmin = -5;
int gmax = 4;

double function(pa_num* pa)
{
	double ret;
	ret = 1.f / (p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= (1.f / P)) ? 0.f : ret;
}

// only for alfa = 2
double op_vlad(pa_num *pa)
{
	double ret = -1;
	ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa));
//	return (p_norm(pa) <= power(P, -gmax)) ? \
//			power(P, 3 * gmax) : ret;
	return (p_norm(pa) <= power(P, -gmax)) ? \
			0.0 : ret;
}

double d_x(pa_num *pa)
{
	return 1;
}

int main()
{
	double (*pfunc)(pa_num* pnum) = NULL;
	int gamma;
	pa_num *n = NULL;
	complex ret;
	double num;
	PADIC_ERR err = ESUCCESS;

	n = (pa_num *)malloc(sizeof(pa_num));
	if (n == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}

	printf("Test#6: Wavelet integrals\n");

	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);

	pfunc = op_vlad;
	num = integral(pfunc, gmin, gmax);
	printf(">>>>> Vladimirov's operator: %g\n\n", num);

	pfunc = d_x;
	num = integral(pfunc, gmin, gmax);
	printf(">>>>> Integral dx: %g\n\n", num);


	pfunc = function;
	gamma = 0;
	printf("Current gamma: %d\n", gamma);
	err = init_pa_num(n, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init of number\n");
		exit(err);
	}
	ret = wavelet_integral(pfunc, n, gamma, 1, gmin, gmax);
	printf(">>>>> Result2: %f + i * %f\n\n", crealf(ret), cimagf(ret));
	free_pa_num(n);
	return 0;
}

