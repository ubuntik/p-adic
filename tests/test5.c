#include "../src/p-adic.h"

double wrapped_indicator(pa_num *pa)
{
	pa_num *n = NULL;
	int gamma = -1;
	double ret = 0;
	PADIC_ERR err = ESUCCESS;

	n = (pa_num *)malloc(sizeof(pa_num));
	if (n == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(n, -2, -2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(n, -2, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
	}
	ret = (double)indicator(pa, n, gamma);
	free_pa_num(n);

	return ret;
}

double function(pa_num* pa)
{
	double ret;
	ret = wrapped_indicator(pa) * p_norm(pa);
	return ret;
}

double sp_point(pa_num* pa)
{
	double ret;
	ret = p_norm(pa) * p_norm(pa);
	return 1 / ret;
}

int main()
{
	double res;
	double res1, res2, res3;
	double (*pfunc)(pa_num* pnum);
	int gmax, gmin;

	printf("Test#5: Integrals\n");

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

	printf("\n>>> p-norm function <<<\n");
	pfunc = p_norm;

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

	printf("\n>>> integral with special point <<<\n");
	pfunc = sp_point;

	gmin = -4;
	gmax = 0;
	res1 = integral(pfunc, gmin, gmax);

	printf("Zp: %f\n", res1);

	return 0;
}

