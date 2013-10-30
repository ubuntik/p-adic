#include "../src/cauchy_problem.h"
#include <errno.h>

// old values 0 and -5
#define G_MAX (0)
#define G_MIN (-5)

static const int gmin = G_MIN;
static const int gmax = G_MAX;

// function should be got fron func_args
// only for ALPHA = 2
double function(pa_num* pa)
{
	double ret = 0;
	if (pa == NULL) {
		fprintf(stderr, "Involid pointer\n");
		return -1;
	}
	ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= power(P, -gmax)) ? \
			power(P, 3*gmax) : ret;
			// p, - - (alfa + 1) * gmax
}

double wrapped_indicator(pa_num *pa)
{
	pa_num *n = NULL;
	int gamma = 0;
	double ret = 0;
	PADIC_ERR err = ESUCCESS;

	n = (pa_num *)malloc(sizeof(pa_num));
	if (n == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(n, 0, 0);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	ret = (double)indicator(pa, n, gamma);
	free_pa_num(n);

	return ret;
}

int main()
{
	PADIC_ERR err = ESUCCESS;

	err = solve_problem(function, wrapped_indicator, gmin, gmax);
	if (err != ESUCCESS)
		printf("FAILED\n");
	else
		printf("OK\n");

	return 0;
}

