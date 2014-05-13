#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cauchy.h>

//#define G_MIN (-7)
//#define G_CHY (0)
//#define G_MAX (2)

#define G_MIN (-2)
#define G_CHY (0)
#define G_MAX (1)

static const int gmin = G_MIN;
static const int gmax = G_MAX;
static const int gchy = G_CHY;

// function should be got fron func_args
double function(pa_num* pa)
{
	double ret = 0;
	if (pa == NULL) {
		fprintf(stderr, "Involid pointer\n");
		return -1;
	}

/* Alpha = 1 */
/*
	ret = 1.0 / (p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= power(P, -gchy)) ? \
			power(P, 2 * gchy) : ret;
*/
/* Alpha = 2 */

	ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= power(P, -gchy)) ? \
			power(P, 3 * gchy) : ret;

/* Alpha = 3 */
/*
	ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= power(P, -gchy)) ? \
			power(P, 4 * gchy) : ret;
*/
/* Alpha = 4 */
/*
	ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa) * p_norm(pa) * p_norm(pa));
	return (p_norm(pa) <= power(P, -gchy)) ? \
			power(P, 5 * gchy) : ret;
*/
			// p, - (alfa + 1) * (-gmax)
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

	int ini_gamma = 0;
	pa_num *ini_n = NULL;

	if (gmax <= gchy) {
		printf("gmax should be more then gchy\n");
		return -1;
	}

	ini_n = (pa_num *)malloc(sizeof(pa_num));
	if (ini_n == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return -1;
	}

	err = init_pa_num(ini_n, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}

	err = solve_problem(function, wrapped_indicator, gmin, gmax, gchy, ini_gamma, ini_n);
	if (err != ESUCCESS)
		printf("FAILED\n");
	else
		printf("OK\n");

	return 0;
}

