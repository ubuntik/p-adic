#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cauchy.h>

#define G_MIN (-2)
#define G_CHY (0)
#define G_MAX (1)

static const int gmin = G_MIN;
static const int gmax = G_MAX;
static const int gchy = G_CHY;

//global variables of initial ball
pa_num *ini_n = NULL;
int ini_gamma = 0;

static PADIC_ERR set_ini_n()
{
	PADIC_ERR err = ESUCCESS;

	ini_n = (pa_num *)malloc(sizeof(pa_num));
	if (ini_n == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	err = init_pa_num(ini_n, 0, 0);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		return err;
	}

// MAY be set_x_by_gamma(...)

	return err;
}

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
	double ret = 0;
/*
	PADIC_ERR err = ESUCCESS;

	err = set_ini_n();
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init <n> number\n");
		return err;
	}
*/

	ret = (double)indicator(pa, ini_n, ini_gamma);
	return ret;
}

int main()
{
	PADIC_ERR err = ESUCCESS;

	if (gmax <= gchy) {
		printf("gmax should be more then gchy\n");
		return -1;
	}

	err = set_ini_n();
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init <n> number\n");
		return err;
	}

	err = count_S_t(function, wrapped_indicator, gmin, gmax, gchy, ini_n, ini_gamma);
	if (err != ESUCCESS)
		printf("FAILED\n");
	else
		printf("OK\n");

	return 0;
}

