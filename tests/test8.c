#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cauchy.h>

#define G_MIN (-5)
#define G_CHY (0)
#define G_MAX (1)

#define ALPHA 2
#define BETA 0

static const int gmin = G_MIN;
static const int gmax = G_MAX;
static const int gchy = G_CHY;

//global variables of initial ball
pa_num *ini_n = NULL;
int ini_gamma = 0;

// function should be got fron func_args

double delta_func(pa_num *x, pa_num *y)
{
	double ret = 0;
	double norm_x = 0;
	double norm_y = 0;
	int gamma_cut = gchy;

	if (x == NULL || y == NULL) {
		fprintf(stderr, "Involid pointer\n");
		return -1;
	}

	norm_x = p_norm(x);
	norm_y = p_norm(y);

	if (norm_x == 0 && norm_y == 0)
		return 1;
	if (norm_y == 0)
		norm_y = -gamma_cut;
	if (norm_x == 0)
		norm_x = -gamma_cut;

	ret = norm_x / norm_y;

	return power(ret, -BETA);
}

double function(pa_num *x, pa_num *y)
{
	double ret = 0;
	int gamma_cut = gchy;
	PADIC_ERR err = ESUCCESS;
	pa_num *pa = NULL;

	if (x == NULL || y == NULL) {
		fprintf(stderr, "Involid pointer\n");
		return -1;
	}

	pa = (pa_num *)malloc(sizeof(pa_num));
	if (pa == NULL) {
		fprintf(stderr, "Not enough memory");
		return -1;
	}

	err = sub(pa, x, y);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed rho_forward\n");
		return -1;
	}

	if (ALPHA == 1) {
		ret = 1.0 / (p_norm(pa) * p_norm(pa));
		return (p_norm(pa) <= power(P, -gamma_cut)) ? \
			power(P, 2 * gamma_cut) : ret;
	} else if (ALPHA == 1.5) {
		ret = 1.0 / (p_norm(pa) * p_norm(pa) * sqrt(p_norm(pa)));
		return (p_norm(pa) <= power(P, -gamma_cut)) ? \
			power(P, 2.5 * gamma_cut) : ret;
	} else if (ALPHA == 2) {
		ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa));
		return (p_norm(pa) <= power(P, -gamma_cut)) ? \
			power(P, 3 * gamma_cut) : ret;
	} else if (ALPHA == 3) {
		ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa) * p_norm(pa));
		return (p_norm(pa) <= power(P, -gamma_cut)) ? \
			power(P, 4 * gamma_cut) : ret;
	} else if (ALPHA == 4) {
		ret = 1.0 / (p_norm(pa) * p_norm(pa) * p_norm(pa) * p_norm(pa) * p_norm(pa));
		return (p_norm(pa) <= power(P, -gamma_cut)) ? \
			power(P, 5 * gamma_cut) : ret;
	} else {
		fprintf(stderr, "UNDEFINED Alpha!\n");
		return 1;
	}
	// p, - (alfa + 1) * (-gmax)
}

double rho_bw(pa_num *x, pa_num *y)
{
	double norm_x = 0;
	double norm_y = 0;

	norm_x = p_norm(x);
	norm_y = p_norm(y);

	return function(x, y) * ((norm_y >= norm_x) ? 1 : delta_func(x, y));
}

double rho_fw(pa_num *x, pa_num *y)
{
	double norm_x = 0;
	double norm_y = 0;

	norm_x = p_norm(x);
	norm_y = p_norm(y);

/*
	fprintf(stderr, ">>> norm x = %g, y = %g\n", norm_x, norm_y);
	fprintf(stderr, ">>> function = %g\n", function(x, y));
	fprintf(stderr, ">>> delta_func = %g\n", delta_func(y, x));
*/

	return function(x, y) * ((norm_x >= norm_y) ? 1 : delta_func(y, x));
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

	err = solve_problem(rho_bw, rho_fw, wrapped_indicator,
//	err = solve_problem(function, function, wrapped_indicator,
				gmin, gmax, gchy);
	if (err != ESUCCESS)
		printf("FAILED\n");
	else
		printf("OK\n");

	return 0;
}

