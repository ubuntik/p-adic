#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cauchy.h>

#define G_MIN (0)
#define G_CHY (1)
#define G_MAX (2)

#define ALPHA 2

static const int gmin = G_MIN;
static const int gmax = G_MAX;
static const int gchy = G_CHY;

//global variables of initial ball
pa_num *ini_n = NULL;
int ini_gamma = 0;

// function should be got fron func_args

double rho_asterisk(pa_num *z)
{
	return exp(-2 * p_norm(z));
}

double rho_a(pa_num *x, pa_num *y)
{
	PADIC_ERR err = ESUCCESS;
	pa_num *z = NULL;

	z = (pa_num *)malloc(sizeof(pa_num));
	if (z == NULL) {
		fprintf(stderr, "Not enough memory");
		return -1;
	}

	err = sub(z, x, y);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed rho_backward\n");
		return -1;
	}

	return exp(-2*p_norm(z));
}


double rho_bw(pa_num *x, pa_num *y)
{
	PADIC_ERR err = ESUCCESS;
	pa_num *z = NULL;

	z = (pa_num *)malloc(sizeof(pa_num));
	if (z == NULL) {
		fprintf(stderr, "Not enough memory");
		return -1;
	}

	err = sub(z, x, y);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed rho_backward\n");
		return -1;
	}

	//return rho_asterisk(z);
	return exp(p_norm(y)) * rho_asterisk(z);
}

double rho_fw(pa_num *x, pa_num *y)
{
	PADIC_ERR err = ESUCCESS;
	pa_num *z = NULL;

	z = (pa_num *)malloc(sizeof(pa_num));
	if (z == NULL) {
		fprintf(stderr, "Not enough memory");
		return -1;
	}

	err = sub(z, x, y);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed rho_forward\n");
		return -1;
	}

	return exp(p_norm(x)) * rho_asterisk(z);
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

//	err = solve_problem(rho_bw, rho_fw, wrapped_indicator,
	err = solve_problem(function, function, wrapped_indicator,
				gmin, gmax, gchy);
	if (err != ESUCCESS)
		printf("FAILED\n");
	else
		printf("OK\n");

	return 0;
}

