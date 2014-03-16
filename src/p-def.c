#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <strings.h>
#include <string.h>

#include <p-def.h>

// len = g_max - g_min + 1

double power(double base, double exponent)
{
	double half_pow = 0.0, int_part = 0.0;
	double fr_part = modf(exponent, &int_part);

	if ((fr_part != 0) && (fr_part != 0.5)) {
		return pow(base, exponent);
	} else if (fr_part == 0.5) {
		return (sqrt(base) * power(base, int_part));
	}

	if (exponent == 0)
		return (double)1;
	else if (exponent < 0)
		return (double)1 / power(base, -exponent);
	else if (fmod(exponent, 2) == 0) {
		half_pow = power(base, exponent / (double)2);
		return half_pow * half_pow;
	} else
		return base * power(base, exponent - 1);
}

long qspace_sz(int g_min, int g_max)
{
	long ret = 1;
	int i = 0;
	int xsz = g_max - g_min;

	for (i = 0; i < xsz; i++)
		ret = ret * P;
	return ret;
}

PADIC_ERR init_pa_num(pa_num* pa, int g_min, int g_max)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}
	pa->g_min = g_min;
	pa->g_max = g_max;
	pa->sign = POS;
	pa->x = (int *)calloc(pa->g_max - pa->g_min + 1, sizeof(int));
	if (pa->x == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	return ESUCCESS;
}

void free_pa_num(pa_num* pa)
{
	if (pa == NULL)
		fprintf(stderr, "Invalid pointer\n");
	free((void *)pa->x);
	free((void *)pa);
}

int get_x_by_gamma(pa_num* pa, int gamma)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}
	if (gamma < pa->g_min)
		return 0;
	else if (gamma > pa->g_max)
		return 0;
	return pa->x[gamma - pa->g_min];
}

PADIC_ERR set_x_by_gamma(pa_num* pa, int gamma, int x)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (gamma < pa->g_min || gamma > pa->g_max) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", pa->g_min);
		fprintf(stderr, "gamma should be less than %d\n", pa->g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}
	if (x < 0 || x >= P) {
		fprintf(stderr, "Invalid value of coefficient\n");
		return EINVCOEFF;
	}
	pa->x[gamma - pa->g_min] = x;
	return ESUCCESS;
}

PADIC_ERR gen_quotient_space(pa_num **qs, int g_min, int g_max)
{
	size_t qs_sz = 0;
	int div = 0, n = 0, i = 0, j = 0, k = 0, l = 0, ngrp = 0;
	PADIC_ERR err = ESUCCESS;
	int min = g_min - g_max - 1;
	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	qs_sz = (size_t)qspace_sz(g_min, g_max);

	for (i = 0; i < qs_sz; i++) {
		qs[i] = (pa_num *)malloc(sizeof(pa_num));
		if (qs[i] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
		err = init_pa_num(qs[i], min, -1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid init of number\n");
			return err;
		}
	}

	div = 1;
	for (k = -1; k >= min; k--) {
		n = 0;
		ngrp = qs_sz / (div * P) ;
		for (j = 0; j < ngrp; j++) {
			for (i = 0; i < P; i++) {
				for (l = 0; l < div; l++) {
					err = set_x_by_gamma(qs[n++], k, i);
					if (err != ESUCCESS) {
						fprintf(stderr,
							"Involid setting\n");
						return err;
					}
				}
			}
		}
		div = div * P;
	}

	return ESUCCESS;
}

// TODO: # 3.. 6:+: 1, 2, 1, 0
// according to flag: PA_LONG_PRINT

PADIC_ERR p_gamma_pa_num(pa_num *res, pa_num *pa, int gamma)
{
	int i = 0;
	PADIC_ERR err = ESUCCESS;

	if (res == NULL || pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	err = init_pa_num(res, pa->g_min + gamma, pa->g_max + gamma);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}
	res->sign = pa->sign;

	for (i = res->g_min; i <= res->g_max; i++) {
		err = set_x_by_gamma(res, i, get_x_by_gamma(pa, i - gamma));
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}
	}
	return ESUCCESS;
}


void print_pa_num(pa_num *pa)
{
	int i = 0;
	fprintf(stdout, "# %d..%d:%s: ",
				pa->g_min, pa->g_max, pa->sign ? "-" : "+");
	for (i = pa->g_min; i <= pa->g_max; i++)
		fprintf(stdout, "%d ", get_x_by_gamma(pa, i));
	fprintf(stdout, "\n");
	fflush(stdout);
}

double padic2double(pa_num *pa)
{
	int i = 0;
	double ret = 0, ppow = 0;

	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}

	ppow = power(P, pa->g_min);
	for (i = pa->g_min; i <= pa->g_max; i++) {
		ret += (double)get_x_by_gamma(pa, i) * ppow;
		ppow = ppow * P;
	}

	if (pa->sign == NEG)
		return -ret;

	return ret;
}

PADIC_ERR get_pa_tree(pa_tree *tree, int g_min, int g_max)
{
	pa_num *pgnum = NULL;
	int i = 0, qs_sz = 0, gamma = 0;
	pa_num **qs = NULL;
	PADIC_ERR err = ESUCCESS;
	int cnt = 0;

	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	tree->tree_sz = (power(P, g_max - g_min + 1) - 1) / (P - 1);

	tree->pa_nodes = (pa_node **)malloc(tree->tree_sz * sizeof(pa_node *));
	if (tree->pa_nodes == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}

	for (i = 0; i < tree->tree_sz; i++) {
		tree->pa_nodes[i] = (pa_node *)malloc(sizeof(pa_node));
		if (tree->pa_nodes[i] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
	}

	// first ball
	tree->pa_nodes[0]->data = 0.f;
	tree->pa_nodes[0]->parent = NULL;

	cnt = 1;

	for (gamma = g_min + 1; gamma <= g_max; gamma++) {
		qs_sz = qspace_sz(g_min, gamma);
		qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
		if (qs == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
		err = gen_quotient_space(qs, g_min, gamma);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid generating qspace\n");
			return err;
		}

		for (i = 0; i < qs_sz; i++) {
			pgnum = (pa_num *)malloc(sizeof(pa_num));
			if (pgnum == NULL) {
				fprintf(stderr, "Cannot alloc memory\n");
				return EINVPNTR;
			}
			err = p_gamma_pa_num(pgnum, qs[i], gamma);
			if (err != ESUCCESS) {
				fprintf(stderr, "Involid mult on pgamma\n");
				return err;
			}

			tree->pa_nodes[cnt]->data =
				padic2double(pgnum);
			tree->pa_nodes[cnt]->parent =
				tree->pa_nodes[(cnt - 1) / P];
			cnt++;
		}
	}
	return ESUCCESS;
}

void free_tree(pa_tree *tree)
{
	int i = 0;

	for (i = 0; i < tree->tree_sz; i++) {
		free(tree->pa_nodes[i]);
	}
	free(tree);
}

PADIC_ERR print_tree(pa_tree *tree, char* file_name)
{
	char header[] = "digraph \"p-adic_tree\" {\nn000;\nn000 [label=\"0\"];\n";
	unsigned int i = 0;
	char* my_str;
	int fd;

	my_str = (char *)malloc(64);

	if ((fd = open(file_name, O_WRONLY | O_CREAT, 0666)) < 0) {
		fprintf(stderr, "Can't create/open file %s\n", file_name);
		return EINVPNTR;
	}

	if (write(fd, (void *)header, strlen(header)) < 0) {
		fprintf(stderr, "Can't write to file %s\n", file_name);
		return EINVPNTR;
	}

	for (i = 1; i < tree->tree_sz; i++) {
		bzero((void *)my_str, 64);
		snprintf(my_str, 42, "n%03u -> n%03u;\nn%03u [label=\"%3.7f\"];\n",
				(i - 1) / P, i, i, tree->pa_nodes[i]->data);
		if (write(fd, (void *)my_str, strlen(my_str)) < 0) {
			fprintf(stderr, "Can't write to file %s\n", file_name);
			return EINVPNTR;
		}
	}
	bzero((void *)my_str, 64);
	snprintf(my_str, 2, "}\n");
	if (write(fd, (void *)my_str, strlen(my_str)) < 0) {
		fprintf(stderr, "Can't write to file %s\n", file_name);
		return EINVPNTR;
	}

	close(fd);
	return ESUCCESS;
}

