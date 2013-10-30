/* library for p-adic arithmetic
 * at ultrametric space
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <fcntl.h>

#define _Complex double

#define P 3
#define PI 3.14159265358979323846f

//Workaround:
#ifndef NAN
#define NAN (0.0/0.0)
#endif

//Workaround:
#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif


typedef enum {ESUCCESS = 0,
	EGAMMAOUT = 1,
	EINVPNTR = 2,
	EINVCOEFF = 3,
	EINVJ = 4,
	EMALLOC = 5
	} PADIC_ERR;

enum PADIC_SIGN {POS = 0,
		NEG = 1
		};

struct pa_num {
	int g_min;
	int g_max;
	int err;
	int sign;
	int *x;
};

typedef struct pa_num pa_num;

typedef struct pa_node pa_node;

struct pa_node {
	double data;
	pa_node *parent;
};

struct pa_tree {
	int tree_sz;
	pa_node **pa_nodes;
};

typedef struct pa_tree pa_tree;

double power(double base, double exponent);

PADIC_ERR init_pa_num(pa_num *pa, int gmin, int gmax);

void free_pa_num(pa_num *pa);

long qspace_sz(int g_min, int g_max);

int get_x_by_gamma(pa_num* pa, int gamma);

PADIC_ERR set_x_by_gamma(pa_num* pa, int gamma, int x);

PADIC_ERR __extend_number(pa_num *ext_pa, pa_num *pa, int g_min, int g_max);

PADIC_ERR __sign_sub(pa_num *res, pa_num *pn1, pa_num *pn2);

int arith_compare(pa_num *pa1, pa_num *pa2);

int reverse_sign(int sign);

PADIC_ERR sub(pa_num *res, pa_num *pa1, pa_num *pa2);

PADIC_ERR add(pa_num *res, pa_num *pa1, pa_num *pa2);

PADIC_ERR gen_quotient_space(pa_num **qs, int g_min, int g_max);

PADIC_ERR p_gamma_pa_num(pa_num *res, pa_num *pa, int gamma);

void print_pa_num(pa_num *pa);

double p_norm(pa_num *pa);

int indicator(pa_num *x, pa_num *n, int gamma);

double from_canonic_to_double(pa_num *pa);

PADIC_ERR get_fractional_part(pa_num *fnum, pa_num *pa);

complex character(pa_num *pa);

complex wavelet(pa_num *x, pa_num *n, int gamma, int j);

PADIC_ERR jmult(pa_num *res, pa_num *pa1, int j);

double integral(double (*func)(pa_num* pnum), int g_min, int g_max);

complex wavelet_integral(double (*func)(pa_num *pnum), pa_num *n, int gamma,
		int j, int g_min, int g_max);

double integral_B_x(double (*func)(pa_num *pnum), pa_num *x, int g_min, int g_max);

complex wavelet_integral_C_gnj_x(double (*func)(pa_num *pnum), pa_num *x,
		pa_num *n, int gamma, int j, int g_min, int g_max);

complex wavelet_integral_Agnj(double (*func)(pa_num *pnum), pa_num *n, int gamma,
		int j, int g_min, int g_max);

PADIC_ERR get_pa_tree(pa_tree *tree, int g_min, int g_max);

void free_tree(pa_tree *tree);

PADIC_ERR print_tree(pa_tree *tree, char* file_name);

