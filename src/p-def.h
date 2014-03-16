/* library for p-adic arithmetic
 * at ultrametric space
 */
#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>

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


typedef enum {
	ESUCCESS = 0,
	EGAMMAOUT = 1,
	EINVPNTR = 2,
	EINVCOEFF = 3,
	EINVJ = 4,
	EMALLOC = 5
} PADIC_ERR;

enum PADIC_SIGN {
	POS = 0,
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

struct pa_node {
	double data;
	struct pa_node *parent;
};

typedef struct pa_node pa_node;

struct pa_tree {
	int tree_sz;
	pa_node **pa_nodes;
};

typedef struct pa_tree pa_tree;

/* basic functions */

double power(double base, double exponent);

PADIC_ERR init_pa_num(pa_num *pa, int gmin, int gmax);

void free_pa_num(pa_num *pa);

int get_x_by_gamma(pa_num* pa, int gamma);

PADIC_ERR set_x_by_gamma(pa_num* pa, int gamma, int x);

long qspace_sz(int g_min, int g_max);

PADIC_ERR gen_quotient_space(pa_num **qs, int g_min, int g_max);

PADIC_ERR p_gamma_pa_num(pa_num *res, pa_num *pa, int gamma);

void print_pa_num(pa_num *pa);

double padic2double(pa_num *pa);

PADIC_ERR get_pa_tree(pa_tree *tree, int g_min, int g_max);

void free_tree(pa_tree *tree);

PADIC_ERR print_tree(pa_tree *tree, char* file_name);

