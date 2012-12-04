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


enum PADIC_ERR {EINIT = 0,
		EGOUT = 1
		};

enum PADIC_SIGN {POS = 0,
		NEG = 1
		};

struct pa_num {
	int g_min;
	int g_max;
	int *x;
	int err;
	int sign;
};

typedef struct pa_num pa_num;

pa_num* init_pa_num(int gmin, int gmax);

void free_pa_num(pa_num *pa);

long qspace_sz(int g_min, int g_max);

int get_x_by_gamma(pa_num* pa, int gamma);

int set_x_by_gamma(pa_num* pa, int gamma, int x);

int arith_compare(pa_num *pa1, pa_num *pa2);

int reverse_sign(int sign);

pa_num* sub(pa_num *pa1, pa_num *pa2);

pa_num* add(pa_num *pa1, pa_num *pa2);

pa_num **gen_quotient_space(int g_min, int g_max);

pa_num* p_gamma_pa_num(pa_num *pa, int gamma);

void print_pa_num(pa_num *pa);

double p_norm(pa_num *pa);

int indicator(pa_num *x, pa_num *n, int gamma);

double from_canonic_to_float(pa_num *pa);

pa_num* get_fractional_part(pa_num *pa);

complex character(pa_num *pa);

complex wavelet(pa_num *x, pa_num *n, int gamma, int j);

pa_num* jmult(pa_num *pa1, int j);

double integral(double (*func)(pa_num* pnum), int g_min, int g_max);

complex wavelet_integral(double (*func)(pa_num *pnum), pa_num *n, int gamma, int j, int g_min, int g_max);

/* Not implemented yet */
pa_num* mult(pa_num *pa1, pa_num *pa2);

/*TODO*/

