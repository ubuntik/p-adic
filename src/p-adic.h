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


typedef enum {ESUCCESS = 0,
	EGAMMAOUT = 1,
	EINVPNTR = 2,
	EINVCOEFF = 3,
	EINVJ = 4
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

float power(float base, int exponent);

PADIC_ERR init_pa_num(pa_num *pa, int gmin, int gmax);

void free_pa_num(pa_num *pa);

long qspace_sz(int g_min, int g_max);

int get_x_by_gamma(pa_num* pa, int gamma);

PADIC_ERR set_x_by_gamma(pa_num* pa, int gamma, int x);

int arith_compare(pa_num *pa1, pa_num *pa2);

int reverse_sign(int sign);

PADIC_ERR sub(pa_num *res, pa_num *pa1, pa_num *pa2);

PADIC_ERR add(pa_num *res, pa_num *pa1, pa_num *pa2);

PADIC_ERR gen_quotient_space(pa_num **qs, int g_min, int g_max);

PADIC_ERR p_gamma_pa_num(pa_num *res, pa_num *pa, int gamma);

void print_pa_num(pa_num *pa);

float p_norm(pa_num *pa);

int indicator(pa_num *x, pa_num *n, int gamma);

float from_canonic_to_float(pa_num *pa);

PADIC_ERR get_fractional_part(pa_num *fnum, pa_num *pa);

complex character(pa_num *pa);

complex wavelet(pa_num *x, pa_num *n, int gamma, int j);

PADIC_ERR jmult(pa_num *res, pa_num *pa1, int j);

float integral(float (*func)(pa_num* pnum), int g_min, int g_max);

complex wavelet_integral(float (*func)(pa_num *pnum), pa_num *n, int gamma,
		int j, int g_min, int g_max);

float integral_B_x(float (*func)(pa_num *pnum), pa_num *x, int g_min, int g_max);

complex wavelet_integral_C_gnj_x(float (*func)(pa_num *pnum), pa_num *x,
		pa_num *n, int gamma, int j, int g_min, int g_max);

