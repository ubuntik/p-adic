/* library for p-adic arithmetic
 * at ultrametric space
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <math.h>
#include <complex.h>

#define P 2
#define PI 3.14159265358979323846f

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

long fspace_sz(int g_min, int g_max);

int get_x_by_gamma(pa_num* pa, int gamma);

int set_x_by_gamma(pa_num* pa, int gamma, int x);

pa_num* __extend_number(pa_num *pa, int g_min, int g_max);

pa_num* minus(pa_num *pa1, pa_num *pa2);

pa_num **gen_factor_space(int g_min, int g_max);

pa_num* p_gamma_pa_num(pa_num *pa, int gamma);

void print_pa_num(pa_num *pa);

int p_norma(pa_num *pa);

int indicator(pa_num *x, pa_num *n, int gamma);

float from_canonic_to_float(pa_num *pa);

pa_num* get_fractional_part(pa_num *pa);

complex character(pa_num *pa);

complex wavelet(pa_num *x, pa_num *n, int gamma, int j);

pa_num* smult(pa_num *pa1, int j);

float integral(float (*func)(pa_num* pnum), int g_min, int g_max);

/* Not implemented yet */
pa_num* mult(pa_num *pa1, pa_num *pa2);

/*TODO*/

