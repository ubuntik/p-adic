/* library for p-adic arithmetic
 * at ultrametric space
 */

#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>

#include <p-arithm.h>

double p_norm(pa_num *pa);

int indicator(pa_num *x, pa_num *n, int gamma);

complex double character(pa_num *pa);

complex double wavelet(pa_num *x, pa_num *n, int gamma, int j);

double integral(double (*func)(pa_num* pnum), int g_min, int g_max);

complex double wavelet_integral(double (*func)(pa_num *pnum), pa_num *n, int gamma,
		int j, int g_min, int g_max);

