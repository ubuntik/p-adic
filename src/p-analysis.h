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

complex character(pa_num *pa);

complex wavelet(pa_num *x, pa_num *n, int gamma, int j);

double integral(double (*func)(pa_num* pnum), int g_min, int g_max);

complex wavelet_integral(double (*func)(pa_num *pnum), pa_num *n, int gamma,
		int j, int g_min, int g_max);

double integral_B_x(double (*func)(pa_num *pnum), pa_num *x, int g_min, int g_max);

complex wavelet_integral_C_gnj_x(double (*func)(pa_num *pnum), pa_num *x,
		pa_num *n, int gamma, int j, int g_min, int g_max);

complex wavelet_integral_Agnj(double (*func)(pa_num *pnum), pa_num *n, int gamma,
		int j, int g_min, int g_max);

