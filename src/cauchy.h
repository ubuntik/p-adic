#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>

#include <p-analysis.h>

PADIC_ERR solve_problem(
                double (*rho_bw)(pa_num *x, pa_num *y),
                double (*rho_fw)(pa_num *x, pa_num *y),
                double (*start_cond)(pa_num *pnum),
                int gmin, int gmax, int gchy,
		pa_num *x0);


PADIC_ERR count_S_t(
                double (*rho_bw)(pa_num *x, pa_num *y),
                double (*rho_fw)(pa_num *x, pa_num *y),
                double (*start_cond)(pa_num *pnum),
                int gmin, int gmax, int gchy,
		pa_num *ini_n, int ini_gamma);

