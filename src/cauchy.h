#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>

#include <p-analysis.h>

PADIC_ERR solve_problem(
                double (*rho)(pa_num *pnum),
                double (*start_cond)(pa_num *pnum),
                int gmin, int gmax, int gchy,
		int ini_gamma, pa_num *ini_n);
