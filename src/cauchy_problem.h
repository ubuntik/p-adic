#include "../src/p-adic.h"

PADIC_ERR solve_problem(
                double (*rho)(pa_num *pnum),
                double (*start_cond)(pa_num *pnum),
                int gmin, int gmax);
