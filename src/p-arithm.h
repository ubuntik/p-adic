/* library for p-adic arithmetic
 * at ultrametric space
 */

#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>

#include <p-def.h>

int arith_compare(pa_num *pa1, pa_num *pa2);

int reverse_sign(int sign);

PADIC_ERR sub(pa_num *res, pa_num *pa1, pa_num *pa2);

PADIC_ERR add(pa_num *res, pa_num *pa1, pa_num *pa2);

PADIC_ERR get_fractional_part(pa_num *fnum, pa_num *pa);

PADIC_ERR jmult(pa_num *res, pa_num *pa1, int j);

