#include "p-adic.h"

// len = g_max - g_min + 1

long qspace_sz(int g_min, int g_max)
{
	long ret = 1;
	int i = 0;
	int xsz = g_max - g_min;

	for (i = 0; i < xsz; i++)
		ret = ret * P;
	return ret;
}

PADIC_ERR init_pa_num(pa_num* pa, int g_min, int g_max)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}
	pa->g_min = g_min;
	pa->g_max = g_max;
	pa->sign = POS;
	pa->x = (int *)calloc(pa->g_max - pa->g_min + 1, sizeof(int));
	if (pa->x == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
//	bzero((void *)pa->x, sizeof(int) * (pa->g_max - pa->g_min + 1));
	return ESUCCESS;
}

void free_pa_num(pa_num* pa)
{
	if (pa == NULL)
		fprintf(stderr, "Invalid pointer\n");
	free((void *)pa->x);
	free((void *)pa);
}

int get_x_by_gamma(pa_num* pa, int gamma)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}
	if (gamma < pa->g_min)
		return 0;
	else if (gamma > pa->g_max)
		return 0;
	return pa->x[gamma - pa->g_min];
}

PADIC_ERR set_x_by_gamma(pa_num* pa, int gamma, int x)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (gamma < pa->g_min || gamma > pa->g_max) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", pa->g_min);
		fprintf(stderr, "gamma should be less than %d\n", pa->g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}
	if (x < 0 || x >= P) {
		fprintf(stderr, "Invalid value of coefficient\n");
		return EINVCOEFF;
	}
	pa->x[gamma - pa->g_min] = x;
	return ESUCCESS;
}

PADIC_ERR __extend_number(pa_num *ext_pa, pa_num *pa, int g_min, int g_max)
{
	PADIC_ERR err = ESUCCESS;
	int i = 0;
	int xtmp = 0;

	if (pa == NULL || ext_pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (g_min > pa->g_min) {
		fprintf(stderr, "Invalid Value g_min: %d\n", g_min);
		fprintf(stderr, "g_min should be less than %d\n", pa->g_min);
		fflush(stderr);
		return EGAMMAOUT;
	}
	if (g_max < pa->g_max) {
		fprintf(stderr, "Invalid Value g_max: %d\n", g_max);
		fprintf(stderr, "g_max should be greater than %d\n", pa->g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	err = init_pa_num(ext_pa, g_min, g_max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init of number\n");
		exit(err);
	}

	for (i = g_min; i <= g_max; i++) {
		xtmp = get_x_by_gamma(pa, i);
		if (xtmp < 0 || xtmp >= P) {
			fprintf(stderr, "Invalid value of coefficient\n");
			return EINVCOEFF;
		}
		err = set_x_by_gamma(ext_pa, i, xtmp);
		if (err != ESUCCESS) {
			fprintf(stderr, "Cannot set value of coefficient\n");
			return err;
		}
	}

	return ESUCCESS;
}

int arith_compare(pa_num *pa1, pa_num *pa2)
{
	int max =0 , min = 0, i = 0;
	if (pa1 == NULL || pa2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	min = (pa1->g_min < pa2->g_min) ? pa1->g_min : pa2->g_min;
	max = (pa1->g_max > pa2->g_max) ? pa1->g_max : pa2->g_max;

	for (i = max; i >= min; i--) {
		if (get_x_by_gamma(pa1, i) > get_x_by_gamma(pa2, i))
			return 1;
		else if (get_x_by_gamma(pa1, i) < get_x_by_gamma(pa2, i))
			return -1;
		else
			continue;
	}
	return 0;
}

int reverse_sign(int sign)
{
	return (sign == POS) ? NEG : POS;
}

int __s_get_x_by_gamma(pa_num* pa, int gamma)
{
	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}
	if (gamma < pa->g_min)
		return pa->sign ? (P - 1) : 0;
	else if (gamma > pa->g_max)
		return 0;
	return pa->x[gamma - pa->g_min];
}

PADIC_ERR __sign_sub(pa_num *res, pa_num *pn1, pa_num *pn2)
{
	int in_mind = 0;
	int coeff = 0, i = 0;
	PADIC_ERR err = ESUCCESS;

	if (res == NULL || pn1 == NULL || pn2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	err = init_pa_num(res, (pn1->g_min < pn2->g_min) ? pn1->g_min : pn2->g_min,
			(pn1->g_max > pn2->g_max) ? (pn1->g_max + 1) : (pn2->g_max + 1));
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}


	for (i = res->g_min; i <= res->g_max; i++) {
		(in_mind) ? (coeff = __s_get_x_by_gamma(pn1, i) - \
				__s_get_x_by_gamma(pn2, i) - 1) : \
			(coeff = __s_get_x_by_gamma(pn1, i) - \
				__s_get_x_by_gamma(pn2, i));
		if (coeff >= 0) {
			set_x_by_gamma(res, i, coeff);
			in_mind = 0;
		} else {
			set_x_by_gamma(res, i, coeff + P);
			in_mind = 1;
		}
	}
	return ESUCCESS;
}

PADIC_ERR __dummy_sub(pa_num *res, pa_num *pn1, pa_num *pn2)
{
	int in_mind = 0;
	int coeff = 0, i = 0;

	if (res == NULL || pn1 == NULL || pn2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	for (i = res->g_min; i <= res->g_max; i++) {
		(in_mind) ? (coeff = get_x_by_gamma(pn1, i) - \
				get_x_by_gamma(pn2, i) - 1) : \
			(coeff = get_x_by_gamma(pn1, i) - \
				get_x_by_gamma(pn2, i));
		if (coeff >= 0) {
			set_x_by_gamma(res, i, coeff);
			in_mind = 0;
		} else {
			set_x_by_gamma(res, i, coeff + P);
			in_mind = 1;
		}
	}
	return ESUCCESS;
}

PADIC_ERR __dummy_add(pa_num *res, pa_num *pn1, pa_num *pn2)
{
	int in_mind = 0;
	int coeff = 0, i = 0;

	if (res == NULL || pn1== NULL || pn2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	for (i = res->g_min; i <= res->g_max; i++) {
		(in_mind) ? (coeff = get_x_by_gamma(pn1, i) + \
				get_x_by_gamma(pn2, i) + 1) : \
			(coeff = get_x_by_gamma(pn1, i) + \
				get_x_by_gamma(pn2, i));
		if (coeff < P) {
			set_x_by_gamma(res, i, coeff);
			in_mind = 0;
		} else {
			set_x_by_gamma(res, i, coeff - P);
			in_mind = 1;
		}
	}
	return ESUCCESS;
}

PADIC_ERR __do_compact(pa_num *shrt, pa_num *pa)
{
	int i = 0, min = 0, max = 0;
	PADIC_ERR err = ESUCCESS;

	if (shrt == NULL || pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	for (i = pa->g_min; i <= pa->g_max; i++) {
		if (get_x_by_gamma(pa, i) != 0)
			break;
	}
	min = i;

	/* it means that pa is null */
	if (min >= pa->g_max) {
		err = init_pa_num(shrt, 0, 0);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid init of number\n");
			return err;
		}
	}
	for (i = pa->g_max; i >= pa->g_min; i--) {
		if (get_x_by_gamma(pa, i) != 0)
			break;
	}
	max = i;

	err = init_pa_num(shrt, min, max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}
	shrt->sign = pa->sign;
	memcpy((void *)shrt->x, (void *)&(pa->x[min - pa->g_min]), \
		(max - min + 1) * sizeof(int));

	return ESUCCESS;
}

PADIC_ERR add(pa_num *res, pa_num *pa1, pa_num *pa2)
{
	int min = 0, max = 0;
	int cmp = 0;
	pa_num *ext_pa = NULL;
	PADIC_ERR err = ESUCCESS;

	if (res == NULL || pa1 == NULL || pa2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	cmp = arith_compare(pa1, pa2);

	if ((cmp == 0) && (pa1->sign != pa2->sign)) {
		err = init_pa_num(res, 0, 0);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid init of number\n");
			return err;
		}
		return ESUCCESS;
	}

	min = (pa1->g_min < pa2->g_min) ? pa1->g_min : pa2->g_min;
	max = (pa1->g_max > pa2->g_max) ? pa1->g_max : pa2->g_max;

	ext_pa = (pa_num *)malloc(sizeof(pa_num));
	if (ext_pa == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}

	/* +1 for number's overhead */
	err = init_pa_num(ext_pa, min, max + 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}

	if (pa1->sign == pa2->sign) {
		res->sign = pa1->sign;
		err = __dummy_add(ext_pa, pa1, pa2);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid dummy add\n");
			return err;
		}
	} else {
		res->sign = (cmp > 0) ? pa1->sign : reverse_sign(pa1->sign);
		err = (cmp > 0) ? __dummy_sub(ext_pa, pa1, pa2) : \
						__dummy_sub(ext_pa, pa2, pa1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid dummy sub\n");
			return err;
		}
	}
	err = __do_compact(res, ext_pa);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid compacting number\n");
		return err;
	}
	free_pa_num(ext_pa);
	return ESUCCESS;
}

PADIC_ERR sub(pa_num *res, pa_num *pa1, pa_num *pa2)
{
	int min = 0, max = 0;
	int cmp = 0;
	pa_num *ext_pa = NULL;
	PADIC_ERR err = ESUCCESS;

	if (res == NULL || pa1 == NULL || pa2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	cmp = arith_compare(pa1, pa2);

	if ((cmp == 0) && (pa1->sign == pa2->sign)) {
		err = init_pa_num(res, 0, 0);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid init of number\n");
			return err;
		}
		return ESUCCESS;
	}

	min = (pa1->g_min < pa2->g_min) ? pa1->g_min : pa2->g_min;
	max = (pa1->g_max > pa2->g_max) ? pa1->g_max : pa2->g_max;

	ext_pa = (pa_num *)malloc(sizeof(pa_num));
	if (ext_pa == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}

	/* +1 for number's overhead */
	err = init_pa_num(ext_pa, min, max + 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}

	if (pa1->sign != pa2->sign) {
		ext_pa->sign = pa1->sign;
		err = __dummy_add(ext_pa, pa1, pa2);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid dummy add\n");
			return err;
		}
	} else {
		ext_pa->sign = (cmp > 0) ? pa1->sign : reverse_sign(pa1->sign);
		err = (cmp > 0) ? __dummy_sub(ext_pa, pa1, pa2) : \
						__dummy_sub(ext_pa, pa2, pa1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid dummy sub\n");
			return err;
		}
	}
	err = __do_compact(res, ext_pa);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid compacting number\n");
		return err;
	}
	free_pa_num(ext_pa);
	return ESUCCESS;
}

PADIC_ERR gen_quotient_space(pa_num **qs, int g_min, int g_max)
{
	size_t qs_sz = 0;
	int div = 0, n = 0, i = 0, j = 0, k = 0, l = 0, ngrp = 0;
	PADIC_ERR err = ESUCCESS;
	int min = g_min - g_max - 1;
	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	qs_sz = (size_t)qspace_sz(g_min, g_max);

	for (i = 0; i < qs_sz; i++) {
		qs[i] = (pa_num *)malloc(sizeof(pa_num));
		if (qs[i] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
		err = init_pa_num(qs[i], min, -1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid init of number\n");
			return err;
		}
	}

	div = 1;
	for (k = -1; k >= min; k--) {
		n = 0;
		ngrp = qs_sz / (div * P) ;
		for (j = 0; j < ngrp; j++) {
			for (i = 0; i < P; i++) {
				for (l = 0; l < div; l++) {
					err = set_x_by_gamma(qs[n++], k, i);
					if (err != ESUCCESS) {
						fprintf(stderr,
							"Involid setting\n");
						return err;
					}
				}
			}
		}
		div = div * P;
	}

	return ESUCCESS;
}

PADIC_ERR p_gamma_pa_num(pa_num *res, pa_num *pa, int gamma)
{
	int i = 0;
	PADIC_ERR err = ESUCCESS;

	if (res == NULL || pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	err = init_pa_num(res, pa->g_min + gamma, pa->g_max + gamma);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}
	res->sign = pa->sign;

	for (i = res->g_min; i <= res->g_max; i++) {
		err = set_x_by_gamma(res, i, get_x_by_gamma(pa, i - gamma));
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}
	}
	return ESUCCESS;
}

void print_pa_num(pa_num *pa)
{
	int i = 0;
	if (pa->sign == NEG)
		fprintf(stdout, "Negative number:\t");
	for (i = pa->g_min; i <= pa->g_max; i++) {
		fprintf(stdout, "x[%d] = %d\t", i, get_x_by_gamma(pa, i));
	}
	fprintf(stdout, "\n");
	fflush(stdout);
}

float power(float base, int exponent)
{
	float half_pow = 0;
	if (exponent == 0)
		return (float)1;
	else if (exponent < 0)
		return (float)1 / power(base, -exponent);
	else if (fmod(exponent, 2) == 0) {
		half_pow = power(base, exponent / (float)2);
		return half_pow * half_pow;
	}
	else
		return base * power(base, exponent - 1);
}

float p_norm(pa_num *pa)
{
	int i = 0, res = INT_MAX;

	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	for (i = pa->g_min; i <= pa->g_max; i++) {
		if (get_x_by_gamma(pa, i) == 0) {
			continue;
		} else {
			res = i;
			break;
		}
	}

	/* it means, we check all coeffs in an array and all of them are null*/
	if (res == INT_MAX)
		return (float)0;

	return power((float)P, -res);
}

int indicator(pa_num *x, pa_num *n, int gamma)
{
	pa_num *mult = NULL, *subtr = NULL;
	int norm = 0;
	PADIC_ERR err = ESUCCESS;

	if (x == NULL || n == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}

	mult = (pa_num *)malloc(sizeof(pa_num));
	if (mult == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return -1;
	}
	subtr = (pa_num *)malloc(sizeof(pa_num));
	if (subtr == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return -1;
	}

	err = p_gamma_pa_num(mult, x, gamma);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid multiplication on p-gamma\n");
		return -1;
	}
	err = sub(subtr, mult, n);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid subtraction\n");
		return -1;
	}
	norm = p_norm(subtr);

	/* |p^gamma*x - n| <= 1 => indicator = 1
	 * |kern| <= 1 => gamma >= 0
	 */

	free_pa_num(subtr);
	free_pa_num(mult);
	return (norm > 1) ? 0 : 1;
}

float from_canonic_to_float(pa_num *pa)
{
	int i = 0;
	float ret = 0;
	float ppow = 0;

	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}

	ppow = power(P, pa->g_min);
	for (i = pa->g_min; i <= pa->g_max; i++) {
		ret += (float)get_x_by_gamma(pa, i) * ppow;
		ppow = ppow * P;
	}

	if (pa->sign == NEG)
		return -ret;

	return ret;
}

PADIC_ERR get_fractional_part(pa_num *res, pa_num *pa)
{
	int i = 0;
	PADIC_ERR err = ESUCCESS;

	if (res == NULL || pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	if (pa->g_min >= 0) {
		err = init_pa_num(res, 0, 0);
		if (err != ESUCCESS)
			fprintf(stderr, "Invalid init of number\n");
		return err;
	}

	err = init_pa_num(res, pa->g_min, -1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}

	res->sign = pa->sign;

	for (i = pa->g_min; i < 0; i++) {
		err = set_x_by_gamma(res, i, get_x_by_gamma(pa, i));
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}
	}
	return ESUCCESS;
}

complex character(pa_num *pa)
{
	pa_num *fnum = NULL;
	complex ret = I;
	PADIC_ERR err = ESUCCESS;

	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return ret;
	}

	fnum = (pa_num *)malloc(sizeof(pa_num));
	if (fnum == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return ret;
	}

	err = get_fractional_part(fnum, pa);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid getting fractional part\n");
		return ret;
	}

	ret = cexpf(2 * PI * I *
			from_canonic_to_float(fnum));
	free_pa_num(fnum);
	return ret;
}

complex wavelet(pa_num *x, pa_num *n, int gamma, int j)
{
	pa_num *kern = NULL, *jkern = NULL;
	pa_num *mult = NULL, *subtr = NULL;
	complex ret = I;
	PADIC_ERR err = ESUCCESS;

	if (x == NULL || n == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return ret;
	}
	if (j <= 0 || j >= P) {
		fprintf(stderr, "Invalid value of j\n");
		return ret;
	}

	kern = (pa_num *)malloc(sizeof(pa_num));
	jkern = (pa_num *)malloc(sizeof(pa_num));
	mult = (pa_num *)malloc(sizeof(pa_num));
	subtr = (pa_num *)malloc(sizeof(pa_num));
	if (kern == NULL || jkern == NULL || mult == NULL || subtr == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return ret;
	}

	err = p_gamma_pa_num(mult, x, gamma);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid multiplication on p-gamma\n");
		return ret;
	}
	err = sub(subtr, mult, n);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid subtraction\n");
		return ret;
	}
	err = p_gamma_pa_num(kern, subtr, -1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid multiplication on p-gamma\n");
		return ret;
	}
	err = jmult(jkern, kern, j);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid multiplication on j\n");
		return ret;
	}

	if (indicator(x, n, gamma))
		ret = crealf(character(jkern)) + I * cimagf(character(jkern));
	else
		ret = 0 + I * 0;

	free_pa_num(jkern);
	free_pa_num(kern);
	free_pa_num(mult);
	free_pa_num(subtr);

	return ret;
}

/* workaround for wavelet: suppose (0 < j < P) and (pa > 0) */
PADIC_ERR jmult(pa_num *res, pa_num *pa, int j)
{
	PADIC_ERR err = ESUCCESS;
	int i = 0, tmp = 0, in_mind = 0;

	if (res == NULL || pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (j <= 0 || j >= P) {
		fprintf(stderr, "Invalid value of j\n");
		return EINVJ;
	}

	err = init_pa_num(res, pa->g_min, pa->g_max + 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid init of number\n");
		return err;
	}

	for (i = res->g_min; i <= res->g_max; i++) {
		tmp = get_x_by_gamma(pa, i) * j + in_mind;
		(tmp >= P) ? set_x_by_gamma(res, i, tmp - P) : \
			set_x_by_gamma(res, i, tmp);
		in_mind = (tmp >= P) ? 1 : 0;
	}

	return ESUCCESS;
}

#if 0
pa_num* mult(pa_num *pa1, pa_num *pa2)
{
	pa_num *ret;

	ret = init_pa_num(pa1->g_min + pa2->g_min, pa1->g_max + pa2->g_max);
	fprintf(stderr, "Not implemented yet\n");
	fflush(stderr);

	return 0;
}
#endif

float integral(float (*func)(pa_num *pnum), int g_min, int g_max)
{
	float ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL;
	float (*pfunc)(pa_num *pnum) = NULL;

	if (func == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	qs_sz = (size_t)qspace_sz(g_min, g_max);
	qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
	if (qs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	err = gen_quotient_space(qs, g_min, g_max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		return err;
	}

	pfunc = func;

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}

		if ( pfunc(pa) != INFINITY ) {
			ret += (float)pfunc(pa);
			free_pa_num(qs[i]);
			free_pa_num(pa);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%f is special point!\n", from_canonic_to_float(pa));
		spoint = (pa_num *)malloc(sizeof(pa_num));
		if (spoint == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = __extend_number(spoint, pa, pa->g_min,
							pa->g_max + 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid extending number\n");
			return err;
		}
		err = set_x_by_gamma(spoint, pa->g_max + 1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}
		printf("Take x = %f\n", from_canonic_to_float(spoint));
		print_pa_num(spoint);
		ret += (float)pfunc(spoint);
		free_pa_num(qs[i]);
		free_pa_num(pa);
		free_pa_num(spoint);
	}
	free(qs);
	return ret * power((float)P, -g_max);
}

complex wavelet_integral(float (*func)(pa_num *pnum), pa_num *n, int gamma, \
						int j, int g_min, int g_max)
{
	complex ret = I;
	PADIC_ERR err = ESUCCESS;
	float img = 0, rez = 0, fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL;
	float (*pfunc)(pa_num *pnum);

	if (n == NULL || func == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (j <= 0 || j >= P) {
		fprintf(stderr, "Invalid value of j\n");
		return EINVJ;
	}
	if (gamma < g_min || gamma > g_max) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", g_min);
		fprintf(stderr, "gamma should be less than %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	pfunc = func;
	/* for integration on the "level - 1" */
	qs_sz = (size_t)qspace_sz(g_min, g_max);
	qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
	if (qs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	err = gen_quotient_space(qs, g_min, g_max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		return err;
	}

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}
		/* workaround for special point (f.i. x = 0) */
		if ( pfunc(pa) != INFINITY ) {
			fun = pfunc(pa);
			rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
			img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
			free_pa_num(pa);
			free_pa_num(qs[i]);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%f is special point!\n", from_canonic_to_float(pa));

		spoint = (pa_num *)malloc(sizeof(pa_num));
		if (spoint == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = __extend_number(spoint, pa, pa->g_min,
							pa->g_max + 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid extending number\n");
			return err;
		}
		err = set_x_by_gamma(spoint, pa->g_max + 1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}

		printf("Take x = %f\n", from_canonic_to_float(spoint));
		print_pa_num(spoint);

		fun = pfunc(spoint);
		rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
		img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(qs[i]);
	}
	ret = rez * power((float)P, -g_max) + I * img * power((float)P, -g_max);
	free(qs);
	return ret;

}

float integral_B_x(float (*func)(pa_num *pnum), pa_num *x, int g_min, int g_max)
{
	float ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL, *res= NULL;
	float (*pfunc)(pa_num *pnum) = NULL;

	if (func == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}

	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	pfunc = func;

	qs_sz = (size_t)qspace_sz(g_min, g_max);
	qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
	if (qs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	err = gen_quotient_space(qs, g_min, g_max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		return err;
	}

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}

		res = (pa_num *)malloc(sizeof(pa_num));
		if (res == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = sub(res, pa, x);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid subtraction\n");
			return err;
		}

		if ( pfunc(res) != INFINITY ) {
			ret += (float)pfunc(res);
			free_pa_num(qs[i]);
			free_pa_num(pa);
			free_pa_num(res);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%f is special point!\n", from_canonic_to_float(pa));

		spoint = (pa_num *)malloc(sizeof(pa_num));
		if (spoint == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = __extend_number(spoint, pa, pa->g_min,
							pa->g_max + 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid extending number\n");
			return err;
		}
		err = set_x_by_gamma(spoint, pa->g_max + 1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}

		free_pa_num(res);
		res = (pa_num *)malloc(sizeof(pa_num));
		if (res == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = sub(res, spoint, x);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid subtraction\n");
			return err;
		}

		ret += (float)pfunc(res);
		free_pa_num(qs[i]);
		free_pa_num(pa);
		free_pa_num(spoint);
	}
	free(qs);
	return ret * power((float)P, -g_max);
}

complex wavelet_integral_C_gnj_x(float (*func)(pa_num *pnum), pa_num *x,
		pa_num *n, int gamma, int j, int g_min, int g_max)
{
	complex ret = I;
	PADIC_ERR err = ESUCCESS;
	float img = 0, rez = 0, fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL, *res = NULL;
	float (*pfunc)(pa_num *pnum);

	if (n == NULL || func == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (j <= 0 || j >= P) {
		fprintf(stderr, "Invalid value of j\n");
		return EINVJ;
	}
	if (gamma < g_min || gamma > g_max) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", g_min);
		fprintf(stderr, "gamma should be less than %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	pfunc = func;
	/* for integration on the "level - 1" */
	qs_sz = (size_t)qspace_sz(g_min, g_max);
	qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
	if (qs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	err = gen_quotient_space(qs, g_min, g_max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		return err;
	}

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}

		res = (pa_num *)malloc(sizeof(pa_num));
		if (res == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = __sign_sub(res, x, pa);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid subtraction\n");
			return err;
		}

		/* workaround for special point (f.i. x = 0) */
		if ( pfunc(res) != INFINITY ) {
			fun = pfunc(res);
			rez += crealf(wavelet(res, n, -gamma, j)) * fun;
			img += cimagf(wavelet(res, n, -gamma, j)) * fun;
			free_pa_num(pa);
			free_pa_num(res);
			free_pa_num(qs[i]);
			continue;
		}

		printf("Warning!!! Special point has found!\n");

		spoint = (pa_num *)malloc(sizeof(pa_num));
		if (spoint == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = __extend_number(spoint, res, res->g_min,
							res->g_max + 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid extending number\n");
			return err;
		}
		err = set_x_by_gamma(spoint, res->g_max + 1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}

		printf("Take x = %f\n", from_canonic_to_float(spoint));
		print_pa_num(spoint);

		fun = pfunc(spoint);
		rez += crealf(wavelet(res, n, -gamma, j)) * fun;
		img += cimagf(wavelet(res, n, -gamma, j)) * fun;
		free_pa_num(res);
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(qs[i]);
	}
	ret = rez * power((float)P, -g_max) + I * img * power((float)P, -g_max);
	free(qs);
	return ret;

}

complex wavelet_integral_Agnj(float (*func)(pa_num *pnum), pa_num *n, int gamma, \
						int j, int g_min, int g_max)
{
	// start condition + sopryazhenniy wavelet
	complex ret = I;
	PADIC_ERR err = ESUCCESS;
	float img = 0, rez = 0, fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL;
	float (*pfunc)(pa_num *pnum);

	if (n == NULL || func == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return EINVPNTR;
	}
	if (j <= 0 || j >= P) {
		fprintf(stderr, "Invalid value of j\n");
		return EINVJ;
	}
	if (gamma < g_min || gamma > g_max) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", g_min);
		fprintf(stderr, "gamma should be less than %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	pfunc = func;
	/* for integration on the "level - 1" */
	qs_sz = (size_t)qspace_sz(g_min, g_max);
	qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
	if (qs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	err = gen_quotient_space(qs, g_min, g_max);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		return err;
	}

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}
		/* workaround for special point (f.i. x = 0) */
		if ( pfunc(pa) != INFINITY ) {
			fun = pfunc(pa);
			rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
			img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
			free_pa_num(pa);
			free_pa_num(qs[i]);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%f is special point!\n", from_canonic_to_float(pa));

		spoint = (pa_num *)malloc(sizeof(pa_num));
		if (spoint == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = __extend_number(spoint, pa, pa->g_min,
							pa->g_max + 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid extending number\n");
			return err;
		}
		err = set_x_by_gamma(spoint, pa->g_max + 1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid setting coeff\n");
			return err;
		}

		printf("Take x = %f\n", from_canonic_to_float(spoint));
		print_pa_num(spoint);

		fun = pfunc(spoint);
		rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
		img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(qs[i]);
	}
	ret = rez * power((float)P, -g_max) + I * img * power((float)P, -g_max);
	free(qs);
	return ret;

}

PADIC_ERR get_pa_tree(pa_tree *tree, int g_min, int g_max)
{
	pa_num *pgnum = NULL;
	int i = 0, qs_sz = 0, gamma = 0;
	pa_num **qs = NULL;
	PADIC_ERR err = ESUCCESS;
	int cnt = 0;

	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	tree->tree_sz = (power(P, g_max - g_min + 1) - 1) / (P - 1);

	tree->pa_nodes = (pa_node **)malloc(tree->tree_sz * sizeof(pa_node *));
	if (tree->pa_nodes == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}

	for (i = 0; i < tree->tree_sz; i++) {
		tree->pa_nodes[i] = (pa_node *)malloc(sizeof(pa_node));
		if (tree->pa_nodes[i] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
	}

	// first ball
	tree->pa_nodes[0]->data = 0.f;
	tree->pa_nodes[0]->parent = NULL;

	cnt = 1;

	for (gamma = g_min + 1; gamma <= g_max; gamma++) {
		qs_sz = qspace_sz(g_min, gamma);
		qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
		if (qs == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
		err = gen_quotient_space(qs, g_min, gamma);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid generating qspace\n");
			return err;
		}

		for (i = 0; i < qs_sz; i++) {
			pgnum = (pa_num *)malloc(sizeof(pa_num));
			if (pgnum == NULL) {
				fprintf(stderr, "Cannot alloc memory\n");
				return EINVPNTR;
			}
			err = p_gamma_pa_num(pgnum, qs[i], gamma);
			if (err != ESUCCESS) {
				fprintf(stderr, "Involid mult on pgamma\n");
				return err;
			}

			tree->pa_nodes[cnt]->data =
				from_canonic_to_float(pgnum);
			tree->pa_nodes[cnt]->parent =
				tree->pa_nodes[(cnt - 1) / P];
			cnt++;
		}
	}
	return ESUCCESS;
}

void free_tree(pa_tree *tree)
{
	int i = 0;

	for (i = 0; i < tree->tree_sz; i++) {
		free(tree->pa_nodes[i]);
	}
	free(tree);
}

PADIC_ERR print_tree(pa_tree *tree, char* file_name)
{
	char header[] = "digraph \"p-adic_tree\" {\nn000;\nn000 [label=\"0\"];\n";
	unsigned int i = 0;
	char* my_str;
	int fd;

	my_str = (char *)malloc(64);

	if ((fd = open(file_name, O_WRONLY | O_CREAT, 0666)) < 0) {
		fprintf(stderr, "Can't create/open file %s\n", file_name);
		return EINVPNTR;
	}

	if (write(fd, (void *)header, strlen(header)) < 0) {
		fprintf(stderr, "Can't write to file %s\n", file_name);
		return EINVPNTR;
	}

	for (i = 1; i < tree->tree_sz; i++) {
		bzero((void *)my_str, 64);
		snprintf(my_str, 38, "n%03u -> n%03u;\nn%03u [label=\"%3.3f\"];\n",
				(i - 1) / P, i, i, tree->pa_nodes[i]->data);
		if (write(fd, (void *)my_str, strlen(my_str)) < 0) {
			fprintf(stderr, "Can't write to file %s\n", file_name);
			return EINVPNTR;
		}
	}
	bzero((void *)my_str, 64);
	snprintf(my_str, 2, "}\n");
	if (write(fd, (void *)my_str, strlen(my_str)) < 0) {
		fprintf(stderr, "Can't write to file %s\n", file_name);
		return EINVPNTR;
	}

	close(fd);
	return ESUCCESS;
}

