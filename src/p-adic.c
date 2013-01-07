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
	pa->x = (int *)malloc(sizeof(int) * (pa->g_max - pa->g_min + 1));
	if (pa->x == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	bzero((void *)pa->x, sizeof(int) * (pa->g_max - pa->g_min + 1));
	return ESUCCESS;
}

void free_pa_num(pa_num* pa)
{
	if (pa == NULL)
		fprintf(stderr, "Invalid pointer\n");
	free(pa->x);
	free(pa);
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

pa_num** gen_quotient_space(int g_min, int g_max)
{
	size_t qs_sz;
	pa_num **ret;
	int div, n, i, j, k, l, ngrp;
	int min = g_min - g_max - 1;
	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return NULL;
	}

	qs_sz = (size_t)qspace_sz(g_min, g_max);
	ret = (pa_num **)malloc(qs_sz * sizeof(pa_num*));

	for (i = 0; i < qs_sz; i++) {
		ret[i] = (pa_num *)malloc(sizeof(pa_num));
		if (ret[i] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return NULL;
		}
		init_pa_num(ret[i], min, -1);
	}

	div = 1;
	for (k = -1; k >= min; k--) {
		n = 0;
		ngrp = qs_sz / (div * P) ;
		for (j = 0; j < ngrp; j++) {
			for (i = 0; i < P; i++) {
				for (l = 0; l < div; l++) {
					set_x_by_gamma(ret[n++], k, i);
				}
			}
		}
		div = div * P;
	}

	return ret;
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

double power(double base, double exponent)
{
	double half_pow = 0;
	if (exponent == 0)
		return (double)1;
	else if (exponent < 0)
		return (double)1 / power(base, -exponent);
	else if (fmod(exponent, 2) == 0) {
		half_pow = power(base, exponent / (double)2);
		return half_pow * half_pow;
	}
	else
		return base * power(base, exponent - 1);
}

double p_norm(pa_num *pa)
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
		return 0;

	return power(P, -res);
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

double from_canonic_to_double(pa_num *pa)
{
	int i = 0;
	double ret = 0;
	double ppow = 0;

	if (pa == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}

	ppow = power(P, pa->g_min);
	for (i = pa->g_min; i <= pa->g_max; i++) {
		ret += (double)get_x_by_gamma(pa, i) * ppow;
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

	ret = cexp(2 * PI * I * \
			from_canonic_to_double(fnum));
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
	ret = character(jkern) * indicator(x, n, gamma);

	free_pa_num(jkern);
	free_pa_num(kern);
	free_pa_num(mult);
	free_pa_num(subtr);

//	return creal(ret) + I * cimag(ret);
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

double integral(double (*func)(pa_num *pnum), int g_min, int g_max)
{
	double ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **fs = NULL;
	pa_num *pa = NULL, *spoint = NULL;
	double (*pfunc)(pa_num *pnum) = NULL;

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
	fs = gen_quotient_space(g_min, g_max);
	pfunc = func;

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, fs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}

		if ( pfunc(pa) != INFINITY ) {
			ret += (double)pfunc(pa);
			free_pa_num(fs[i]);
			free_pa_num(pa);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%f is special point!\n", from_canonic_to_double(pa));
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
		printf("Take x = %f\n", from_canonic_to_double(spoint));
		print_pa_num(spoint);
		ret += (double)pfunc(spoint);
		free_pa_num(fs[i]);
		free_pa_num(pa);
		free_pa_num(spoint);
	}
	free(fs);
	return ret * (double)power(P, -g_max);
}

complex wavelet_integral(double (*func)(pa_num *pnum), pa_num *n, int gamma, \
						int j, int g_min, int g_max)
{
	complex ret = I;
	PADIC_ERR err = ESUCCESS;
	double img = 0, rez = 0, fun, wav1, wav2, res1, res2, fpa;
	int qs_sz, i;
	pa_num **fs;
	pa_num *pa, *spoint;
	double (*pfunc)(pa_num *pnum);

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
	fs = gen_quotient_space(g_min, g_max);

	for (i = 0; i < qs_sz; i++) {
		pa = (pa_num *)malloc(sizeof(pa_num));
		if (pa == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return ret;
		}
		err = p_gamma_pa_num(pa, fs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return ret;
		}
		print_pa_num(pa);
		printf("%f is represetatives!\n\n", from_canonic_to_double(pa));
		/* workaround for special point (f.i. x = 0) */
		if ( pfunc(pa) != INFINITY ) {
			fun = pfunc(pa);
			fpa = from_canonic_to_double(pa);
			printf("func(%f) = %f\n", fpa, fun);
			wav1 = creal(wavelet(pa, n, gamma, j));
			wav2 = cimag(wavelet(pa, n, gamma, j));
			printf("pa = %f -- wav1 = %f -- wav2 = %f\n\n", fpa, wav1, wav2);
			res1 = creal(wavelet(pa, n, gamma, j)) * fun;
			res2 = cimag(wavelet(pa, n, gamma, j)) * fun;
			printf("pa %f -- res = %f + i * %f\n\n", fpa, res1, res2);
			rez += res1;
			img += res2;
			free_pa_num(pa);
			free_pa_num(fs[i]);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%f is special point!\n", from_canonic_to_double(pa));

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

		printf("Take x = %f\n", from_canonic_to_double(spoint));
		print_pa_num(spoint);

		printf("Check wavelet\n");
		wav1 = creal(wavelet(pa, n, gamma, j));
		wav2 = cimag(wavelet(pa, n, gamma, j));
		printf("wav from pa: %f -- %f\n", wav1, wav2);

		wav1 = creal(wavelet(spoint, n, gamma, j));
		wav2 = cimag(wavelet(spoint, n, gamma, j));
		printf("wav from pa: %f -- %f\n", wav1, wav2);

		fun = pfunc(spoint);
		fpa = from_canonic_to_double(spoint);
		printf("func(%f) = %f\n", fpa, fun);
		wav1 = creal(wavelet(pa, n, gamma, j));
		wav2 = cimag(wavelet(pa, n, gamma, j));
		printf("pa = %f -- wav1 = %f -- wav2 = %f\n\n", fpa, wav1, wav2);
		res1 = creal(wavelet(pa, n, gamma, j)) * fun;
		res2 = cimag(wavelet(pa, n, gamma, j)) * fun;
		printf("pa %f -- res = %f + i * %f\n\n", fpa, res1, res2);
		rez += res1;
		img += res2;
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(fs[i]);
	}
	ret = rez * (double)power(P, -g_max) + I * img * (double)power(P, -g_max);
	free(fs);
	return ret;

}

