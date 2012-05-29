#include "p-adic.h"

// len = g_max - g_min + 1

long fspace_sz(int g_min, int g_max)
{
	long ret = 1;
	int i;
	int xsz = g_max - g_min;

	for (i = 0; i < xsz; i++)
		ret = ret * P;
	return ret;
}

pa_num* init_pa_num(int g_min, int g_max)
{
	pa_num *ret;
	ret = (pa_num*)malloc(sizeof(pa_num));
	ret->g_min = g_min;
	ret->g_max = g_max;
	ret->err = EINIT;
	ret->sign = POS;
	ret->x = (int *)malloc(sizeof(int)*(g_max - g_min + 1));
	bzero((void *)ret->x, sizeof(int)*(g_max - g_min + 1));
	return ret;
}

void free_pa_num(pa_num *pa)
{
	free(pa->x);
	free(pa);
}

int get_x_by_gamma(pa_num* pa, int gamma)
{
	if (gamma < pa->g_min)
		return 0;
	else if (gamma > pa->g_max)
		return (pa->sign) ? (P - 1) : 0;
	return pa->x[gamma - pa->g_min];
}

int set_x_by_gamma(pa_num* pa, int gamma, int x)
{
	if (gamma < pa->g_min || gamma > pa->g_max) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", pa->g_min);
		fprintf(stderr, "gamma should be less than %d\n", pa->g_max);
		return pa->err = EGOUT;
	}
	pa->x[gamma - pa->g_min] = x;
	return EINIT;
}

pa_num* __extend_number(pa_num *pa, int g_min, int g_max)
{
	int i;
	int xtmp;
	pa_num *ret;

	ret = init_pa_num(g_min, g_max);

	if (g_min > pa->g_min) {
		printf("Invalid Value g_min: %d\n", g_min);
		printf("g_min should be greater than %d\n", pa->g_min);
		ret->err = EGOUT;
		exit(-1);
	}

	if (g_max < pa->g_max) {
		printf("Invalid Value g_max: %d\n", g_max);
		printf("g_max should be less than %d\n", pa->g_max);
		ret->err = EGOUT;
		exit(-1);
	}

	for (i = g_min; i <= g_max; i++) {
		xtmp = get_x_by_gamma(pa, i);
		set_x_by_gamma(ret, i, xtmp);
	}

	return ret;
}

pa_num* minus(pa_num *pa1, pa_num *pa2)
{
	pa_num *ret, *shrt, *pn1, *pn2;
	int min, max;
	int coeff, i, skip;
	int in_mind = 0;

	min = (pa1->g_min < pa2->g_min) ? pa1->g_min : pa2->g_min;
	max = (pa1->g_max > pa2->g_max) ? pa1->g_max : pa2->g_max;
	/* +1 for number's sign */
	ret = init_pa_num(min, max + 1);

	pn1 = __extend_number(pa1, ret->g_min, ret->g_max);
	pn2 = __extend_number(pa2, ret->g_min, ret->g_max);

	in_mind = 0;
	for (i = ret->g_min; i <= ret->g_max; i++) {
		(in_mind) ? (coeff = get_x_by_gamma(pn1, i) - get_x_by_gamma(pn2, i) - 1) : \
			(coeff = get_x_by_gamma(pn1, i) - get_x_by_gamma(pn2, i));
		if (coeff >= 0) {
			set_x_by_gamma(ret, i, coeff);
			in_mind = 0;
		} else {
			set_x_by_gamma(ret, i, coeff + P);
			in_mind = 1;
		}
	}

	if (get_x_by_gamma(ret, ret->g_max) == (P - 1)) {
		ret->sign = NEG;
	}

	for (i = ret->g_min; i <= ret->g_max; i++) {
		if (get_x_by_gamma(ret, i) != 0)
			break;
	}
	min = i;

	if (min >= ret->g_max) {
		free_pa_num(ret);
		free_pa_num(pn1);
		free_pa_num(pn2);
		return shrt = init_pa_num(0, 0);
	}

	skip = (ret->sign == POS) ? 0 : (P - 1);

	for (i = ret->g_max; i >= ret->g_min; i--) {
		if (get_x_by_gamma(ret, i) != skip)
			break;
	}
	max = (skip) ? i + 1 : i;

	shrt = init_pa_num(min, max);
	shrt->sign = ret->sign;
	memcpy((void *)shrt->x, (void *)&(ret->x[min - ret->g_min]), \
		(max - min + 1) * sizeof(int));

	free_pa_num(ret);
	free_pa_num(pn1);
	free_pa_num(pn2);
	return shrt;
}

pa_num** gen_factor_space(int g_min, int g_max)
{
	size_t fs_sz;
	pa_num **ret;
	int div, n, i, j, k, l, ngrp;
	int min = g_min - g_max - 1;

	fs_sz = (size_t)fspace_sz(g_min, g_max);
	ret = (pa_num **)malloc(fs_sz * sizeof(pa_num*));

	for (i = 0; i < fs_sz; i++) {
		ret[i] = init_pa_num(min, -1);
	}

	div = 1;
	for (k = -1; k >= min; k--) {
		n = 0;
		ngrp = fs_sz / (div * P) ;
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

pa_num* p_gamma_pa_num(pa_num *pa, int gamma)
{
	pa_num *ret;
	int i;

	ret = init_pa_num(pa->g_min + gamma, pa->g_max + gamma);

	for (i = ret->g_min; i <= ret->g_max; i++) {
		set_x_by_gamma(ret, i, get_x_by_gamma(pa, i - gamma));
	}

	return ret;
}

void print_pa_num(pa_num *pa)
{
	int i;
	for (i = pa->g_min; i <= pa->g_max; i++) {
		printf("x[%d] = %d\t", i, get_x_by_gamma(pa, i));
	}
	printf("\n");
}

int p_norma(pa_num *pa)
{
	int i, res;

	for (i = pa->g_min; i <= pa->g_max; i++) {
		if (get_x_by_gamma(pa, i) == 0) {
			continue;
		} else {
			res = i;
			break;
		}
	}

	return (i == pa->g_max + 1) ? 0 : pow(P, -res);
}

int indicator(pa_num *x, pa_num *n, int gamma)
{
	pa_num *mult, *subtr;
	int norma;

	mult = p_gamma_pa_num(x, gamma);
	subtr = minus(mult, n);
	norma = p_norma(subtr);

	/* |p^gamma*x - n| <= 1 => indicator = 1
	 * |kern| <= 1 => gamma >= 0
	 */

	free_pa_num(subtr);
	free_pa_num(mult);
	return (norma > 1) ? 0 : 1;
}

float from_canonic_to_float(pa_num *pa)
{
	int i;
	float ret = 0;
	float ppow;

	ppow = (float)pow(P, pa->g_min);
	for (i = pa->g_min; i <= pa->g_max; i++) {
		ret += (float)get_x_by_gamma(pa, i) * ppow;
		ppow = ppow * P;
	}

	return ret;
}

pa_num* get_fractional_part(pa_num *pa)
{
	pa_num *ret;
	int i;

	ret = init_pa_num(pa->g_min, -1);

	for (i = pa->g_min; i < 0; i++) {
		set_x_by_gamma(ret, i, get_x_by_gamma(pa, i));
	}
	return ret;
}

complex character(pa_num *pa)
{
	pa_num *fnum;
	complex ret;

	fnum = get_fractional_part(pa);
	ret = cexpf(2.0 * PI * I * \
			from_canonic_to_float(fnum));
	free_pa_num(fnum);
	return ret;
}

complex wavelet(pa_num *x, pa_num *n, int gamma, int j)
{
	pa_num *kern;
	pa_num *jkern;
	pa_num *mult, *subtr;
	complex ret;

	mult = p_gamma_pa_num(x, gamma);
	subtr = minus(mult, n);
	kern = p_gamma_pa_num(subtr, -1);
	jkern = smult(kern, j);
	ret = character(jkern) * indicator(x, n, gamma);
	free_pa_num(jkern);
	free_pa_num(kern);
	free_pa_num(mult);
	free_pa_num(subtr);

	return crealf(ret) + I * cimagf(ret);
}

/* workaround for wavelet: suppose (0 < j < P) and (pa > 0) */
pa_num* smult(pa_num *pa, int j)
{
	pa_num *ret;
	int i, tmp, in_mind = 0;

	ret = init_pa_num(pa->g_min, pa->g_max + 1);

	for (i = ret->g_min; i <= ret->g_max; i++) {
		tmp = get_x_by_gamma(pa, i) * j + in_mind;
		(tmp >= P) ? set_x_by_gamma(ret, i, tmp - P) : \
			set_x_by_gamma(ret, i, tmp);
		in_mind = (tmp >= P) ? 1 : 0;
	}

	return ret;
}

pa_num* mult(pa_num *pa1, pa_num *pa2)
{
	pa_num *ret;

	ret = init_pa_num(pa1->g_min + pa2->g_min, pa1->g_max + pa2->g_max);
	printf("Not implemented yet\n");

	return 0;
}

float integral(float (*func)(pa_num *pnum), int g_min, int g_max)
{
	float ret;
	int fs_sz, i;
	pa_num **fs;
	pa_num *pa;
	float (*pfunc)(pa_num *pnum);

	fs_sz = (size_t)fspace_sz(g_min, g_max);
	fs = gen_factor_space(g_min, g_max);
	pfunc = func;

	ret = 0;
	for (i = 0; i < fs_sz; i++) {
		pa = p_gamma_pa_num(fs[i], -g_min);
		ret += (float)pfunc(pa);
		free_pa_num(pa);
		free_pa_num(fs[i]);
	}
	free(fs);
	return ret * (float)pow(P, g_min);
}

complex wavelet_integral(pa_num *n, int gamma, int j, int g_min, int g_max)
{
	complex ret;
	float img = 0.f, rez = 0.f, tmp1, tmp2;
	int fs_sz, i, gmin;
	pa_num **fs;
	pa_num *pa;

	/* for integration on the "level - 1" */
	gmin = gamma - 1;
	fs_sz = (size_t)fspace_sz(gmin, g_max);
	fs = gen_factor_space(gmin, g_max);

	for (i = 0; i < fs_sz; i++) {
		pa = p_gamma_pa_num(fs[i], -gmin);
		tmp1 = crealf(wavelet(pa, n, gamma, j));
		tmp2 = cimagf(wavelet(pa, n, gamma, j));
		printf("iter %d -- wv = %f + i * %f\n\n", i, tmp1, tmp2);
		rez += tmp1;
		img += tmp2;
		free_pa_num(pa);
		free_pa_num(fs[i]);
	}
	ret = rez * (float)pow(P, g_min) + I * img * (float)pow(P, g_min);
	free(fs);
	return ret;

}

#if 0
pa_num* get_invert(pa_num *pa)
{
	pa_num *ret, *nil;

	nil = init_pa_num(0, 0);
	ret = minus(nil, pa);
	print_pa_num(ret);

	return ret;
}
#endif

