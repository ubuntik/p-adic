#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <p-arithm.h>

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
	memcpy((void *)shrt->x, (void *)(&pa->x[min - pa->g_min]),
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

// TODO: for DEBUG -- add args and res print for add() and sub()

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

