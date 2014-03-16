#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <p-arithm.h>

int arith_compare(pa_num *pa1, pa_num *pa2)
{
	int max =0 , min = 0, i = 0, ret = INT_MAX;
	if (pa1 == NULL || pa2 == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return INT_MAX;
	}

	min = (pa1->g_min < pa2->g_min) ? pa1->g_min : pa2->g_min;
	max = (pa1->g_max > pa2->g_max) ? pa1->g_max : pa2->g_max;

	ret = 0;
	for (i = max; i >= min; i--) {
		if (get_x_by_gamma(pa1, i) > get_x_by_gamma(pa2, i))
			ret = 1;
		else if (get_x_by_gamma(pa1, i) < get_x_by_gamma(pa2, i))
			ret = -1;
		else
			continue;
	}

	if (LOG_LEVEL >= 2) {
		char *s = "=";
		if (ret != 0)
			s = (ret > 0) ? ">" : "<";
		fprintf(stdout, "arith_compare: p1 %s p2\n", s);
		fprintf(stdout, "arith_compare: ");
		print_pa_num(pa1);
		fprintf(stdout, "arith_compare: ");
		print_pa_num(pa2);
	}

	return ret;
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

	if (LOG_LEVEL >= 2) {
		fprintf(stdout, "add: ");
		print_pa_num(pa1);
		fprintf(stdout, "add: ");
		print_pa_num(pa2);
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
		res->sign = (cmp < 0) ? pa1->sign : reverse_sign(pa1->sign);
		err = (cmp < 0) ? __dummy_sub(ext_pa, pa1, pa2) : \
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

	if (LOG_LEVEL >= 2) {
		fprintf(stdout, "add: %g + %g = %g\n",
			padic2double(pa1), padic2double(pa2), padic2double(res));
		fprintf(stdout, "add: ");
		print_pa_num(res);
	}

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

	if (LOG_LEVEL >= 2) {
		fprintf(stdout, "sub: ");
		print_pa_num(pa1);
		fprintf(stdout, "sub: ");
		print_pa_num(pa2);
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
		ext_pa->sign = (cmp < 0) ? pa1->sign : reverse_sign(pa1->sign);
		err = (cmp < 0) ? __dummy_sub(ext_pa, pa1, pa2) : \
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

	if (LOG_LEVEL >= 2) {
		fprintf(stdout, "sub: %g - %g = %g\n",
			padic2double(pa1), padic2double(pa2), padic2double(res));
		fprintf(stdout, "sub: ");
		print_pa_num(res);
	}

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

	if (LOG_LEVEL >= 2) {
		fprintf(stdout, "get_fractional_part: {%g} = %g\n",
			padic2double(pa), padic2double(res));
		fprintf(stdout, "get_fractional_part: ");
		print_pa_num(pa);
		fprintf(stdout, "get_fractional_part: ");
		print_pa_num(res);
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

	if (LOG_LEVEL >= 2) {
		fprintf(stdout, "jmult: %d*%g = %g\n",
			j, padic2double(pa), padic2double(res));
		fprintf(stdout, "jmult: ");
		print_pa_num(pa);
		fprintf(stdout, "jmult: ");
		print_pa_num(res);
	}

	return ESUCCESS;
}

