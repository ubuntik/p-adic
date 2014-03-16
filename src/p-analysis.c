#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <p-analysis.h>

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

#if 0
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
#endif	// if 0

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
	if (res == INT_MAX) {
		if (LOG_LEVEL >= 1) {
			fprintf(stdout, "p_norm: |0| = 0\n");
			fprintf(stdout, "p_norm: ");
			print_pa_num(pa);
		}
		return (double)0;
	}

	if (LOG_LEVEL >= 1) {
		fprintf(stdout, "p_norm: P^%d = %g\n", -res, power((double)P, -res));
		fprintf(stdout, "p_norm: ");
                print_pa_num(pa);
        }
	return power((double)P, -res);
}

int indicator(pa_num *x, pa_num *n, int gamma)
{
	pa_num *mult = NULL, *subtr = NULL;
	double norm = 0;
	PADIC_ERR err = ESUCCESS;

	if (x == NULL || n == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return -1;
	}

	if (LOG_LEVEL >= 1) {
		fprintf(stdout, "p_norm: x = %g, n = %g, gamma = %d\n",
			padic2double(x), padic2double(n), gamma);
		fprintf(stdout, "p_norm: ");
		print_pa_num(x);
		fprintf(stdout, "p_norm: ");
		print_pa_num(n);
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

	if (LOG_LEVEL >= 1) {
		fprintf(stdout, "indicator: |p^gamma*x - n| <= 1 => indicator = 1\n");
		fprintf(stdout, "indicator: |%g| = %g => ind = %d\n",
				padic2double(subtr), norm, (norm <= 1) ? 1 : 0);
		fprintf(stdout, "indicator: ");
		print_pa_num(subtr);
        }

	free_pa_num(subtr);
	free_pa_num(mult);
	return (norm <= 1) ? 1 : 0;
}

complex double character(pa_num *pa)
{
	pa_num *fnum = NULL;
	complex double ret = I;
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

	ret = cexp(2 * PI * I *
			padic2double(fnum));

	if (LOG_LEVEL >= 1) {
		fprintf(stdout, "character: char(%g) = %g + i*%g\n",
			padic2double(pa), creal(ret), cimag(ret));
		fprintf(stdout, "character: 2pi{x} = %g\n", 2 * PI * padic2double(fnum));
		fprintf(stdout, "character: ");
		print_pa_num(pa);
        }
	free_pa_num(fnum);
	return ret;
}


complex double wavelet(pa_num *x, pa_num *n, int gamma, int j)
{
	pa_num *kern = NULL, *jkern = NULL;
	pa_num *mult = NULL, *subtr = NULL;
	complex double ret = I;
	PADIC_ERR err = ESUCCESS;

	if (x == NULL || n == NULL) {
		fprintf(stderr, "Invalid pointer\n");
		return ret;
	}
	if (j <= 0 || j >= P) {
		fprintf(stderr, "Invalid value of j\n");
		return ret;
	}

	if (LOG_LEVEL >= 1) {
		fprintf(stdout, "wavelet: x = %g, n = %g, gamma = %d, j = %d\n",
			padic2double(x), padic2double(n), gamma, j);
		fprintf(stdout, "wavelet: ");
		print_pa_num(x);
		fprintf(stdout, "wavelet: ");
		print_pa_num(n);
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

	if (indicator(x, n, gamma)) {
		ret = creal(character(jkern)) * power(P, -(double)gamma / 2) +
			I * cimag(character(jkern)) * power(P, -(double)gamma / 2);
	} else
		ret = 0 + I * 0;

	if (LOG_LEVEL >= 1)
		fprintf(stdout, "wavelet: W(%g) = %g + i*%g\n",
			padic2double(jkern), creal(ret), cimag(ret));

	free_pa_num(jkern);
	free_pa_num(kern);
	free_pa_num(mult);
	free_pa_num(subtr);

	return ret;
}

PADIC_ERR __get_spoint(pa_num *pa, pa_num *spoint)
{
	PADIC_ERR err = ESUCCESS;

	fprintf(stderr, "Warning!!! Special point has found: %g!\n",
				padic2double(pa));
	print_pa_num(pa);

	err = __extend_number(spoint, pa, pa->g_min, pa->g_max + 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid extending number\n");
		return err;
	}
	err = set_x_by_gamma(spoint, pa->g_max + 1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid setting coeff\n");
		return err;
	}

	fprintf(stderr, "Take x = %g\n", padic2double(spoint));
	print_pa_num(spoint);

	return err;
}

double integral(double (*func)(pa_num *pnum), int g_min, int g_max)
{
	double ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
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
			return EMALLOC;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return err;
		}

		if ( pfunc(pa) == INFINITY ) {
			pa_num *spoint = NULL;
			spoint = (pa_num *)malloc(sizeof(pa_num));
			if (spoint == NULL) {
				fprintf(stderr, "Cannot alloc memory\n");
				return err;
			}
			err = __get_spoint(pa, spoint);
			if (err != ESUCCESS) {
				fprintf(stderr, "Failed to get special point\n");
				return err;
			}
			ret += (double)pfunc(spoint);
			free_pa_num(spoint);
		} else
			ret += (double)pfunc(pa);

		free_pa_num(qs[i]);
		free_pa_num(pa);
	}
	free(qs);
	return ret * power((double)P, -g_max);
}

complex double wavelet_integral(double (*func)(pa_num *pnum), pa_num *n, int gamma, \
						int j, int g_min, int g_max)
{
	complex double ret = I;
	PADIC_ERR err = ESUCCESS;
	double img = 0, rez = 0, fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
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

		if ( pfunc(pa) == INFINITY ) {
			pa_num *spoint = NULL;
			spoint = (pa_num *)malloc(sizeof(pa_num));
			if (spoint == NULL) {
				fprintf(stderr, "Cannot alloc memory\n");
				return err;
			}
			err = __get_spoint(pa, spoint);
			if (err != ESUCCESS) {
				fprintf(stderr, "Failed to get special point\n");
				return err;
			}
			fun = pfunc(spoint);
			rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
			img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
			free_pa_num(spoint);
		} else {
			fun = pfunc(pa);
			rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
			img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
		}
		free_pa_num(pa);
		free_pa_num(qs[i]);
	}
	ret = rez * power((double)P, -g_max) + I * img * power((double)P, -g_max);
	free(qs);
	return ret;
}

