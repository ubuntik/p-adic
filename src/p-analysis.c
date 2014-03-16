#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <p-analysis.h>

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

	free_pa_num(jkern);
	free_pa_num(kern);
	free_pa_num(mult);
	free_pa_num(subtr);

	return ret;
}

double integral(double (*func)(pa_num *pnum), int g_min, int g_max)
{
	double ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
//	pa_num *spoint = NULL;
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

		if ( pfunc(pa) != INFINITY ) {

//			print_pa_num(pa);
//			printf("|%g| = %g\n", padic2double(pa), p_norm(pa));

			ret += (double)pfunc(pa);
			free_pa_num(qs[i]);
			free_pa_num(pa);
			continue;
		}

		fprintf(stderr, "Warning!!! Special point has found!\n");
		print_pa_num(pa);
#if 0
		fprintf(stderr, "%g is special point!\n", padic2double(pa));
		spoint = (pa_num *)malloc(sizeof(pa_num));
		if (spoint == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return err;
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
		fprintf(stderr, "Take x = %g\n", padic2double(spoint));
		print_pa_num(spoint);
		ret += (double)pfunc(spoint);
		free_pa_num(qs[i]);
		free_pa_num(pa);
		free_pa_num(spoint);
#endif
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
//	pa_num *spoint = NULL;
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
		printf("%g is special point!\n", padic2double(pa));
#if 0
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

		printf("Take x = %g\n", padic2double(spoint));
		print_pa_num(spoint);

		fun = pfunc(spoint);
		rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
		img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(qs[i]);
#endif
	}
	ret = rez * power((double)P, -g_max) + I * img * power((double)P, -g_max);
	free(qs);
	return ret;

}

double integral_B_x(double (*func)(pa_num *pnum), pa_num *x, int g_min, int g_max)
{
	double ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *res= NULL;
//	pa_num *spoint = NULL;
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
			ret += (double)pfunc(res);
			free_pa_num(qs[i]);
			free_pa_num(pa);
			free_pa_num(res);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
		print_pa_num(pa);
		printf("%g is special point!\n", padic2double(pa));
#if 0
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

		ret += (double)pfunc(res);
		free_pa_num(qs[i]);
		free_pa_num(pa);
		free_pa_num(spoint);
#endif
	}
	free(qs);
	return ret * power((double)P, -g_max);
}

complex double wavelet_integral_C_gnj_x(double (*func)(pa_num *pnum), pa_num *x,
		pa_num *n, int gamma, int j, int g_min, int g_max)
{
	complex double ret = I;
	PADIC_ERR err = ESUCCESS;
	double img = 0, rez = 0, fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *res = NULL;
//	pa_num *spoint = NULL;
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
#if 0
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

		printf("Take x = %g\n", padic2double(spoint));
		print_pa_num(spoint);

		fun = pfunc(spoint);
		rez += crealf(wavelet(res, n, -gamma, j)) * fun;
		img += cimagf(wavelet(res, n, -gamma, j)) * fun;
		free_pa_num(res);
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(qs[i]);
#endif
	}
	ret = rez * power((double)P, -g_max) + I * img * power((double)P, -g_max);
	free(qs);
	return ret;

}

complex double wavelet_integral_Agnj(double (*func)(pa_num *pnum), pa_num *n, int gamma, \
						int j, int g_min, int g_max)
{
	// start condition + sopryazhenniy wavelet
	complex double ret = I;
	PADIC_ERR err = ESUCCESS;
	double img = 0, rez = 0, fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
//	pa_num *spoint = NULL;
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
		printf("%g is special point!\n", padic2double(pa));

#if 0
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

		printf("Take x = %g\n", padic2double(spoint));
		print_pa_num(spoint);

		fun = pfunc(spoint);
		rez += crealf(wavelet(pa, n, -gamma, j)) * fun;
		img += cimagf(wavelet(pa, n, -gamma, j)) * fun;
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(qs[i]);
#endif
	}
	// sopryazhennie
	ret = rez * power((double)P, -g_max) - I * img * power((double)P, -g_max);
	free(qs);
	return ret;

}

