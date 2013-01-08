#include "../src/p-adic.h"

#define G_MAX (0)
#define G_MIN (-1)

static int gmin = G_MIN;
static int gmax = G_MAX;

// only for ALPHA = 2
double function(pa_num* pa)
{
	double ret = 0;
	if (pa == NULL) {
		fprintf(stderr, "Involid pointer\n");
		return -1;
	}
	ret = 1.0 / (p_norm(pa) * p_norm(pa));
	return (p_norm(pa) < power((double)P, (double)-G_MAX)) ? \
			power((double)P, (double)(2 * G_MAX)) : ret;
}

double integral_B_x(pa_num *x, int g_min, int g_max)
{
	double ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL, *res= NULL;
	double (*pfunc)(pa_num *pnum) = NULL;

	if (g_max < g_min) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", g_min);
		fprintf(stderr, "less or equal than g_max = %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	pfunc = function;

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

//		print_pa_num(res);
//		printf("%f is sub!\n", from_canonic_to_double(res));
//
		if ( pfunc(res) != INFINITY ) {
			ret += (double)pfunc(res);
			free_pa_num(qs[i]);
			free_pa_num(pa);
			free_pa_num(res);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
//		print_pa_num(pa);
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
//		printf("Take x = %f\n", from_canonic_to_double(spoint));
//		print_pa_num(spoint);

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
//		print_pa_num(res);
//		printf("%f is sub!\n", from_canonic_to_double(res));

		ret += (double)pfunc(res);
		free_pa_num(qs[i]);
		free_pa_num(pa);
		free_pa_num(spoint);
	}
	free(qs);
	return ret * (double)power(P, -g_max);
}

complex wavelet_integral_C_gnj_x(pa_num *x, pa_num *n, int gamma, int j, \
							int g_min, int g_max)
{
	complex ret = I;
	PADIC_ERR err = ESUCCESS;
	double img = 0, rez = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *spoint = NULL, *res = NULL;
	double (*pfunc)(pa_num *pnum);

	if (n == NULL) {
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

	pfunc = function;
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
//		print_pa_num(pa);
//		printf("%f is represetatives!\n\n", from_canonic_to_double(pa));

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

//		print_pa_num(res);
//		printf("%f is sub!\n", from_canonic_to_double(res));
//
		/* workaround for special point (f.i. x = 0) */
		if ( pfunc(res) != INFINITY ) {
			rez += creal(wavelet(res, n, gamma, j)) * pfunc(res);
			img += cimag(wavelet(res, n, gamma, j)) * pfunc(res);
			free_pa_num(pa);
			free_pa_num(res);
			free_pa_num(qs[i]);
			continue;
		}

		printf("Warning!!! Special point has found!\n");
//		print_pa_num(pa);
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
//		printf("Take x = %f\n", from_canonic_to_double(spoint));
//		print_pa_num(spoint);

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
//		print_pa_num(res);
//		printf("%f is sub!\n", from_canonic_to_double(res));

		rez += creal(wavelet(res, n, gamma, j)) * pfunc(res);
		img += cimag(wavelet(res, n, gamma, j)) * pfunc(res);
		free_pa_num(pa);
		free_pa_num(spoint);
		free_pa_num(res);
		free_pa_num(qs[i]);
	}
	ret = rez * (double)power(P, -g_max) + I * img * (double)power(P, -g_max);
	free(qs);
	return ret;
}

void count_integral_C(pa_num *x)
{
}

void do_for_n(pa_num *x, int gamma, pa_num *n)
{
	int j = 0;
	complex Cgnj_x = I;

	for (j = 1; j < P; j++) {
		printf("j = %d\n", j);
		Cgnj_x = wavelet_integral_C_gnj_x(x, n, gamma, j, gmin,gmax);
		printf("Cgnj(x): %g + i * %g\n\n", creal(Cgnj_x), cimag(Cgnj_x));
	}
}

void do_for_gamma(pa_num *x, int gamma)
{
	int i = 0;
	pa_num *n = NULL;
	PADIC_ERR err = ESUCCESS;
	int ns_sz = 0;
	pa_num **ns = NULL;

	ns_sz = (size_t)qspace_sz(gmin, gamma);
	ns = (pa_num **)malloc(ns_sz * sizeof(pa_num*));
	if (ns == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = gen_quotient_space(ns, gmin, gamma);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating qspace\n");
		exit(err);
	}

	for (i = 0; i < ns_sz; i++) {
		n = (pa_num *)malloc(sizeof(pa_num));
		if (n == NULL) {
			fprintf(stderr,"Cannot alloc memory\n");
			exit(-1);
		}
		err = p_gamma_pa_num(n, ns[i], gamma);
		if (err != ESUCCESS) {
			fprintf(stderr,"Involid mult pgamma\n");
			exit(err);
		}
		printf("n = %f\n", from_canonic_to_double(n));
		print_pa_num(n);

		do_for_n(x, gamma, n);

		free_pa_num(n);
		free_pa_num(ns[i]);
	}
	free(ns);
}

void do_for_x(pa_num *x)
{
	int gamma = 0;
	double bx = 0;

	// count B(x) integral
	bx = integral_B_x(x, gmin, gmax);
	printf("B(x) = %f\n", bx);

	// count Cgnj(x) integral
	for (gamma = gmin; gamma <= gmax; gamma++) {
		printf("gamma = %d\n", gamma);
		do_for_gamma(x, gamma);
	}
}

int main()
{
	int xs_sz = 0, i = 0;
	pa_num *x = NULL;
	pa_num **xs = NULL;
	PADIC_ERR err = ESUCCESS;

	printf("Test#6: Wavelet integrals\n");
	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);

	// x is in B_gamma_max
	xs_sz = (size_t)qspace_sz(gmin, gmax);
	printf("x space sz = %d\n", xs_sz);

	xs = (pa_num **)malloc(xs_sz * sizeof(pa_num*));
	if (xs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = gen_quotient_space(xs, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		exit(err);
	}

	for (i = 0; i < xs_sz; i++) {
		x = (pa_num *)malloc(sizeof(pa_num));
		if (x == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			exit(-1);
		}
		err = p_gamma_pa_num(x, xs[i], gmax);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid mult p gamma\n");
			exit(err);
		}
		
		printf("\nx = %f\n", from_canonic_to_double(x));
		print_pa_num(x);

		do_for_x(x);

		free_pa_num(x);
		free_pa_num(xs[i]);
	}
	free(xs);
	//ret = wavelet_integral(pfunc, n, gamma, 1, gmin, gmax);
	return 0;
}

