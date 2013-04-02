#include "../src/p-adic.h"

#define G_MAX (1)
#define G_MIN (-1)

static int gmin = G_MIN;
static int gmax = G_MAX;

// only for ALPHA = 2
float function(pa_num* pa)
{
	float ret = 0;
	if (pa == NULL) {
		fprintf(stderr, "Involid pointer\n");
		return -1;
	}
	ret = 1.0 / (p_norm(pa) * p_norm(pa));
	return (p_norm(pa) < power(P, -G_MAX)) ? \
			power(P, -G_MAX) : ret;
}

void do_for_n(pa_num *x, int gamma, pa_num *n)
{
	int j = 0;
	complex Cgnj_x = I;
	complex Wgnj_x = I;
	float (*pfunc)(pa_num *pnum);
	pfunc = function;

	for (j = 1; j < P; j++) {
		Cgnj_x = wavelet_integral_C_gnj_x(pfunc, x, n, gamma, j, gmin, gmax);
		Wgnj_x = wavelet_integral(pfunc, n, gamma, j, gmin, gmax);
		printf("%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\n", from_canonic_to_float(x),
				from_canonic_to_float(n), gamma, j, crealf(Cgnj_x), cimagf(Cgnj_x),
				crealf(Wgnj_x), cimagf(Wgnj_x));
	}
}

void do_for_gamma(pa_num *x, int gamma)
{
	int i = 0;
	PADIC_ERR err = ESUCCESS;
	int ns_sz = 0;
	pa_num **ns = NULL;

	ns_sz = (size_t)qspace_sz(gamma, gmax);
	ns = (pa_num **)malloc(ns_sz * sizeof(pa_num*));
	if (ns == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = gen_quotient_space(ns, gamma, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating qspace\n");
		exit(err);
	}

	for (i = 0; i < ns_sz; i++) {
		do_for_n(x, gamma, ns[i]);
		free_pa_num(ns[i]);
	}
	free(ns);
}

void do_for_x(pa_num *x)
{
	int gamma = 0;
	float bx = 0;
	float (*pfunc)(pa_num *pnum);
	pfunc = function;

	// count B(x) integral
	bx = integral_B_x(pfunc, x, gmin, gmax);
//	printf("B(x) = %f\n", bx);

	// count Cgnj(x) integral
	for (gamma = gmin; gamma < gmax; gamma++) {
		do_for_gamma(x, gamma);
	}
}

int main()
{
	int xs_sz = 0, i = 0;
	pa_num *x = NULL;
	pa_num **xs = NULL;
	PADIC_ERR err = ESUCCESS;

//	printf("Test#6: Wavelet integrals\n");
//	printf("Parameters: g_max = %d, g_min = %d\n", gmax, gmin);

	// x is in B_gamma_max
	xs_sz = (size_t)qspace_sz(gmin, gmax);
//	printf("x space sz = %d\n", xs_sz);
	printf("x\tn\tgamma\tj\tCgnj(x)\n");

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

		do_for_x(x);

		free_pa_num(x);
		free_pa_num(xs[i]);
	}
	free(xs);
	return 0;
}

