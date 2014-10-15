#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include <cauchy.h>

#define TIME 1000000
#define ACCURACY 0.0000000001
#define PRNT_CMX(x) (fabs(x) > ACCURACY ? x : 0)
#define PP(a) power((double)P, (double)a)

typedef struct Couchy_coeffs {
	int gamma;
	pa_num *n;
	int j;
	complex double Wgnj_x;
	complex double Cgnj_x;
	complex double Bgnj_x;
	complex double Agnj_x;
	complex double Lgnj_x;
} Couchy_coeffs;

// Bgnj_x = Integral(rho_backward(x,y)*WAVEgnj(y)dy) by Qp
// Cgnj_x = Integral(rho_forward(x,y)*WAVEgnj(x)dy) by Qp
// Bgnj_x - Cgnj_x

int cnt = 0;

double (*rho)(pa_num *x) = NULL;
double (*rho_backward)(pa_num *x, pa_num *y) = NULL;
double (*rho_forward)(pa_num *x, pa_num *y) = NULL;
double (*start_cond)(pa_num *pnum) = NULL;
int g_min = 0;
int g_max = -1;
int g_chy = -1;
pa_num **ns = NULL;

//#define ACC TRUE

#ifdef ACC
#define ALPHA 2
double *lambda = NULL;
#endif

PADIC_ERR do_for_j(int gamma, pa_num *n, int j, pa_num *x, Couchy_coeffs *array)
{
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
	double B_img = 0, B_rez = 0, B_fun = 0;
	double C_img = 0, C_rez = 0, C_fun = 0;
	double A_img = 0, A_rez = 0, A_fun = 0;
	complex double Bgnj_x = I, Cgnj_x = I, Agnj_x = I;

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
			return EMALLOC;
		}
		err = p_gamma_pa_num(pa, qs[i], g_max);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return err;
		}

		if ((rho_backward(x, pa) == INFINITY) || (rho_forward(x, pa) == INFINITY)) {
			fprintf(stderr, "AAAA!!! We are going to infinity!\n");
			return EINVPNTR;
		} else {
			B_fun = rho_backward(x, pa);
			C_fun = rho_forward(x, pa);
			A_fun = start_cond(pa);

			if (creal(array[cnt].Wgnj_x) != 0) {
				C_rez += creal(array[cnt].Wgnj_x) * C_fun;
				C_img += cimag(array[cnt].Wgnj_x) * C_fun;
			}
			B_rez += creal(wavelet(pa, n, -gamma, j)) * B_fun;
			B_img += cimag(wavelet(pa, n, -gamma, j)) * B_fun;
			A_rez += creal(wavelet(pa, n, -gamma, j)) * A_fun;
			A_img += cimag(wavelet(pa, n, -gamma, j)) * A_fun;
		}
		free_pa_num(pa);
		free_pa_num(qs[i]);
	}

	Bgnj_x = B_rez * power((double)P, -g_max) + I * B_img * power((double)P, -g_max);
	Cgnj_x = C_rez * power((double)P, -g_max) + I * C_img * power((double)P, -g_max);
	Agnj_x = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);

	array[cnt].Bgnj_x = Bgnj_x;
	array[cnt].Cgnj_x = Cgnj_x;
	array[cnt].Agnj_x = Agnj_x;

	if (fabs(creal(array[cnt].Wgnj_x)) > ACCURACY)
		array[cnt].Lgnj_x = (array[cnt].Bgnj_x - array[cnt].Cgnj_x) / array[cnt].Wgnj_x;
	else
		array[cnt].Lgnj_x = 0;

	if (fabs(cimag(array[cnt].Lgnj_x)) > ACCURACY) {
		fprintf(stderr, ">>> L <<< Something wrong!!!\n");
		//return EINVCOEFF;
	}

	return err;
}


PADIC_ERR do_for_n(int gamma, pa_num *n, pa_num *x, Couchy_coeffs *array)
{
	int j = 0;
	PADIC_ERR err = ESUCCESS;

	if (start_cond == NULL || rho_backward == NULL || rho_forward == NULL) {
		fprintf(stderr, "Invalid pointer to function\n");
		return EINVPNTR;
	}
	if (x == NULL || n == NULL || array == NULL) {
		fprintf(stderr, "Invalid pointer to data\n");
		return EINVPNTR;
	}

	if (gamma < g_min || gamma > g_chy) {
		fprintf(stderr, "Invalid Value gamma: %d\n", gamma);
		fprintf(stderr, "gamma should be greater than %d\n", g_min);
		fprintf(stderr, "gamma should be less than %d\n", g_max);
		fflush(stderr);
		return EGAMMAOUT;
	}

	for (j = 1; j < P; j++) {
		array[cnt].gamma = gamma;
		array[cnt].n = n;
		array[cnt].j = j;
		array[cnt].Wgnj_x = wavelet(x, n, -gamma, j);

		err = do_for_j(gamma, n, j, x, array);
		if (err != ESUCCESS)
			return err;
		fprintf(stdout, "%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\n",
			gamma, padic2double(n),
			PRNT_CMX(creal(array[cnt].Wgnj_x)), PRNT_CMX(cimag(array[cnt].Wgnj_x)),
			PRNT_CMX(creal(array[cnt].Bgnj_x)), PRNT_CMX(cimag(array[cnt].Bgnj_x)),
			PRNT_CMX(creal(array[cnt].Cgnj_x)), PRNT_CMX(cimag(array[cnt].Cgnj_x)),
			PRNT_CMX(creal(array[cnt].Lgnj_x)), PRNT_CMX(cimag(array[cnt].Lgnj_x)));
		cnt++;
	}
	return err;
}

PADIC_ERR get_integrals(pa_num *x, Couchy_coeffs *array)
{
	int gamma = 0, ns_sz = 0, i = 0;
	PADIC_ERR err = ESUCCESS;

	if (start_cond == NULL || rho_backward == NULL || rho_forward == NULL) {
		fprintf(stderr, "Invalid pointer to function\n");
		return EINVPNTR;
	}
	if (x == NULL || array == NULL) {
		fprintf(stderr, "Invalid pointer to data\n");
		return EINVPNTR;
	}

	// count Bgnj(x), Cgnj(x) and Agnj_x integrals
	for (gamma = g_min; gamma <= g_chy; gamma++) {
		ns_sz = (size_t)qspace_sz(g_min, gamma);
		ns = (pa_num **)malloc(ns_sz * sizeof(pa_num*));
		if (ns == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}

		err = gen_quotient_space(ns, g_min, gamma);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid generating qspace\n");
			return err;
		}
		for (i = 0; i < ns_sz; i++) {
			err = do_for_n(gamma, ns[i], x, array);
			if (err != ESUCCESS)
				return err;
//			free_pa_num(ns[i]);
		}
//		free(ns);
	}
	return err;
}

PADIC_ERR count_matrix(int array_sz, double **matrix,
			Couchy_coeffs *array, pa_num *x)
{
	int i = 0, j = 0, k = 0;
	int qs_sz = -1;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
	PADIC_ERR err = ESUCCESS;
	double A_img = 0, A_rez = I;
	double B_img = I, B_rez = I;
	complex double A = 0;
	complex double B = 0;
	double C = 0;


	for (i = 0; i < cnt; i++) {
		for (j = 0; j < cnt; j++) {
			A_rez += creal(wavelet(pa, array[j].n, -array[j].gamma, array[j].j)
				* array[i].Bgnj_x);
			A_img += cimag(wavelet(pa, array[j].n, -array[j].gamma, array[j].j)
				* array[i].Bgnj_x);
			B_rez += creal(wavelet(pa, array[j].n, -array[j].gamma, array[j].j)
				* array[i].Cgnj_x);
			B_img += cimag(wavelet(pa, array[j].n, -array[j].gamma, array[j].j)
				* array[i].Cgnj_x);
			free_pa_num(pa);

			A = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);
			B = B_rez * power((double)P, -g_max) - I * B_img * power((double)P, -g_max);

			fprintf(stdout, ">>> %.04g\t%.04g\t%.04g\t%.04g <<<\n",
			PRNT_CMX(creal(A)), PRNT_CMX(cimag(A)), PRNT_CMX(creal(B)), PRNT_CMX(cimag(B)));
/*
			if (fabs(cimag(A)) > ACCURACY)
				fprintf(stderr, ">>> A <<< Something wrong!!!\n");
			if (fabs(cimag(B)) > ACCURACY)
				fprintf(stderr, ">>> B <<< Something wrong!!!\n");
*/
			matrix[i][j] = creal(A) - creal(B);
			fprintf(stdout, "%g\t", PRNT_CMX(matrix[i][j]));
		}
		fprintf(stdout, "\n");
	}
	return err;
}

FILE *__prepare_file(const char *sfx)
{
	int fd = -1;
	char *srv_str = NULL;
	FILE *output;

	srv_str = (char *)calloc(pathconf(".", _PC_PATH_MAX), sizeof(char));
	if (srv_str == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return NULL;
	}

	snprintf(srv_str, 19, "./res/%05d-%s.dat\n", getpid(), sfx);

	if ((fd = open(srv_str, O_WRONLY | O_CREAT | O_EXCL, 0666)) < 0) {
		perror("Cannot open file");
		return NULL;
	}
	free(srv_str);

	if ((output = fdopen(fd, "w")) == NULL) {
		perror("Cannot open file");
		return NULL;
	}
	return output;
}

#ifdef ACC
void __get_lambdas()
{
	int l_sz = g_chy - g_min + 1;
	int R = -g_min;
	int r = -g_chy;
	int gamma = 0;
	lambda = (double *)malloc(l_sz * sizeof(double));
	for (gamma = r; gamma <= R; gamma++) {
		lambda[gamma - r] =	- PP(-ALPHA * gamma) *
			((double)(PP(ALPHA) - (1.0/P)) / (double)(PP(ALPHA) - 1))
					+ PP(-ALPHA * R) *
			((double)(1 - (1.0/P)) / (double)(PP(ALPHA) - 1));
	}
}

int __p_norm_power(pa_num *pa)
{
	int i = 0, res = INT_MAX;
	for (i = pa->g_min; i <= pa->g_max; i++) {
		if (get_x_by_gamma(pa, i) == 0) {
			continue;
		} else {
			res = i;
			break;
		}
	}
	if (res == INT_MAX)
		return 0;
	return -res;
}

double __get_accurate_value(pa_num *x, double t)
{
	double ret = 0;
	int nu = __p_norm_power(x);
	int R = -g_min;
	int r = -g_chy;
	int gamma = 0;

	__get_lambdas();
	for (gamma = r; gamma <= R; gamma++)
		ret += PP(-gamma) * exp(lambda[gamma - r] * t);
	ret *= (P - 1);
	ret += PP(-R);

	if (nu <= 0)
		return ret;

	if (((nu + 1) < r) || ((nu + 1) > R)) {
		fprintf(stderr, ">>> SMTH WRONG with Accurate solution nu = %d\n", nu);
		return -1;
	}
	ret = 0;
	for (gamma = nu + 1; gamma <= R; gamma++)
		ret += PP(-gamma) * exp(lambda[gamma - r] * t);
	ret *= (P - 1);
	ret += PP(-R);
	ret -= PP(-nu) * exp(lambda[nu - r] * t);
	return ret;
}
#endif // ACC

PADIC_ERR solve_problem(
		double (*rho_bw)(pa_num *x, pa_num *y),
		double (*rho_fw)(pa_num *x, pa_num *y),
		double (*start_cond0)(pa_num *pnum),
		int gmin, int gmax, int gchy,
		pa_num *x0)
{
	int i = 0, j = 0;
	PADIC_ERR err = ESUCCESS;
	FILE *output, *lnout;
	complex double Sum_Phi = 0, Phi_gnj_t = 0;
	int array_sz = 0;
	double t = 0, step = 0.001, max = 1;
	Couchy_coeffs *array = NULL;
	double **matrix = NULL;

	if (gmax < gmin) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", gmin);
		fprintf(stderr, "less or equal than g_max = %d\n", gmax);
		fflush(stderr);
		return EGAMMAOUT;
	}

	if (start_cond0 == NULL || rho_bw == NULL || rho_fw == NULL) {
		fprintf(stderr, "Invalid pointer to function\n");
		return EINVPNTR;
	}

	if ((output = __prepare_file("pl")) == NULL) {
		fprintf(stderr, "Cannot open file\n");
		return EINVPNTR;
	}

	if ((lnout = __prepare_file("ln")) == NULL) {
		fprintf(stderr, "Cannot open file\n");
		return EINVPNTR;
	}

	for (j = gmin; j <= gchy; j++)
		array_sz += qspace_sz(j + 1, gchy + 1) * (P - 1);

	array = (Couchy_coeffs *)malloc(array_sz * sizeof(Couchy_coeffs));
	if (array == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	matrix = (double **)malloc(array_sz * sizeof(double *));
	if (matrix == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	g_min = gmin;
	g_max = gmax;
	g_chy = gchy;
	rho_backward = rho_bw;
	rho_forward = rho_fw;
	start_cond = start_cond0;

	fprintf(stdout, "Get x = %g\n", padic2double(x0));
	print_pa_num(x0);

	fprintf(stdout, "gamma\tn\tWgnj(x)\t\tBgnj(x)\t\tCgnj(x)\t\tL(x)\n");

	err = get_integrals(x0, array);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Integrals\n");
		return err;
	}

	if (cnt != array_sz) {
		fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
		return EINVPNTR;
	}

	for (j = 0; j < array_sz; j++) {
		matrix[j] = (double *)malloc(array_sz * sizeof(double));
		if (matrix[j] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}
	}

	err = count_matrix(array_sz, matrix, array, x0);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Matrix\n");
		return err;
	}

/*
	step = 0.001;
	max = 1;
	for (t = step; t <= max; t += step) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(creal(array[i].Lgnj_x) * t);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		if (fabs(cimag(Sum_Phi)) > ACCURACY)
			fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
#ifdef ACC
		double acc = __get_accurate_value(x0, t);
		fprintf(output, "%.03g\t%g\t%g\n", t, creal(Sum_Phi), acc);
		fprintf(lnout, "%.03g\t%g\t%g\n", log(t), log(creal(Sum_Phi)), log(acc));
#else
		fprintf(output, "%.03g\t%g\n", t, creal(Sum_Phi));
		fprintf(lnout, "%.03g\t%g\n", log(t), log(creal(Sum_Phi)));
#endif
		if ((max - t < step) && (max < TIME)) {
			step *= 10;
			max *= 10;
		}
	}
*/
	fflush(output);
	return 0;
}

