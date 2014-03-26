#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include <cauchy.h>

typedef struct Couchy_coeffs {
	int gamma;
	pa_num *n;
	int j;
	complex double Wgnj_x;
	complex double Cgnj_x;
	complex double Bgnj_x;
	complex double Agnj_x;
} Couchy_coeffs;

// Bgnj_x = Integral(rho(x-y)*WAVEgnj(y)dy) by Qp
// Cgnj_x = Integral(rho(x-y)*WAVEgnj(x)dy) by Qp
// Bgnj_x - Cgnj_x

// how many integrals we have to count (j * gamma * n)
int array_sz = 0;
Couchy_coeffs *array = NULL;
int cnt = 0;

pa_num *x = NULL;
double (*rho)(pa_num *pnum) = NULL;
double (*start_cond)(pa_num *pnum) = NULL;
int g_min = 0;
int g_max = -1;
int g_chy = -1;

PADIC_ERR do_for_j(int gamma, pa_num *n, int j)
{
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *res = NULL;
	double B_img = 0, B_rez = 0, B_fun = 0;
	double C_img = 0, C_rez = 0, C_fun = 0;
	double A_img = 0, A_rez = 0, A_fun = 0;
	complex double Bgnj_x = I;
	complex double Cgnj_x = I;
	complex double Agnj_x = I;


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

		// x - y
		res = (pa_num *)malloc(sizeof(pa_num));
		if (res == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}
		err = sub(res, x, pa);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid subtraction\n");
			return err;
		}

		if ((rho(pa) == INFINITY) || (rho(res) == INFINITY)) {
			fprintf(stderr, "AAAA!!! We are going to infinity!\n");
			return EINVPNTR;
		} else {
			B_fun = rho(res);
			C_fun = rho(res);
			A_fun = start_cond(pa);

			if (creal(array[cnt].Wgnj_x) != 0) {
				B_rez += creal(array[cnt].Wgnj_x) * B_fun;
				B_img += cimag(array[cnt].Wgnj_x) * B_fun;
			}
			C_rez += creal(wavelet(pa, n, -gamma, j)) * C_fun;
			C_img += cimag(wavelet(pa, n, -gamma, j)) * C_fun;
			A_rez += creal(wavelet(pa, n, -gamma, j)) * A_fun;
			A_img += cimag(wavelet(pa, n, -gamma, j)) * A_fun;
		}
		free_pa_num(pa);
		free_pa_num(res);
		free_pa_num(qs[i]);
	}

	Bgnj_x = B_rez * power((double)P, -g_max) + I * B_img * power((double)P, -g_max);
	Cgnj_x = C_rez * power((double)P, -g_max) + I * C_img * power((double)P, -g_max);
	Agnj_x = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);

/*
	if (cimag(Bgnj_x) > 0.000001) {
		fprintf(stderr, ">>> C <<< Something wrong!!!\n");
		return EINVCOEFF;
	}
	if (cimag(Cgnj_x) > 0.000001) {
		fprintf(stderr, ">>> C <<< Something wrong!!!\n");
		return EINVCOEFF;
	}
	if (cimag(Agnj_x) > 0.000001) {
		fprintf(stderr, ">>> A <<< Something wrong!!!\n");
		return EINVCOEFF;
	}

	array[cnt].Bgnj_x = creal(Bgnj_x);
	array[cnt].Cgnj_x = creal(Cgnj_x);
	array[cnt].Agnj_x = creal(Agnj_x);
*/

	array[cnt].Bgnj_x = Bgnj_x;
	array[cnt].Cgnj_x = Cgnj_x;
	array[cnt].Agnj_x = Agnj_x;
	return err;
}


PADIC_ERR do_for_n(int gamma, pa_num *n)
{
	int j = 0;
	PADIC_ERR err = ESUCCESS;

	if (start_cond == NULL || rho == NULL) {
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

/*
		if (cimag(array[cnt].Wgnj_x) > 0.000001) {
			fprintf(stderr, ">>> W >>> Something wrong!!!\n");
			return EINVCOEFF;
		}
*/

		err = do_for_j(gamma, n, j);
		if (err != ESUCCESS)
			return err;
		fprintf(stdout, "%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\n",
			gamma, padic2double(n),
			creal(array[cnt].Wgnj_x), cimag(array[cnt].Wgnj_x),
			creal(array[cnt].Bgnj_x), cimag(array[cnt].Bgnj_x),
			creal(array[cnt].Cgnj_x), cimag(array[cnt].Cgnj_x),
			creal(array[cnt].Bgnj_x - array[cnt].Cgnj_x),
			cimag(array[cnt].Bgnj_x - array[cnt].Cgnj_x));
		cnt++;
	}
	return err;
}

PADIC_ERR get_integrals()
{
	int gamma = 0, ns_sz = 0, i = 0;
	PADIC_ERR err = ESUCCESS;
	pa_num **ns = NULL;

	if (start_cond == NULL || rho == NULL) {
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
			err = do_for_n(gamma, ns[i]);
			if (err != ESUCCESS)
				return err;
			free_pa_num(ns[i]);
		}
		free(ns);
	}
	return err;
}

PADIC_ERR solve_problem(
		double (*rho0)(pa_num *pnum),
		double (*start_cond0)(pa_num *pnum),
		int gmin, int gmax, int gchy)
{
	int fd = -1, i = 0;
	char *srv_str = NULL;
	pa_num *x0 = NULL;
	PADIC_ERR err = ESUCCESS;
	double Sum_Phi = 0, Phi_gnj_t = 0;
	float t = 0.01;
	int array_sz = 0, gamma = 0;
	FILE *output;

	if (gmax < gmin) {
		fprintf(stderr, "Invalid gammas' values:\n");
		fprintf(stderr, "g_min = %d should be ", gmin);
		fprintf(stderr, "less or equal than g_max = %d\n", gmax);
		fflush(stderr);
		return EGAMMAOUT;
	}

	if (start_cond0 == NULL || rho0 == NULL) {
		fprintf(stderr, "Invalid pointer to function\n");
		return EINVPNTR;
	}

	// prepare file for results
	srv_str = (char *)calloc(pathconf(".", _PC_PATH_MAX), sizeof(char));
	if (srv_str == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}

	snprintf(srv_str, 16, "./res/%05d.xls\n", getpid());

	if ((fd = open(srv_str, O_WRONLY | O_CREAT | O_EXCL, 0666)) < 0) {
		perror("Cannot open file");
		exit(errno);
	}
	free(srv_str);

	if ((output = fdopen(fd, "w")) == NULL) {
		perror("Cannot open file");
		exit(errno);
	}

	for (gamma = gmin; gamma <= gchy; gamma++)
		array_sz += qspace_sz(gamma + 1, gchy + 1) * (P - 1);

	array = (Couchy_coeffs *)malloc(array_sz * sizeof(Couchy_coeffs));
	if (array == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}

	// get x
	x0 = (pa_num *)malloc(sizeof(pa_num));
	if (x0 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}

	err = init_pa_num(x0, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}

/*
	err = set_x_by_gamma(x0, -1, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
*/

	fprintf(stdout, "Get x = %g\n", padic2double(x0));
	print_pa_num(x0);

	g_min = gmin;
	g_max = gmax;
	g_chy = gchy;
	x = x0;
	rho = rho0;
	start_cond = start_cond0;

	fprintf(stdout, "gamma\tn\tWgnj(x)\t\t\tBgnj(x)\t\t\tCgnj(x)\t\t\tB-C(x)\n");

	err = get_integrals();
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Integrals\n");
		exit(err);
	}

	if (cnt != array_sz) {
		fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
		exit(-1);
	}

	for (t = 0.01; t < 1; t += 0.01) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(-(array[i].Bgnj_x - array[i].Cgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, Sum_Phi);
	}

	for (t = 1; t < 10000; t += 1) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(-(array[i].Bgnj_x - array[i].Cgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, Sum_Phi);
	}


	for (t = 10000; t < 150000; t += 10) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(-(array[i].Bgnj_x - array[i].Cgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, Sum_Phi);
	}

	for (t = 150000; t < 10000000; t += 100) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(-(array[i].Bgnj_x - array[i].Cgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, Sum_Phi);
	}


	fflush(output);
	free_pa_num(x0);
	close(fd);

	return 0;
}

