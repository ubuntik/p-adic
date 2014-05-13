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
	complex double Lgnj_x;
} Couchy_coeffs;

// Bgnj_x = Integral(rho(x-y)*WAVEgnj(y)dy) by Qp
// Cgnj_x = Integral(rho(x-y)*WAVEgnj(x)dy) by Qp
// Bgnj_x - Cgnj_x

// how many integrals we have to count (j * gamma * n)
Couchy_coeffs **array = NULL;
int cnt = 0;

double (*rho)(pa_num *pnum) = NULL;
double (*start_cond)(pa_num *pnum) = NULL;
int g_min = 0;
int g_max = -1;
int g_chy = -1;

PADIC_ERR do_for_j(int gamma, pa_num *n, int j, pa_num *x, Couchy_coeffs *array)
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
				C_rez += creal(array[cnt].Wgnj_x) * C_fun;
				C_img += cimag(array[cnt].Wgnj_x) * C_fun;
			}
			B_rez += creal(wavelet(pa, n, -gamma, j)) * B_fun;
			B_img += cimag(wavelet(pa, n, -gamma, j)) * B_fun;
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

	array[cnt].Bgnj_x = Bgnj_x;
	array[cnt].Cgnj_x = Cgnj_x;
	array[cnt].Agnj_x = Agnj_x;

	if (fabs(creal(array[cnt].Wgnj_x)) > 0.0000000001)
		array[cnt].Lgnj_x = (array[cnt].Bgnj_x - array[cnt].Cgnj_x) / array[cnt].Wgnj_x;
	else
		array[cnt].Lgnj_x = 0;

	if (fabs(cimag(array[cnt].Lgnj_x)) > 0.0000000001) {
		fprintf(stderr, ">>> L <<< Something wrong!!!\n");
		return EINVCOEFF;
	}

	return err;
}


PADIC_ERR do_for_n(int gamma, pa_num *n, pa_num *x, Couchy_coeffs *array)
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

		err = do_for_j(gamma, n, j, x, array);
		if (err != ESUCCESS)
			return err;
		fprintf(stdout, "%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\n",
			gamma, padic2double(n),
			fabs(creal(array[cnt].Wgnj_x)) > 0.0000000001 ? creal(array[cnt].Wgnj_x) : 0,
			fabs(cimag(array[cnt].Wgnj_x)) > 0.0000000001 ? cimag(array[cnt].Wgnj_x) : 0,
			fabs(creal(array[cnt].Bgnj_x)) > 0.0000000001 ? creal(array[cnt].Bgnj_x) : 0,
			fabs(cimag(array[cnt].Bgnj_x)) > 0.0000000001 ? cimag(array[cnt].Bgnj_x) : 0,
			fabs(creal(array[cnt].Cgnj_x)) > 0.0000000001 ? creal(array[cnt].Cgnj_x) : 0,
			fabs(cimag(array[cnt].Cgnj_x)) > 0.0000000001 ? cimag(array[cnt].Cgnj_x) : 0,
			fabs(creal(array[cnt].Lgnj_x)) > 0.0000000001 ? creal(array[cnt].Lgnj_x) : 0,
			fabs(cimag(array[cnt].Lgnj_x)) > 0.0000000001 ? cimag(array[cnt].Lgnj_x) : 0);

		cnt++;
	}
	return err;
}

PADIC_ERR get_integrals(pa_num *x, Couchy_coeffs *array)
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
			err = do_for_n(gamma, ns[i], x, array);
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
		int gmin, int gmax, int gchy,
		int ini_gamma, pa_num *ini_n)
{
	int fd = -1, i = 0, j = 0;
	char *srv_str = NULL;
	pa_num *x0 = NULL;
	PADIC_ERR err = ESUCCESS;
	complex double Sum_Phi = 0, Phi_gnj_t = 0;
	float t = 0.01;
	int array_sz = 0, gamma = 0;
	FILE *output;
	pa_num **xs = NULL;
	int x_sz = -1;

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

	x_sz = (size_t)qspace_sz(ini_gamma, gmax);
	xs = (pa_num **)malloc(x_sz * sizeof(pa_num*));
	if (xs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	err = gen_quotient_space(xs, ini_gamma, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating qspace\n");
		return err;
	}

	array = (Couchy_coeffs **)malloc(x_sz * sizeof(Couchy_coeffs *));
	if (array == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	for (i = 0; i < x_sz; i++) {
		x0 = (pa_num *)malloc(sizeof(pa_num));
		if (x0 == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}

		err = p_gamma_pa_num(x0, xs[i], gmax);
		if (err != ESUCCESS) {
			fprintf(stderr, "Invalid multiplication on p-gamma\n");
			return err;
		}

		fprintf(stdout, "Get x = %g\n", padic2double(x0));
		print_pa_num(x0);

		array[i] = (Couchy_coeffs *)malloc(array_sz * sizeof(Couchy_coeffs));
		if (array[i] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			exit(-1);
		}

		g_min = gmin;
		g_max = gmax;
		g_chy = gchy;
		rho = rho0;
		start_cond = start_cond0;

		fprintf(stdout, "gamma\tn\tWgnj(x)\t\tBgnj(x)\t\tCgnj(x)\t\tL(x)\n");

		err = get_integrals(x0, array[i]);
		if (err != ESUCCESS) {
			fprintf(stderr, "Failed to count Integrals\n");
			exit(err);
		}

		if (cnt != array_sz) {
			fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
			exit(-1);
		}

		cnt = 0;
	}

	for (t = 0.01; t < 1; t += 0.01) {
		Sum_Phi = x_sz * power(P, (double)gmin);
		for (j = 0; j < x_sz; j++) {
			for (i = 0; i < array_sz; i++) {
				Phi_gnj_t = array[j][i].Agnj_x * exp(creal(array[j][i].Lgnj_x) * t);
				Sum_Phi += Phi_gnj_t * array[j][i].Wgnj_x;
			}

			if (fabs(cimag(Sum_Phi)) > 0.0000000001) {
				fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
			}
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi) * power((double)P, -gmax));
	}

	for (t = 1; t < 1000; t += 1) {
		Sum_Phi = x_sz * power(P, (double)gmin);
		for (j = 0; j < x_sz; j++) {
			for (i = 0; i < array_sz; i++) {
				Phi_gnj_t = array[j][i].Agnj_x * exp(creal(array[j][i].Lgnj_x) * t);
				Sum_Phi += Phi_gnj_t * array[j][i].Wgnj_x;
			}

			if (fabs(cimag(Sum_Phi)) > 0.0000000001) {
				fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
			}
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi) * power((double)P, -gmax));
	}


	for (t = 1000; t < 10000; t += 10) {
		Sum_Phi = x_sz * power(P, (double)gmin);
		for (j = 0; j < x_sz; j++) {
			for (i = 0; i < array_sz; i++) {
				Phi_gnj_t = array[j][i].Agnj_x * exp(creal(array[j][i].Lgnj_x) * t);
				Sum_Phi += Phi_gnj_t * array[j][i].Wgnj_x;
			}

			if (fabs(cimag(Sum_Phi)) > 0.0000000001) {
				fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
			}
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi) * power((double)P, -gmax));
	}

	for (t = 10000; t < 100000; t += 100) {
		Sum_Phi = x_sz * power(P, (double)gmin);
		for (j = 0; j < x_sz; j++) {
			for (i = 0; i < array_sz; i++) {
				Phi_gnj_t = array[j][i].Agnj_x * exp(creal(array[j][i].Lgnj_x) * t);
				Sum_Phi += Phi_gnj_t * array[j][i].Wgnj_x;
			}

			if (fabs(cimag(Sum_Phi)) > 0.0000000001) {
				fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
			}
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi) * power((double)P, -gmax));
	}

#if 0
	for (t = 100000; t < 1000000; t += 1000) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(creal(array[i].Lgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi));
	}

	for (t = 1000000; t < 100000000; t += 10000) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(creal(array[i].Lgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi));
	}


	for (t = 100000000; t < 10000000000; t += 100000) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp(creal(array[i].Lgnj_x) * t);
//			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
			Sum_Phi += Phi_gnj_t * array[i].Wgnj_x;
		}

		fprintf(output, "%.02f\t%g\n", t, creal(Sum_Phi));
	}
#endif


	fflush(output);
	free_pa_num(x0);
	close(fd);

	return 0;
}

