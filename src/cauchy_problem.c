#include "../src/cauchy_problem.h"
#include <errno.h>

typedef struct Couchy_coeffs {
	int gamma;
	pa_num *n;
	int j;
	complex Wgnj_x;
	double Cgnj_x;
	double Agnj_x;
	} Couchy_coeffs;

// how many integrals we have to count (j * gamma * n)
int array_sz = 0;
Couchy_coeffs *array = NULL;
int cnt = 0;

pa_num *x = NULL;
double (*rho)(pa_num *pnum) = NULL;
double (*start_cond)(pa_num *pnum) = NULL;
int g_min = 0;
int g_max = -1;

pa_num* special_poin_workarount(pa_num *spoint, pa_num *x)
{
	pa_num *ext_pa = NULL, *res = NULL;
	PADIC_ERR err = ESUCCESS;

	if (spoint == NULL) {
		fprintf(stderr, "Invalid argument\n");
		return spoint;
	}

	fprintf(stderr, "Warning!!! Special point has found!\n");
	print_pa_num(spoint);
	fprintf(stderr, "%g is special point!\n", from_canonic_to_double(spoint));

	ext_pa = (pa_num *)malloc(sizeof(pa_num));
	if (ext_pa == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return NULL;
	}
	err = __extend_number(ext_pa, spoint, spoint->g_min,
						spoint->g_max + 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid extending number\n");
		return NULL;
	}

	err = set_x_by_gamma(ext_pa, spoint->g_max + 1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid setting coeff\n");
		return NULL;
	}

	res = (pa_num *)malloc(sizeof(pa_num));
	if (res == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return NULL;
	}

	err = sub(res, ext_pa, x);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid subtraction\n");
		return NULL;
	}

	free_pa_num(spoint);
	free_pa_num(ext_pa);

	return res;
}

double __integral_B_x()
{
	double ret = 0;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *res= NULL;

	if (rho == NULL) {
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

	for (i = 0; i < qs_sz; i++) {
		/* get (P^gamma * x - n) */
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

		if ( rho(res) == INFINITY )
			res = special_poin_workarount(res, x);

		if (res == NULL) {
			fprintf(stderr, "Failed to workaround special point\n");
			return EINVPNTR;
		}

		ret += (double)rho(res);
		free_pa_num(res);
		free_pa_num(qs[i]);
		free_pa_num(pa);
	}
	free(qs);
	return ret * power((double)P, -g_max);
}

PADIC_ERR do_for_j(int gamma, pa_num *n, int j)
{
	PADIC_ERR err = ESUCCESS;
	double C_img = 0, C_rez = 0, C_fun = 0;
	double A_img = 0, A_rez = 0, A_fun = 0;
	int qs_sz = 0, i = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL, *res = NULL;
	complex Cgnj_x = I;
	complex Agnj_x = I;

	if (start_cond == NULL || rho == NULL) {
		fprintf(stderr, "Invalid pointer to function\n");
		return EINVPNTR;
	}

	if (x == NULL || n == NULL || array == NULL) {
		fprintf(stderr, "Invalid pointer to data\n");
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

		res = (pa_num *)malloc(sizeof(pa_num));
		if (res == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}
		err = __sign_sub(res, x, pa);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid subtraction\n");
			return err;
		}

		if ( rho(res) == INFINITY ) {
			C_fun = rho(special_poin_workarount(res, x));
			A_fun = start_cond(special_poin_workarount(res, x));
		} else {
			C_fun = rho(res);
			A_fun = start_cond(pa);
		}

		C_rez += creal(wavelet(res, n, -gamma, j)) * C_fun;
		C_img += cimag(wavelet(res, n, -gamma, j)) * C_fun;

		A_rez += creal(wavelet(pa, n, -gamma, j)) * A_fun;
		A_img += cimag(wavelet(pa, n, -gamma, j)) * A_fun;

		free_pa_num(pa);
		free_pa_num(res);
		free_pa_num(qs[i]);
	}
	free(qs);

	Cgnj_x = C_rez * power((double)P, -g_max) + I * A_img * power((double)P, -g_max);
	Agnj_x = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);

	if (cimag(Cgnj_x) > 0.000001) {
		fprintf(stderr, ">>> C <<< Something wrong!!!\n");
		return EINVCOEFF;
	}
	if (cimag(Agnj_x) > 0.000001) {
		fprintf(stderr, ">>> A <<< Something wrong!!!\n");
		return EINVCOEFF;
	}

	array[cnt].Cgnj_x = creal(Cgnj_x);
	array[cnt].Agnj_x = creal(Agnj_x);
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

	if (gamma < g_min || gamma > g_max) {
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
		array[cnt].Wgnj_x = wavelet(x, n, gamma, j);
		if (cimag(array[cnt].Wgnj_x) > 0.000001) {
			fprintf(stderr, ">>> W >>> Something wrong!!!");
			return EINVCOEFF;
		}

		err = do_for_j(gamma, n, j);
		if (err != ESUCCESS)
			return err;
		cnt++;
	}
	return err;
}

PADIC_ERR get_integrals()
{
	int gamma = 0;
	PADIC_ERR err = ESUCCESS;
	int ns_sz = 0;
	int i = 0;
	pa_num **ns = NULL;

	if (start_cond == NULL || rho == NULL) {
		fprintf(stderr, "Invalid pointer to function\n");
		return EINVPNTR;
	}
	if (x == NULL || array == NULL) {
		fprintf(stderr, "Invalid pointer to data\n");
		return EINVPNTR;
	}

	// count Cgnj(x) and Agnj_x integrals
	for (gamma = g_min; gamma < g_max; gamma++) {
		ns_sz = (size_t)qspace_sz(gamma, g_max);

		ns = (pa_num **)malloc(ns_sz * sizeof(pa_num*));
		if (ns == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}

		err = gen_quotient_space(ns, gamma, g_max);
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
		int gmin, int gmax)
{
	int fd = -1;
	char *srv_str = NULL;
	int i = 0;
	pa_num *x0 = NULL;
	PADIC_ERR err = ESUCCESS;
	double Sum_Phi = 0;
	double Phi_gnj_t = 0;
	float t = 0.01;
	double B_x = 0;
	int array_sz = 0;
	int gamma = 0;

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

	snprintf(srv_str, 17, "./res/%05d.xlsx\n", getpid());

#if 0
	if ((fd = open(srv_str, O_WRONLY | O_CREAT | O_EXCL, 0666)) < 0) {
		perror("Cannot open file");
		exit(errno);
	}
	bzero(srv_str, pathconf(".", _PC_PATH_MAX));
#endif

	for (gamma = gmin; gamma < gmax; gamma++)
		array_sz += qspace_sz(gamma, gmax) * (P - 1);

	array = (Couchy_coeffs *)malloc(array_sz * sizeof(Couchy_coeffs));
	if (array == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}

	// get x = 0
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

	g_min = gmin;
	g_max = gmax;
	x = x0;
	rho = rho0;
	start_cond = start_cond0;

	// count B(x) integral
	B_x = __integral_B_x();

	err = get_integrals();
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Integrals\n");
		exit(err);
	}

	if (cnt != array_sz) {
		fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
		exit(-1);
	}

	for (t = 0.01; t < 20; t += 0.01) {
		Sum_Phi = power(P, (double)gmax / 2);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = array[i].Agnj_x * exp((array[i].Cgnj_x - B_x) * t);
			Sum_Phi += Phi_gnj_t * creal(array[i].Wgnj_x);
		}

		printf("%f\t%g\n", t, Sum_Phi);
#if 0
		if (write(fd, (void *)srv_str, 15) < 0) {
			perror("Write to file failed");
			exit(errno);
		}
#endif
		bzero(srv_str, sizeof(*srv_str));
	}
	free_pa_num(x0);
	close(fd);
	free(srv_str);

	return 0;
}

