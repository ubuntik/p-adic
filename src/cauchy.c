#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include <cauchy.h>

#define TIME 100
#define ACCURACY 0.0000000001
#define PRNT_CMX(x) (fabs(x) > ACCURACY ? x : 0)
#define PP(a) power((double)P, (double)a)

typedef struct Matrix_struct {
	pa_num *x;
	int gamma;
	pa_num *n;
	int j;
	complex double Wgnj_x;
	complex double Cgnj_x;
	complex double Bgnj_x;
	complex double Agnj_x;
	complex double Lgnj_x;
} Matrix_struct;

// Bgnj_x = Integral(rho_backward(x,y)*WAVEgnj(y)dy) by Qp
// Cgnj_x = Integral(rho_forward(x,y)*WAVEgnj(x)dy) by Qp
// Bgnj_x - Cgnj_x

int xs_sz = 0;
int array_sz = 0;
double (*rho)(pa_num *x) = NULL;
double (*rho_backward)(pa_num *x, pa_num *y) = NULL;
double (*rho_forward)(pa_num *x, pa_num *y) = NULL;
double (*start_cond)(pa_num *pnum) = NULL;
int g_min = 0;
int g_max = -1;
int g_chy = -1;

extern void dgeevx( char* blanc, char* jobvl, char* jobvr, char* sense,
			int* n, double* a, int* lda, double* wr, double* wi,
			double* vl, int* ldvl, double* vr, int* ldvr,
			int ilo, int ihi, double *scale, double* abnrm, double* rconde, double* rcondv,
			double* work, int* lwork, int* iwork, int* info );
extern void dgeev( char* jobvl, char* jobvr,
			int* n, double* a, int* lda, double* wr, double* wi,
			double* vl, int* ldvl, double* vr, int* ldvr,
			double* work, int* lwork, int* info );
extern void print_eigenvalues( char* desc, int n, double* wr, double* wi );
extern void print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv );

int find_eingens(int array_sz, double* a, double* wr, double *wi, double *vr) {
	int n = array_sz, lda = array_sz, ldvl = array_sz, ldvr = array_sz, info, lwork;
	double wkopt;
	double* work;
	double vl[ldvl*n];
/*	double scale[n], rconde[n], rcondv[n];
	double abnrm = 0;
	int ilo = 1, ihi = n;
	int iwork[2*n - 2];
*/

	lwork = -1;

	dgeev("Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
			&wkopt, &lwork, &info);

//	dgeevx("B", "V", "V", "B", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
//			ilo, ihi, scale, &abnrm, rconde, rcondv, &wkopt, &lwork, iwork, &info);

	lwork = (int)wkopt;
	work = (double*)calloc(lwork, sizeof(double));
//	work = (double*)calloc(2 * n * (n + 6), sizeof(double));

	dgeev("Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
			work, &lwork, &info);

//	dgeevx("B", "V", "V", "B", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
//			ilo, ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork, &info);

	if(info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		return -1;
	}

	free((void*)work);
	return 0;
}

void print_eigenvalues(char* desc, int n, double* wr, double* wi) {
	int j;
	printf("\n %s\n", desc);
	for(j = 0; j < n; j++) {
			fprintf(stdout, " %3.4f", wr[j] );
	}
	printf( "\n" );
}

void print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < n; i++ ) {
		for( j = 0; j < n; j++ ) {
			printf( " %3.4f", v[i+j*ldv] );
		}
		printf( "\n" );
	}
}

PADIC_ERR do_for_j(Matrix_struct *mst, int gamma, pa_num *n, int j, pa_num **xs, int xs_sz)
{
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0, k = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
	double B_img = 0, B_rez = 0, B_fun = 0;
	double C_img = 0, C_rez = 0, C_fun = 0;
	double A_img = 0, A_rez = 0, A_fun = 0;
	complex double Bgnj_x = I, Cgnj_x = I, Agnj_x = I;

	for (i = 0; i < xs_sz; i++) {
		mst[i].n = n;
		mst[i].j = j;
		mst[i].gamma = gamma;
		mst[i].x = xs[i];
		mst[i].Wgnj_x = wavelet(mst[i].x, n, -gamma, j);

		B_img = 0, B_rez = 0, B_fun = 0;
		C_img = 0, C_rez = 0, C_fun = 0;
		A_img = 0, A_rez = 0, A_fun = 0;

/*
		printf(">> x = %g, n = %g, g = %d, W = %g %g\n",
			padic2double(mst[i].x), padic2double(mst[i].n), mst[i].gamma,
			PRNT_CMX(creal(mst[i].Wgnj_x)), PRNT_CMX(cimag(mst[i].Wgnj_x)));
*/

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

	//	W = mst[i].Wgnj_x;
		for (k = 0; k < qs_sz; k++) {
			pa = (pa_num *)malloc(sizeof(pa_num));
			if (pa == NULL) {
				fprintf(stderr, "Cannot alloc memory\n");
				return EMALLOC;
			}
			err = p_gamma_pa_num(pa, qs[k], g_max);
			if (err != ESUCCESS) {
				fprintf(stderr, "Invalid multiplication on p-gamma\n");
				return err;
			}

			if ((rho_backward(mst[i].x, pa) == INFINITY) || (rho_forward(mst[i].x, pa) == INFINITY)) {
				fprintf(stderr, "AAAA!!! We are going to infinity!\n%g %g\n",
						padic2double(mst[i].x), padic2double(pa));
				fprintf(stderr, "rho_backward = %g, rho_forward = %g\n",
						rho_backward(mst[i].x, pa), rho_forward(mst[i].x, pa));

				return EINVPNTR;
			} else {
				B_fun = rho_backward(mst[i].x, pa);
				C_fun = rho_forward(mst[i].x, pa);
				A_fun = start_cond(pa);

				if (fabs(creal(mst[i].Wgnj_x)) > ACCURACY) {
					C_rez += creal(mst[i].Wgnj_x) * C_fun;
					C_img += cimag(mst[i].Wgnj_x) * C_fun;
				}
				B_rez += creal(wavelet(pa, n, -gamma, j)) * B_fun;
				B_img += cimag(wavelet(pa, n, -gamma, j)) * B_fun;
				A_rez += creal(wavelet(pa, n, -gamma, j)) * A_fun;
				A_img += cimag(wavelet(pa, n, -gamma, j)) * A_fun;
			}

/*
			fprintf(stderr, "xi = %g >> B: %g %g\tC: %g %g\n",
					padic2double(pa),
					PRNT_CMX(B_rez), PRNT_CMX(B_img),
					PRNT_CMX(C_rez), PRNT_CMX(C_img));
*/

			free_pa_num(pa);
			free_pa_num(qs[k]);
		}

		Bgnj_x = B_rez * power((double)P, -g_max) + I * B_img * power((double)P, -g_max);
		Cgnj_x = C_rez * power((double)P, -g_max) + I * C_img * power((double)P, -g_max);
		Agnj_x = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);

		mst[i].Bgnj_x = Bgnj_x;
		mst[i].Cgnj_x = Cgnj_x;
		mst[i].Agnj_x = Agnj_x;

		if (fabs(creal(mst[i].Wgnj_x)) > ACCURACY)
			mst[i].Lgnj_x = (Bgnj_x - Cgnj_x) / mst[i].Wgnj_x;
		else
			mst[i].Lgnj_x = 0;

	}
	return err;
}

PADIC_ERR count_matrix(Matrix_struct **mst, int array_sz, double *matrix, int xs_sz)
{
	int i = 0, j = 0, k = 0;
	PADIC_ERR err = ESUCCESS;
	double A_img = 0, A_rez = 0;
	double B_img = 0, B_rez = 0;
	complex double A = 0;
	complex double B = 0;

	fprintf(stdout, "\n\n");
	for (i = 0; i < array_sz; i++) {
		for (j = 0; j < array_sz; j++) {
			for (k = 0; k < xs_sz; k++) {
				A_rez += creal(conj(mst[i][k].Wgnj_x) * mst[j][k].Bgnj_x);
				A_img += cimag(conj(mst[i][k].Wgnj_x) * mst[j][k].Bgnj_x);
				B_rez += creal(conj(mst[i][k].Wgnj_x) * mst[j][k].Cgnj_x);
				B_img += cimag(conj(mst[i][k].Wgnj_x) * mst[j][k].Cgnj_x);
			}

			A = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);
			B = B_rez * power((double)P, -g_max) - I * B_img * power((double)P, -g_max);

			if (fabs(cimag(A)) > ACCURACY)
				fprintf(stderr, ">>> A <<< Something wrong!!!\n");
			if (fabs(cimag(B)) > ACCURACY)
				fprintf(stderr, ">>> B <<< Something wrong!!!\n");

			matrix[j*array_sz+i] = creal(A) - creal(B);
			fprintf(stdout, "%.3g\t", PRNT_CMX(matrix[j*array_sz+i]));
			A_img = 0, A_rez = 0;
			B_img = 0, B_rez = 0;
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

PADIC_ERR solve_problem(
		double (*rho_bw)(pa_num *x, pa_num *y),
		double (*rho_fw)(pa_num *x, pa_num *y),
		double (*start_cond0)(pa_num *pnum),
		int gmin, int gmax, int gchy)
{
	int i = 0, j = 0, gamma = 0;
	PADIC_ERR err = ESUCCESS;
	int array_sz = 0, ns_sz = 0, xs_sz = 0;
	pa_num **ns = NULL;
	pa_num **xs = NULL;
	double *matrix = NULL;
	Matrix_struct **mst = NULL;
	int cnt = 0;
	FILE *output, *lnout;
	complex double Sum_Phi = 0, Phi_gnj_t = 0;
    double t = 0, step = 0.001, max = 1;

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

	for (j = gmin; j <= gchy; j++)
		array_sz += qspace_sz(j + 1, gchy + 1) * (P - 1);

	mst = (Matrix_struct **)malloc(array_sz * sizeof(Matrix_struct*));
	if (mst == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	xs_sz = (size_t)qspace_sz(gmin, gmax);
	xs = (pa_num **)malloc(xs_sz * sizeof(pa_num*));
	if (xs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}


	err = gen_quotient_space_gamma(xs, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating qspace\n");
		return err;
	}

	matrix = (double *)calloc(array_sz*array_sz, sizeof(double));
	if (matrix == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}
	for (j = 0; j < array_sz; j++) {
		mst[j] = (Matrix_struct *)calloc(xs_sz, sizeof(Matrix_struct));
		if (mst[j] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}
	}
	g_min = gmin;
	g_max = gmax;
	g_chy = gchy;
	rho_backward = rho_bw;
	rho_forward = rho_fw;
	start_cond = start_cond0;

	cnt = 0;
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
			for (j = 1; j < P; j++) {
				err = do_for_j(mst[cnt], gamma, ns[i], j, xs, xs_sz);
				if (err != ESUCCESS)
					return err;
				cnt++;
			}
		}
	}


	fprintf(stdout, "x\tgamma\tn\tj\tWgnj(x)\t\tBgnj(x)\t\tCgnj(x)\t\tL(x)\t\tA(x)\n");
	for (i = 0; i < xs_sz; i++) {
		for (j = 0; j < array_sz; j++) {
			fprintf(stdout, "%.04g\t%d\t%.04g\t%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\n",
				padic2double(mst[j][i].x), mst[j][i].gamma, padic2double(mst[j][i].n), mst[j][i].j,
				PRNT_CMX(creal(mst[j][i].Wgnj_x)), PRNT_CMX(cimag(mst[j][i].Wgnj_x)),
				PRNT_CMX(creal(mst[j][i].Bgnj_x)), PRNT_CMX(cimag(mst[j][i].Bgnj_x)),
				PRNT_CMX(creal(mst[j][i].Cgnj_x)), PRNT_CMX(cimag(mst[j][i].Cgnj_x)),
				PRNT_CMX(creal(mst[j][i].Lgnj_x)), PRNT_CMX(cimag(mst[j][i].Lgnj_x)),
				PRNT_CMX(creal(mst[j][i].Agnj_x)), PRNT_CMX(cimag(mst[j][i].Agnj_x)));
		}
	}


	if (cnt != array_sz) {
		fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
		return EINVPNTR;
	}

	err = count_matrix(mst, array_sz, matrix, xs_sz);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Matrix\n");
		return err;
	}

	double wr[array_sz], wi[array_sz], vr[array_sz*array_sz];
	find_eingens(array_sz, matrix, wr, wi, vr);

/*
	double my_vr[] = {
			1,0,0,0,0,0,0,0,
			0,0,1,0,0,0,0,0,
			0,1,0,0,0,0,0,0,
			0,0,0,1,0,0,0,0,
			0,0,0,0,1,0,0,0,
			0,0,0,0,0,1,0,0,
			0,0,0,0,0,0,1,0,
			0,0,0,0,0,0,0,1};

	memcpy(&vr[0], &my_vr[0], array_sz*array_sz*sizeof(double));
*/

	print_eigenvalues("Eigenvalues", array_sz, &wr[0], &wi[0]);
	print_eigenvectors("Right eigenvectors", array_sz, &wi[0], &vr[0], array_sz);

	if ((output = __prepare_file("pl")) == NULL) {
		fprintf(stderr, "Cannot open file\n");
		return EINVPNTR;
	}

	if ((lnout = __prepare_file("ln")) == NULL) {
		fprintf(stderr, "Cannot open file\n");
		return EINVPNTR;
	}

	step = 0.1;
	//max = 10000;
	max = 1;
	int idx_x = 0;

	for (t = step; t <= max; t += step) {
		Sum_Phi = power(P, (double)gmin);
		for (i = 0; i < array_sz; i++) {
			Phi_gnj_t = 0;
			for (j = 0; j < array_sz; j++) {
				if(fabs(wi[j]) > ACCURACY)
					continue;
				Phi_gnj_t += vr[i+j*array_sz] * exp(wr[j] * t);
			}
			Sum_Phi += Phi_gnj_t * mst[i][idx_x].Agnj_x * mst[i][idx_x].Wgnj_x;
		}

		if (fabs(cimag(Sum_Phi)) > ACCURACY)
			fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
		fprintf(output, "%.03g\t%g\n", t, creal(Sum_Phi));
		fprintf(lnout, "%.03g\t%g\n", log(t), log(creal(Sum_Phi)));
		if ((max - t < step) && (max < TIME)) {
			step *= 10;
			max *= 10;
		}
	}


	return 0;
}

PADIC_ERR count_S_t(
		double (*rho_bw)(pa_num *x, pa_num *y),
		double (*rho_fw)(pa_num *x, pa_num *y),
		double (*start_cond0)(pa_num *pnum),
		int gmin, int gmax, int gchy,
		pa_num *ini_n, int ini_gamma)

{
	int i = 0, j = 0, k = 0, gamma = 0;
	PADIC_ERR err = ESUCCESS;
	complex double Sum_Phi = 0, Phi_gnj_t = 0;
	int array_sz = 0, ns_sz = 0, xs_sz = 0;
	FILE *lnout;
	double t = 0, step = 0.001, max = 1;
	double *matrix = NULL;
	Matrix_struct **mst = NULL;
	int cnt = 0;
	pa_num **ns = NULL;
	pa_num **xs = NULL;

(void)ini_n;

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

	for (j = gmin; j <= gchy; j++)
		array_sz += qspace_sz(j + 1, gchy + 1) * (P - 1);

	mst = (Matrix_struct **)malloc(array_sz * sizeof(Matrix_struct*));
	if (mst == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	xs_sz = (size_t)qspace_sz(gmin, gmax);
	xs = (pa_num **)malloc(xs_sz * sizeof(pa_num*));
	if (xs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}


	err = gen_quotient_space_gamma(xs, gmin, gmax);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating qspace\n");
		return err;
	}

	matrix = (double *)calloc(array_sz*array_sz, sizeof(double));
	if (matrix == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}
	for (j = 0; j < array_sz; j++) {
		mst[j] = (Matrix_struct *)calloc(xs_sz, sizeof(Matrix_struct));
		if (mst[j] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}
	}
	g_min = gmin;
	g_max = gmax;
	g_chy = gchy;
	rho_backward = rho_bw;
	rho_forward = rho_fw;
	start_cond = start_cond0;

	cnt = 0;
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
			for (j = 1; j < P; j++) {
				err = do_for_j(mst[cnt], gamma, ns[i], j, xs, xs_sz);
				if (err != ESUCCESS)
					return err;
				cnt++;
			}
		}
	}

	fprintf(stdout, "x\tgamma\tn\tj\tWgnj(x)\t\tBgnj(x)\t\tCgnj(x)\t\tL(x)\n");
	for (i = 0; i < xs_sz; i++) {
		for (j = 0; j < array_sz; j++) {
			fprintf(stdout, "%.04g\t%d\t%.04g\t%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\n",
				padic2double(mst[j][i].x), mst[j][i].gamma, padic2double(mst[j][i].n), mst[j][i].j,
				PRNT_CMX(creal(mst[j][i].Wgnj_x)), PRNT_CMX(cimag(mst[j][i].Wgnj_x)),
				PRNT_CMX(creal(mst[j][i].Bgnj_x)), PRNT_CMX(cimag(mst[j][i].Bgnj_x)),
				PRNT_CMX(creal(mst[j][i].Cgnj_x)), PRNT_CMX(cimag(mst[j][i].Cgnj_x)),
				PRNT_CMX(creal(mst[j][i].Lgnj_x)), PRNT_CMX(cimag(mst[j][i].Lgnj_x)));
		}
	}

	if (cnt != array_sz) {
		fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
		return EINVPNTR;
	}

	err = count_matrix(mst, array_sz, matrix, xs_sz);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Matrix\n");
		return err;
	}

	double wr[array_sz], wi[array_sz], vr[array_sz*array_sz];
	find_eingens(array_sz, matrix, wr, wi, vr);

	print_eigenvalues("Eigenvalues", array_sz, &wr[0], &wi[0]);
	print_eigenvectors("Right eigenvectors", array_sz, &wi[0], &vr[0], array_sz);

	if ((lnout = __prepare_file("ln")) == NULL) {
		fprintf(stderr, "Cannot open file\n");
		return EINVPNTR;
	}

	int s_xs_sz = (size_t)qspace_sz(gchy, gmax);

	step = 0.001;
	max = 1;
	for (t = step; t <= max; t += step) {
		Sum_Phi = s_xs_sz * power(P, (double)gmin);
		for (k = 0; k < s_xs_sz; k++) {
			for (i = 0; i < array_sz; i++) {
				Phi_gnj_t = 0;
				for (j = 0; j < array_sz; j++) {
					Phi_gnj_t += vr[i+j*array_sz] * exp(wr[j] * t);
				}
				Sum_Phi += Phi_gnj_t * mst[i][k].Agnj_x * mst[i][k].Wgnj_x;
			}
			if (fabs(cimag(Sum_Phi)) > ACCURACY)
				fprintf(stderr, ">>> Sum_Phi <<< Something wrong!!!\n");
		}

		fprintf(lnout, "%.03f\t%g\n", log(t), log(creal(Sum_Phi) * power((double)P, -gmax)));

		if ((max - t < step) && (max < TIME)) {
			step *= 10;
			max *= 10;
		}
	}

	fflush(lnout);
	return 0;
}


