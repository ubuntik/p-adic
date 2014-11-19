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
	complex double Wgnj_x;
	complex double Cgnj_x;
	complex double Bgnj_x;
	complex double Agnj_x;
	complex double Lgnj_x;
} Couchy_coeffs;

typedef struct Matrix_struct {
	int gamma;
	pa_num *n;
	int j;
	int xs_sz;
	pa_num **xs;
	Couchy_coeffs *coeff;
} Matrix_struct;

// Bgnj_x = Integral(rho_backward(x,y)*WAVEgnj(y)dy) by Qp
// Cgnj_x = Integral(rho_forward(x,y)*WAVEgnj(x)dy) by Qp
// Bgnj_x - Cgnj_x

double (*rho)(pa_num *x) = NULL;
double (*rho_backward)(pa_num *x, pa_num *y) = NULL;
double (*rho_forward)(pa_num *x, pa_num *y) = NULL;
double (*start_cond)(pa_num *pnum) = NULL;
int g_min = 0;
int g_max = -1;
int g_chy = -1;

PADIC_ERR do_for_j(Matrix_struct *mst)
{
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0, k = 0;
	pa_num **qs = NULL;
	pa_num *pa = NULL;
	pa_num *x = NULL;
	pa_num *n = NULL;
	int j = 0, gamma = 0;
	complex double W = I;
	double B_img = 0, B_rez = 0, B_fun = 0;
	double C_img = 0, C_rez = 0, C_fun = 0;
	double A_img = 0, A_rez = 0, A_fun = 0;
	complex double Bgnj_x = I, Cgnj_x = I, Agnj_x = I;

	mst->coeff = (Couchy_coeffs *)malloc(mst->xs_sz * sizeof(Couchy_coeffs));
	if (mst->coeff == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}

	n = mst->n;
	j = mst->j;
	gamma = mst->gamma;
	for (i = 0; i < mst->xs_sz; i++) {
		x = mst->xs[i];


		fprintf(stdout, "Get x = %g\n", padic2double(x));
		print_pa_num(x);
		fprintf(stdout, "gamma\tn\tj\tWgnj(x)\t\tBgnj(x)\t\tCgnj(x)\t\tL(x)\n");


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

		mst->coeff[i].Wgnj_x = wavelet(x, n, -gamma, j);
		W = mst->coeff[i].Wgnj_x;
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

			if ((rho_backward(x, pa) == INFINITY) || (rho_forward(x, pa) == INFINITY)) {
				fprintf(stderr, "AAAA!!! We are going to infinity!\n");
				return EINVPNTR;
			} else {
				B_fun = rho_backward(x, pa);
				C_fun = rho_forward(x, pa);
				A_fun = start_cond(pa);

				if (fabs(creal(W)) > ACCURACY) {
					C_rez += creal(W) * C_fun;
					C_img += cimag(W) * C_fun;
				}
				B_rez += creal(wavelet(pa, n, -gamma, j)) * B_fun;
				B_img += cimag(wavelet(pa, n, -gamma, j)) * B_fun;
				A_rez += creal(wavelet(pa, n, -gamma, j)) * A_fun;
				A_img += cimag(wavelet(pa, n, -gamma, j)) * A_fun;
			}
			free_pa_num(pa);
			free_pa_num(qs[k]);
		}

		Bgnj_x = B_rez * power((double)P, -g_max) + I * B_img * power((double)P, -g_max);
		Cgnj_x = C_rez * power((double)P, -g_max) + I * C_img * power((double)P, -g_max);
		Agnj_x = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);

		mst->coeff[i].Bgnj_x = Bgnj_x;
		mst->coeff[i].Cgnj_x = Cgnj_x;
		mst->coeff[i].Agnj_x = Agnj_x;

		if (fabs(creal(W)) > ACCURACY)
			mst->coeff[i].Lgnj_x = (Bgnj_x - Cgnj_x) / W;
		else
			mst->coeff[i].Lgnj_x = 0;


		fprintf(stdout, "%d\t%.04g\t%d\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\t%.04g\n",
			gamma, padic2double(n), j,
			PRNT_CMX(creal(mst->coeff[i].Wgnj_x)), PRNT_CMX(cimag(mst->coeff[i].Wgnj_x)),
			PRNT_CMX(creal(mst->coeff[i].Bgnj_x)), PRNT_CMX(cimag(mst->coeff[i].Bgnj_x)),
			PRNT_CMX(creal(mst->coeff[i].Cgnj_x)), PRNT_CMX(cimag(mst->coeff[i].Cgnj_x)),
			PRNT_CMX(creal(mst->coeff[i].Lgnj_x)), PRNT_CMX(cimag(mst->coeff[i].Lgnj_x)));


	}
	return err;
}

PADIC_ERR count_matrix(Matrix_struct *msts, int array_sz,
			double **matrix)
{
	int i = 0, j = 0, k = 0;
	PADIC_ERR err = ESUCCESS;
	double A_img = 0, A_rez = I;
	double B_img = I, B_rez = I;
	complex double A = 0;
	complex double B = 0;

	for (i = 0; i < array_sz; i++) {
		for (j = 0; j < array_sz; j++) {
			for (k = 0; k < msts[i].xs_sz; k++) {
				A_rez +=
					creal(wavelet(msts[i].xs[k], msts[i].n,
		//			creal(wavelet(msts[j].xs[k], msts[i].n,
					-msts[i].gamma, msts[i].j)
					* msts[j].coeff[k].Bgnj_x);
				A_img +=
					cimag(wavelet(msts[i].xs[k], msts[i].n,
		//			cimag(wavelet(msts[j].xs[k], msts[i].n,
					-msts[i].gamma, msts[i].j)
					* msts[j].coeff[k].Bgnj_x);
				B_rez +=
					creal(wavelet(msts[i].xs[k], msts[i].n,
		//			creal(wavelet(msts[j].xs[k], msts[i].n,
					-msts[i].gamma, msts[i].j)
					* msts[j].coeff[k].Cgnj_x);
				B_img +=
					cimag(wavelet(msts[i].xs[k], msts[i].n,
		//			cimag(wavelet(msts[j].xs[k], msts[i].n,
					-msts[i].gamma, msts[i].j)
					* msts[j].coeff[k].Cgnj_x);
			}

			A = A_rez * power((double)P, -g_max) - I * A_img * power((double)P, -g_max);
			B = B_rez * power((double)P, -g_max) - I * B_img * power((double)P, -g_max);

			if (fabs(cimag(A)) > ACCURACY)
				fprintf(stderr, ">>> A <<< Something wrong!!!\n");
			if (fabs(cimag(B)) > ACCURACY)
				fprintf(stderr, ">>> B <<< Something wrong!!!\n");

			matrix[i][j] = creal(A) - creal(B);
			fprintf(stdout, "%g\t", PRNT_CMX(matrix[i][j]));
		}
		fprintf(stdout, "\n");
	}
	return err;
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
	double **matrix = NULL;
	Matrix_struct *mst = NULL;
	int cnt = 0;

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

	mst = (Matrix_struct *)malloc(array_sz * sizeof(Matrix_struct));
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

	matrix = (double **)malloc(array_sz * sizeof(double *));
	if (matrix == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EMALLOC;
	}
	for (j = 0; j < array_sz; j++) {
		matrix[j] = (double *)malloc(array_sz * sizeof(double));
		if (matrix[j] == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EMALLOC;
		}
		mst[j].xs_sz = xs_sz;
		mst[j].xs = xs;
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
				mst[cnt].gamma = gamma;
				mst[cnt].n = ns[i];
				mst[cnt].j = j;

				err = do_for_j(&mst[cnt]);
				if (err != ESUCCESS)
					return err;
				cnt++;
			}
		}
	}

	if (cnt != array_sz) {
		fprintf(stderr, "Wrong calculaton on integrals: #%d, must be #%d\n", cnt, array_sz);
		return EINVPNTR;
	}

	err = count_matrix(mst, array_sz, matrix);
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to count Matrix\n");
		return err;
	}

	return 0;
}

