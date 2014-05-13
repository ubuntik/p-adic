#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <p-analysis.h>

#define G_MAX 1
#define G_MIN (-2)


int main()
{
	pa_num *pa;
	pa_num *x, *n;
	int gamma = -1;
	int ind;
	int i, j;
	complex double res;
	PADIC_ERR err = ESUCCESS;
#if 0
	printf("Test#3: Wavelet functions\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	pa = (pa_num *)malloc(sizeof(pa_num));
	if (pa == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pa, G_MIN, G_MAX);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pa, G_MIN, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pa, G_MIN + 1, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pa, G_MAX - 1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pa, G_MAX, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	print_pa_num(pa);
	printf("Number: %g\n", padic2double(pa));

	res = character(pa);
	printf("Character: %g + i%g\n", creal(res), cimag(res));
	free_pa_num(pa);

	for (i = 0; i < P; i++) {
		x = (pa_num *)malloc(sizeof(pa_num));
		if (x == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			exit(-1);
		}
		err = init_pa_num(x, -1, 0);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid init number\n");
			exit(err);
		}
		err = set_x_by_gamma(x, -1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid setting value\n");
			exit(err);
		}
		err = set_x_by_gamma(x, 0, i);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid setting value\n");
			exit(err);
		}
		print_pa_num(x);
		printf("x >> %g\n", padic2double(x));

		n = (pa_num *)malloc(sizeof(pa_num));
		if (n == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			exit(-1);
		}
		err = init_pa_num(n, -1, -1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid init number\n");
			exit(err);
		}
		err = set_x_by_gamma(n, -1, 1);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid setting value\n");
			exit(err);
		}
		print_pa_num(n);
		printf("n >> %g\n", padic2double(n));

		gamma = 0;

		ind = indicator(x, n, gamma);
		printf("gamma = %d, indicator = %d\n", gamma, ind);

		for (j = 1; j < P; j++) {
			res = wavelet(x, n, gamma, j);
			printf("(j = %d) >> %g + i%g\n", j, creal(res), cimag(res));
		}
		printf("\n");
		free_pa_num(x);
		free_pa_num(n);
	}
#endif

	pa = (pa_num *)malloc(sizeof(pa_num));
	if (pa == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pa, 0, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}

	err = set_x_by_gamma(pa, 0, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	err = set_x_by_gamma(pa, 1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	err = set_x_by_gamma(pa, 2, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}


	printf(">>>> %g + i %g\n", character(pa));

	return 0;
}

