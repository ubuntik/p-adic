#include "../src/p-adic.h"
#define G_MAX 2
#define G_MIN (-1)
#define XSZ (G_MAX - G_MIN)

int main()
{
	pa_num *pn1 = NULL, *pn2 = NULL;
	pa_num *res = NULL, *psub = NULL, *padd= NULL;
	int i = 0;
	PADIC_ERR err = ESUCCESS;

	printf("Test#3\n\n");
	printf(">>> Checking p-adic arithmetic <<<\n\n");

	pn1 = (pa_num *)malloc(sizeof(pa_num));
	if (pn1 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn1, -2,-1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -2, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	printf("Take first p-adic number:\n");
	print_pa_num(pn1);
	printf("pa1 = %f\n\n", from_canonic_to_float(pn1));

	pn2 = (pa_num *)malloc(sizeof(pa_num));
	if (pn2 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn2, -2, -1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn2, -2, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn2, -1, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	printf("Take second p-adic number:\n");
	print_pa_num(pn2);
	printf("pa2 = %f\n\n", from_canonic_to_float(pn2));

	printf("Check their subtraction and addition:\n");

	psub = (pa_num *)malloc(sizeof(pa_num));
	if (psub == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = sub(psub, pn1, pn2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid subtraction\n");
		exit(err);
	}
	printf("pa1 - pa2 = {");
	print_pa_num(psub);
	printf("pa1 - pa2 = %f\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	padd = (pa_num *)malloc(sizeof(pa_num));
	if (padd == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = add(padd, pn1, pn2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid addition\n");
		exit(err);
	}
	printf("pa1 + pa2 = {");
	print_pa_num(padd);
	printf("pa1 + pa2 = %f\n\n", from_canonic_to_float(padd));
	free_pa_num(padd);

	psub = (pa_num *)malloc(sizeof(pa_num));
	if (psub == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = sub(psub, pn2, pn1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid subtraction\n");
		exit(err);
	}
	printf("pa2 - pa1 = {");
	print_pa_num(psub);
	printf("pa2 - pa1 = %f\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	free_pa_num(pn1);
	free_pa_num(pn2);

	pn1 = (pa_num *)malloc(sizeof(pa_num));
	if (pn1 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn1, -2,-1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -2, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	pn1->sign = NEG;

	printf("Take first p-adic number (negative):\n");
	print_pa_num(pn1);
	printf("pa1 = %f\n\n", from_canonic_to_float(pn1));

	pn2 = (pa_num *)malloc(sizeof(pa_num));
	if (pn2 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn2, -2, -1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn2, -2, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn2, -1, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	printf("Take second p-adic number:\n");
	print_pa_num(pn2);
	printf("pa2 = %f\n\n", from_canonic_to_float(pn2));

	printf("Check their subtraction and addition:\n");

	psub = (pa_num *)malloc(sizeof(pa_num));
	if (psub == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = sub(psub, pn1, pn2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid subtraction\n");
		exit(err);
	}
	printf("pa1 - pa2 = {");
	print_pa_num(psub);
	printf("pa1 - pa2 = %f\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	padd = (pa_num *)malloc(sizeof(pa_num));
	if (padd == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = add(padd, pn1, pn2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid addition\n");
		exit(err);
	}
	printf("pa1 + pa2 = {");
	print_pa_num(padd);
	printf("pa1 - pa2 = %f\n\n", from_canonic_to_float(padd));
	free_pa_num(padd);

	psub = (pa_num *)malloc(sizeof(pa_num));
	if (psub == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = sub(psub, pn2, pn1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid subtraction\n");
		exit(err);
	}
	printf("pa2 - pa1 = {");
	print_pa_num(psub);
	printf("pa2 - pa1 = %f\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	padd = (pa_num *)malloc(sizeof(pa_num));
	if (padd == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = add(padd, pn1, pn2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid addition\n");
		exit(err);
	}
	printf("pa2 + pa1 = {");
	print_pa_num(padd);
	printf("pa2 + pa1 = %f\n\n", from_canonic_to_float(padd));
	free_pa_num(padd);

	free_pa_num(pn1);
	free_pa_num(pn2);

	pn1 = (pa_num *)malloc(sizeof(pa_num));
	if (pn1 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn1, -2,-1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -2, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	printf("Take the p-adic number:\n");
	print_pa_num(pn1);
	printf("pa1 = %f\n\n", from_canonic_to_float(pn1));

	printf("Check multiplication of p-adic number and p^gamma:\n");
	printf("pa1 * p^2 = {");

	res = (pa_num *)malloc(sizeof(pa_num));
	if (res == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = p_gamma_pa_num(res, pn1, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Invalid mult on p gamma\n");
		exit(err);
	}
	print_pa_num(res);
	printf("pa1 * p^2 = %f\n\n", from_canonic_to_float(res));

	free_pa_num(res);
	free_pa_num(pn1);

	printf(">>> Checking Indicator function <<< \n\n");

	pn1 = (pa_num *)malloc(sizeof(pa_num));
	if (pn1 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn1, -2,-2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -2, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	printf("Take an abstract p-adic number:\n");
	print_pa_num(pn1);
	printf("x >> %f\n\n", from_canonic_to_float(pn1));

	pn2 = (pa_num *)malloc(sizeof(pa_num));
	if (pn2 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn2, -3, -3);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn2, -3, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	printf("Take a representative:\n");
	print_pa_num(pn2);
	printf("n >> %f\n\n", from_canonic_to_float(pn2));

	printf("Checking indicator for different gamma = [-2; 1]\n");
	for (i = -2; i < 2; i++) {
		printf("gamma = %d, indicator = %d\n", i, indicator(pn1, pn2, i));
	}

	free_pa_num(pn1);
	free_pa_num(pn2);

	pn1 = (pa_num *)malloc(sizeof(pa_num));
	if (pn1 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn1, -3,-1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -3, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -2, 0);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	err = set_x_by_gamma(pn1, -1, 1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}

	printf("Take an abstract p-adic number:\n");
	print_pa_num(pn1);
	printf("x >> %f\n\n", from_canonic_to_float(pn1));

	pn2 = (pa_num *)malloc(sizeof(pa_num));
	if (pn2 == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		exit(-1);
	}
	err = init_pa_num(pn2, -1, -1);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid init number\n");
		exit(err);
	}
	err = set_x_by_gamma(pn2, -1, 2);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid setting value\n");
		exit(err);
	}
	printf("Take a representative:\n");
	print_pa_num(pn2);
	printf("n >> %f\n\n", from_canonic_to_float(pn2));

	printf("Checking indicator for different gamma = [0; 5]\n");
	for (i = 0; i < 5; i++) {
		printf("gamma = %d, indicator = %d\n", i, indicator(pn1, pn2, i));
	}
	free_pa_num(pn1);
	free_pa_num(pn2);

	return 0;

}
