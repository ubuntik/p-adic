#include "../src/p-adic.h"
#define G_MAX 2
#define G_MIN (-1)
#define XSZ (G_MAX - G_MIN)

int main()
{
	pa_num *pn1, *pn2, *res, *psub, *padd;
	int i;

	printf("Test#3\n\n");
	printf(">>> Checking p-adic arithmetic <<<\n\n");

	pn1 = init_pa_num(-2,-1);
	pn1->x[0] = 2;
	pn1->x[1] = 1;

	printf("Take first p-adic number:\n");
	print_pa_num(pn1);
	printf("pa1 = %g\n\n", from_canonic_to_float(pn1));

	pn2 = init_pa_num(-2, -1);
	pn2->x[0] = 1;
	pn2->x[1] = 2;

	printf("Take second p-adic number:\n");
	print_pa_num(pn2);
	printf("pa2 = %g\n\n", from_canonic_to_float(pn2));

	printf("Check their subtraction and addition:\n");

	psub = sub(pn1, pn2);
	printf("pa1 - pa2 = {");
	print_pa_num(psub);
	printf("pa1 - pa2 = %g\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	padd = add(pn1, pn2);
	printf("pa1 + pa2 = {");
	print_pa_num(padd);
	printf("pa1 + pa2 = %g\n\n", from_canonic_to_float(padd));
	free_pa_num(padd);

	psub = sub(pn2, pn1);
	printf("pa2 - pa1 = {");
	print_pa_num(psub);
	printf("pa2 - pa1 = %g\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	free_pa_num(pn1);
	free_pa_num(pn2);

	pn1 = init_pa_num(-2,-1);
	pn1->x[0] = 2;
	pn1->x[1] = 1;
	pn1->sign = NEG;

	printf("Take first p-adic number (negative):\n");
	print_pa_num(pn1);
	printf("pa1 = %g\n\n", from_canonic_to_float(pn1));

	pn2 = init_pa_num(-2, -1);
	pn2->x[0] = 1;
	pn2->x[1] = 2;

	printf("Take second p-adic number:\n");
	print_pa_num(pn2);
	printf("pa2 = %g\n\n", from_canonic_to_float(pn2));

	printf("Check their subtraction and addition:\n");

	psub = sub(pn1, pn2);
	printf("pa1 - pa2 = {");
	print_pa_num(psub);
	printf("pa1 - pa2 = %g\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	padd = add(pn1, pn2);
	printf("pa1 + pa2 = {");
	print_pa_num(padd);
	printf("pa1 - pa2 = %g\n\n", from_canonic_to_float(padd));
	free_pa_num(padd);

	psub = sub(pn2, pn1);
	printf("pa2 - pa1 = {");
	print_pa_num(psub);
	printf("pa2 - pa1 = %g\n\n", from_canonic_to_float(psub));
	free_pa_num(psub);

	padd = add(pn2, pn1);
	printf("pa2 + pa1 = {");
	print_pa_num(padd);
	printf("pa2 + pa1 = %g\n\n", from_canonic_to_float(padd));
	free_pa_num(padd);

	free_pa_num(pn1);
	free_pa_num(pn2);

	pn1 = init_pa_num(-2,-1);
	pn1->x[0] = 2;
	pn1->x[1] = 1;

	printf("Take the p-adic number:\n");
	print_pa_num(pn1);
	printf("pa1 = %g\n\n", from_canonic_to_float(pn1));

	printf("Check multiplication of p-adic number and p^gamma:\n");
	printf("pa1 * p^2 = {");
	res = p_gamma_pa_num(pn1, 2);
	print_pa_num(res);
	printf("pa1 * p^2 = %g\n\n", from_canonic_to_float(res));

	free_pa_num(res);
	free_pa_num(pn1);

	printf(">>> Checking Indicator function <<< \n\n");

	pn1 = init_pa_num(-2,-2);
	pn1->x[0] = 1;
	printf("Take an abstract p-adic number:\n");
	print_pa_num(pn1);
	printf("x >> %g\n\n", from_canonic_to_float(pn1));

	pn2 = init_pa_num(-3, -3);
	pn2->x[0] = 1;
	printf("Take a representative:\n");
	print_pa_num(pn2);
	printf("n >> %g\n\n", from_canonic_to_float(pn2));

	printf("Checking indicator for different gamma = [-2; 1]\n");
	for (i = -2; i < 2; i++) {
		printf("gamma = %d, indicator = %d\n", i, indicator(pn1, pn2, i));
	}

	free_pa_num(pn1);
	free_pa_num(pn2);

	pn1 = init_pa_num(-3,-1);
	pn1->x[0] = 2;
	pn1->x[1] = 0;
	pn1->x[2] = 1;
	printf("Take an abstract p-adic number:\n");
	print_pa_num(pn1);
	printf("x >> %g\n\n", from_canonic_to_float(pn1));

	pn2 = init_pa_num(-1, -1);
	pn2->x[0] = 2;
	printf("Take a representative:\n");
	print_pa_num(pn2);
	printf("n >> %g\n\n", from_canonic_to_float(pn2));

	printf("Checking indicator for different gamma = [0; 5]\n");
	for (i = 0; i < 5; i++) {
		printf("gamma = %d, indicator = %d\n", i, indicator(pn1, pn2, i));
	}
	free_pa_num(pn1);
	free_pa_num(pn2);

	return 0;

}
