#include "../src/p-adic.h"
#define G_MAX 1
#define G_MIN (-1)

int main()
{
	pa_num** qs = NULL;
	int qs_sz = 0, i = 0;
	PADIC_ERR err = ESUCCESS;

	qs_sz = (size_t)qspace_sz(G_MIN, G_MAX);
	qs = (pa_num **)malloc(qs_sz * sizeof(pa_num*));
	if (qs == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}
	err = gen_quotient_space(qs, G_MIN, G_MAX);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating quotient space\n");
		return err;
	}

	printf("Test#1: Generate factor-space\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	for (i = 0; i < qs_sz; i++) {
		printf("Number %d:\n", i);
		printf("Canonical view coefficients:\n");
		print_pa_num(qs[i]);
		printf("%g\n", from_canonic_to_double(qs[i]));
		printf("===============================\n");
		free_pa_num(qs[i]);
	}

	free(qs);

	return 0;
}

