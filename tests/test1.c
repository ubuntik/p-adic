#include "../src/p-adic.h"
#define G_MAX 1
#define G_MIN (-1)

int main()
{
	pa_num** qs;
	int qs_sz, i;

	qs = gen_quotient_space(G_MIN, G_MAX);
	qs_sz = (size_t)qspace_sz(G_MIN, G_MAX);

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

