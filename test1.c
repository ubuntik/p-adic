#include "p-adic.h"
#define G_MAX 1
#define G_MIN (-1)
#define XSZ (G_MAX - G_MIN)

int main()
{
	pa_num** fs;
	int fs_sz, i;

	fs = gen_factor_space(G_MIN, G_MAX);
	fs_sz = (size_t)fspace_sz(G_MIN, G_MAX);

	printf("Test#1: Generate factor-space\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	for (i = 0; i < fs_sz; i++) {
		printf("Number %d:\n", i);
		printf("Canonical view coefficients:\n");
		print_pa_num(fs[i]);
		printf("%g\n", from_canonic_to_float(fs[i]));
		printf("===============================\n");
		free_pa_num(fs[i]);
	}

	free(fs);

	return 0;
}

