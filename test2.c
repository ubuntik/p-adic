#include "p-adic.h"
#define G_MAX (2)
#define G_MIN (1)
#define XSZ (G_MAX - G_MIN)

int main()
{
	pa_num** fs;
	pa_num *pgnum;
	int fs_sz, i, g;

	printf("Test#2: Generate tree\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	for (g = G_MIN; g <= G_MAX; g++) {
		fs_sz = (size_t)fspace_sz(G_MIN, g);

		printf("Gamma = %d\n", g);
		printf("Number of balls = %d\n", fs_sz);

		fs = gen_factor_space(G_MIN, g);

		for (i = 0; i < fs_sz; i++) {
			pgnum = p_gamma_pa_num(fs[i], g);
			printf("Number %d:\n", i);
			printf("Canonical view coefficients:\n");
			print_pa_num(pgnum);
			printf("User friendly view:\n");
			printf("%g\n", from_canonic_to_float(pgnum));
			printf("===============================\n");
			free_pa_num(pgnum);
			free_pa_num(fs[i]);
		}
		printf("###############################\n");

		free(fs);

	}
	return 0;
}

