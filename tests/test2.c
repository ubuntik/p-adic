#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <p-def.h>

#define G_MAX (2)
#define G_MIN (1)
#define XSZ (G_MAX - G_MIN)

int main(int argc, char **argv)
{
	pa_num** qs = NULL;
	pa_num *pgnum = NULL;
	PADIC_ERR err = ESUCCESS;
	int qs_sz = 0, i = 0, g = 0;

	printf("Test#2: Generate tree\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	for (g = G_MIN; g <= G_MAX; g++) {
		printf("Gamma = %d\n", g);

		qs_sz = (size_t)qspace_sz(G_MIN, g);
		printf("Number of balls = %d\n", qs_sz);

		qs = (pa_num **)malloc(qs_sz * sizeof(pa_num *));
		if (qs == NULL) {
			fprintf(stderr, "Cannot alloc memory\n");
			return EINVPNTR;
		}
		err = gen_quotient_space(qs, G_MIN, g);
		if (err != ESUCCESS) {
			fprintf(stderr, "Involid generating quotient space\n");
			return err;
		}

		for (i = 0; i < qs_sz; i++) {
			pgnum = (pa_num *)malloc(sizeof(pa_num));
			if (pgnum == NULL) {
				fprintf(stderr, "Cannot alloc memory\n");
				exit(-1);
			}
			err = p_gamma_pa_num(pgnum, qs[i], g);
			if (err != ESUCCESS) {
				fprintf(stderr, "Invalid mult on p gamma\n");
				exit(err);
			}
			printf("Number %d:\n", i);
			printf("Canonical view coefficients:\n");
			print_pa_num(pgnum);
			printf("User friendly view:\n");
			printf("%g\n", padic2double(pgnum));
			printf("===============================\n");
			free_pa_num(pgnum);
			free_pa_num(qs[i]);
		}
		printf("###############################\n");
		free(qs);
	}
	return 0;
}

