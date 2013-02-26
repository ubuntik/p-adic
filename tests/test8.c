#include "../src/p-adic.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#define G_MAX (0)
#define G_MIN (-2)
#define XSZ (G_MAX - G_MIN)

int main(int argc, char **argv)
{
	pa_tree* tree = NULL;
	PADIC_ERR err = ESUCCESS;

	printf("Test#8: Draw tree\n");
	printf("p = %d; Gamma_min = %d; Gamma_max = %d\n", P, G_MIN, G_MAX);

	tree = (pa_tree *)malloc(sizeof(pa_tree));
	if (tree == NULL) {
		fprintf(stderr, "Cannot alloc memory\n");
		return EINVPNTR;
	}

	err = get_pa_tree(tree, G_MIN, G_MAX);
	if (err != ESUCCESS) {
		fprintf(stderr, "Involid generating tree\n");
		return err;
	}

	printf("Tree size: %d\n", tree->tree_sz);

	err = print_tree(tree, "./tree.gv");
	if (err != ESUCCESS) {
		fprintf(stderr, "Failed to print tree\n");
		return err;
	}

	free_tree(tree);

	return 0;
}

