#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "prototype.h"

static summary_func *
pick_summary_func(const char *func_name)
{
    if (!strcmp(func_name, "D")) {
        return D;
    } else if (!strcmp(func_name, "D2")) {
        return D2;
    } else if (!strcmp(func_name, "r2")) {
        return r2;
    } else if (!strcmp(func_name, "D_prime")) {
        return D_prime;
    } else if (!strcmp(func_name, "r")) {
        return r;
    } else if (!strcmp(func_name, "Dz")) {
        return Dz;
    } else if (!strcmp(func_name, "pi2")) {
        return pi2;
    }
    return NULL;
}

int
main(int argc, char **argv)
{
    char *filename;
    summary_func *func;
    bool print_weights = false;
    double *result;

    if (argc < 3 || argc > 4) {
        printf("usage: %s <summary_func> <input_tree> [print_weights]\n", argv[0]);
        exit(1);
    }

    filename = argv[2];
    if (access(filename, F_OK) != 0) {
        printf("ERROR: %s: File does not exist\n", filename);
        exit(1);
    }

    func = pick_summary_func(argv[1]);
    if (!func) {
        printf("error: unknown summary function '%s'\n", argv[1]);
        exit(1);
    }

    if (argc == 4) {
        if (strlen(argv[3]) != 1 || (argv[3][0] != '0' && argv[3][0] != '1')) {
            printf("error: acceptable values for print_weights are (0, 1), got %s\n",
                argv[3]);
            exit(1);
        }
        if (argv[3][0] == '1') {
            print_weights = true;
        }
    }

    tsk_size_t num_stat;
    return process_tree(func, filename, &result, &num_stat, print_weights);
    tsk_safe_free(result);
}
