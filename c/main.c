#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

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
    if (argc != 3) {
        printf("usage: %s <summary_func> <input_tree>\n", argv[0]);
        exit(1);
    }
    filename = argv[2];

    summary_func *func = pick_summary_func(argv[1]);
    if (!func) {
        printf("error: unknown summary function '%s'\n", argv[1]);
        exit(1);
    }

    return process_tree(func, filename);
}
