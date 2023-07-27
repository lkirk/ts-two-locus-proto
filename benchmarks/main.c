#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>

#include <tskit.h>

static void
usage(char *progname)
{
    fprintf(
        stderr, "usage: %s <subcommand> -t <input_tree> -s <summary_func>\n", progname);
    exit(EXIT_FAILURE);
}

static void
die_usage(char *progname, const char *err, ...)
{
    va_list args;
    va_start(args, err);
    vfprintf(stderr, err, args);
    va_end(args);
    usage(progname);
}

static void
die(const char *err, ...)
{
    va_list args;
    va_start(args, err);
    vfprintf(stderr, err, args);
    va_end(args);
    exit(EXIT_FAILURE);
}

typedef int summary_func(const tsk_treeseq_t *self, tsk_size_t num_sample_sets,
    const tsk_size_t *sample_set_sizes, const tsk_id_t *sample_sets,
    tsk_size_t num_left_windows, const double *left_windows,
    tsk_size_t num_right_windows, const double *right_windows, tsk_flags_t options,
    tsk_size_t *result_size, double **result);

static summary_func *
pick_summary_func(const char *func_name)
{
    if (!strcmp(func_name, "D")) {
        return tsk_treeseq_D;
    } else if (!strcmp(func_name, "D2")) {
        return tsk_treeseq_D2;
    } else if (!strcmp(func_name, "r2")) {
        return tsk_treeseq_r2;
    } else if (!strcmp(func_name, "D_prime")) {
        return tsk_treeseq_D_prime;
    } else if (!strcmp(func_name, "r")) {
        return tsk_treeseq_r;
    } else if (!strcmp(func_name, "Dz")) {
        return tsk_treeseq_Dz;
    } else if (!strcmp(func_name, "pi2")) {
        return tsk_treeseq_pi2;
    }
    return NULL;
}

#define NEW_SUBCOMMAND (1UL << 1)
#define OLD_SUBCOMMAND (1UL << 2)

static void
parse_args(int argc, char **argv, uint32_t *subcommand, char **tree_filename,
    char **summary_func_name)
{
    int flag;
    bool posix_set = false;

    // make getopt behave!
    if (!!getenv("POSIXLY_CORRECT")) {
        posix_set = true;
    } else {
        setenv("POSIXLY_CORRECT", "", 0);
    }
    while (optind < argc) {
        if ((flag = getopt(argc, argv, "t:s:h")) != -1) {
            switch (flag) {
                case 't':
                    *tree_filename = optarg;
                    break;
                case 's':
                    *summary_func_name = optarg;
                    break;
                case 'h':
                    usage(argv[0]);
                    break;
                case '?':
                    die_usage(argv[0], "ERROR: Unrecognized option: %c\n", optopt);
                    break;
            }
        } else {
            if (*subcommand != 0) {
                die_usage(argv[0],
                    "ERROR: subcommand already set to: '%s'"
                    ", subcommands are mutually exclusive\n",
                    *subcommand == NEW_SUBCOMMAND ? "new" : "old");
            }
            if (!strcmp(argv[optind], "new")) {
                *subcommand = NEW_SUBCOMMAND;
            } else if (!strcmp(argv[optind], "old")) {
                *subcommand = OLD_SUBCOMMAND;
            } else {
                die_usage(argv[0], "ERROR: Unrecognized option: %s\n", argv[optind]);
            }
            optind++;
        }
    }

    switch (*subcommand) {
        case 0:
            die_usage(argv[0], "ERROR: one subcommand 'new' or 'old' must be set\n");
            break;
        case NEW_SUBCOMMAND:
            if (!summary_func_name) {
                die_usage(
                    argv[0], "ERROR: summary function (-s) is a required argument\n");
            }
    }

    if (!tree_filename) {
        die_usage(argv[0], "ERROR: tree filename (-t) is a required argument\n");
    }

    if (!posix_set) {
        unsetenv("POSIXLY_CORRECT");
    }
}

int
main(int argc, char **argv)
{
    char *tree_filename = NULL;
    char *summary_func_name = NULL;
    uint32_t subcommand = 0;

    parse_args(argc, argv, &subcommand, &tree_filename, &summary_func_name);

    fprintf(stdout, "Subcommand: '%s'\n", subcommand == NEW_SUBCOMMAND ? "new" : "old");

    summary_func *func;
    if (subcommand == NEW_SUBCOMMAND) {
        if ((func = pick_summary_func(summary_func_name)) == NULL) {
            die("ERROR: unknown summary func name: %s\n", summary_func_name);
        };
        fprintf(stdout, "Summary Function: '%s'\n", summary_func_name);
    }

    if (access(tree_filename, F_OK) != 0) {
        die("ERROR: %s: File does not exist\n", tree_filename);
    }

    fprintf(stdout, "Tree Sequence File: '%s'\n", tree_filename);

    tsk_treeseq_t ts;
    tsk_treeseq_load(&ts, tree_filename, 0);

    tsk_size_t sample_set_sizes[1] = { ts.num_samples };
    tsk_size_t num_sample_sets = 1;
    /* tsk_id_t sample_sets[ts.num_samples]; */
    tsk_id_t *sample_sets = tsk_malloc(ts.num_samples * sizeof(*sample_sets));
    for (tsk_size_t i = 0; i < ts.num_samples; i++) {
        sample_sets[i] = (tsk_id_t) i;
    }

    int ret = 0;
    double *result;
    tsk_size_t result_size;
    tsk_size_t expected_result_size
        = (ts.tables->sites.num_rows * (ts.tables->sites.num_rows + 1) / 2UL);

    switch (subcommand) {
        case NEW_SUBCOMMAND:
            ret = tsk_treeseq_D(&ts, num_sample_sets, sample_set_sizes, sample_sets, 0,
                NULL, 0, NULL, 0, &result_size, &result);
            break;
        case OLD_SUBCOMMAND:
            result = tsk_calloc(expected_result_size, sizeof(*result));
            result_size = 0;
            for (tsk_size_t site = 0; site < ts.tables->sites.num_rows; site++) {
                tsk_ld_calc_t ld_calc;
                tsk_ld_calc_init(&ld_calc, &ts);
                result += result_size;
                ret = tsk_ld_calc_get_r2_array(&ld_calc, (tsk_id_t) site,
                    TSK_DIR_FORWARD, UINT64_MAX, FLT_MAX, result, &result_size);
            }
            break;
        default:
            die("ERROR: unknown subcommand %d\n", subcommand);
    }

    tsk_treeseq_free(&ts);

    if (ret != 0) {
        puts(tsk_strerror(ret));
    }

    return ret;
}
