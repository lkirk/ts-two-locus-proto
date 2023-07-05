#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "prototype.h"

typedef int(summary_func)(tsk_size_t, const double *, tsk_size_t, double *, void *);

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

static tsk_flags_t
pick_norm_strategy(const char *func_name)
{
    if (!strcmp(func_name, "D")) {
        return TSK_TOTAL_WEIGHTED;
    } else if (!strcmp(func_name, "D2")) {
        return TSK_TOTAL_WEIGHTED;
    } else if (!strcmp(func_name, "r2")) {
        return TSK_HAP_WEIGHTED;
    } else if (!strcmp(func_name, "D_prime")) {
        return TSK_HAP_WEIGHTED;
    } else if (!strcmp(func_name, "r")) {
        return TSK_TOTAL_WEIGHTED;
    } else if (!strcmp(func_name, "Dz")) {
        return TSK_TOTAL_WEIGHTED;
    } else if (!strcmp(func_name, "pi2")) {
        return TSK_TOTAL_WEIGHTED;
    }
    return 0;
}

static bool
pick_polarisation(const char *func_name)
{
    if (!strcmp(func_name, "D")) {
        return true;
    } else if (!strcmp(func_name, "D2")) {
        return false;
    } else if (!strcmp(func_name, "r2")) {
        return false;
    } else if (!strcmp(func_name, "D_prime")) {
        return true;
    } else if (!strcmp(func_name, "r")) {
        return true;
    } else if (!strcmp(func_name, "Dz")) {
        return false;
    } else if (!strcmp(func_name, "pi2")) {
        return false;
    }
    return false;
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

    tsk_flags_t options = 0;
    options |= pick_norm_strategy(argv[1]);
    if (pick_polarisation(argv[1])) {
        options |= TSK_STAT_POLARISED;
    }

    tsk_treeseq_t ts;
    tsk_treeseq_load(&ts, filename, 0);

    /* two_locus_stat(&ts); */
    double *sample_weights = tsk_malloc(ts.num_samples * sizeof(*sample_weights));
    for (tsk_size_t i = 0; i < ts.num_samples; i++) {
        sample_weights[i] = 1;
    }

    const tsk_size_t sample_set_sizes[] = { ts.num_samples };
    const tsk_id_t set_indexes[] = { 0 };
    const tsk_size_t num_sample_sets = 1;
    tsk_id_t sample_sets[ts.num_samples];
    for (tsk_id_t i = 0; i < (tsk_id_t) ts.num_samples; i++) {
        sample_sets[i] = i;
    }

    /* double windows[] = { 0, ts.tables->sequence_length }; */

    double *result;
    sample_count_stat_params_t args = { .sample_sets = sample_sets,
        .num_sample_sets = num_sample_sets,
        .sample_set_sizes = sample_set_sizes,
        .set_indexes = set_indexes };

    int ret = two_site_general_stat(&ts, num_sample_sets, sample_weights,
        num_sample_sets, func, &args, 0, NULL, options, &result);

    if (ret != 0) {
        puts(tsk_strerror(ret));
        exit(1);
    }

    tsk_treeseq_free(&ts);
    tsk_safe_free(result);
    tsk_safe_free(sample_weights);
}
