#include <stdlib.h>

#include "prototype.h"

int
main(int argc, char **argv)
{
    char *filename;
    if (argc != 2) {
        printf("usage: %s <input_tree>\n", argv[0]);
        exit(1);
    }
    filename = argv[1];

    tsk_treeseq_t ts;
    tsk_treeseq_load(&ts, filename, 0);

    /* two_locus_stat(&ts); */
    double *sample_weights = tsk_malloc(ts.num_samples * sizeof(*sample_weights));
    for (tsk_size_t i = 0; i < ts.num_samples; i++) {
        sample_weights[i] = 1;
    }
    double result;
    two_site_general_stat(&ts, 1, sample_weights, 1, NULL, NULL, 0, NULL, 0, &result);
    tsk_treeseq_free(&ts);
}
