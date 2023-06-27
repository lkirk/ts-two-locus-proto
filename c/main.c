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

    two_locus_stat(&ts);
    tsk_treeseq_free(&ts);
}
