#include <stdio.h>
#include <stdlib.h>

#include "tskit.h"

int main(int argc, char **argv) {
  char *filename;
  if (argc != 2) {
    printf("usage: %s <input_tree>\n", argv[0]);
    exit(1);
  }
  filename = argv[1];

  tsk_treeseq_t ts;
  tsk_treeseq_load(&ts, filename, 0);
  // tsk_treeseq_print_state(&ts, stdout);

  tsk_id_t samples[ts.num_samples];
  for (int i = 0; i < ts.num_samples; i++) {
    samples[i] = i;
  }
  tsk_size_t sample_set_sizes = ts.num_samples;

  double pi;
  int ret = tsk_treeseq_diversity(&ts, 1, &sample_set_sizes, samples, 0, NULL,
                                  TSK_STAT_SITE, &pi);
  if (ret != 0) {
    puts("error!");
    exit(1);
  }
}
