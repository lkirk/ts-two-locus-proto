#include <stdio.h>
#include <stdlib.h>

#include "tskit.h"
#include "tskit/stats.h"

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

  /* int num_sample_sets = 4; */
  /* tsk_id_t sample_sets[] = {0, 1, 2, 3, 4, 5, 4, 5, 6, 1, 2}; */
  /* tsk_size_t sample_set_sizes[] = {3, 3, 3, 2}; */

  /* /\* int num_sample_sets = 1; *\/ */
  /* /\* tsk_id_t sample_sets[ts.num_samples]; *\/ */
  /* /\* for (int i = 0; i < ts.num_samples; i++) { *\/ */
  /* /\*   sample_sets[i] = i; *\/ */
  /* /\* } *\/ */
  /* /\* tsk_size_t sample_set_sizes[] = {ts.num_samples}; *\/ */

  /* /\* double *pi = calloc(num_sample_sets, sizeof(double)); *\/ */

  /* double pi[num_sample_sets]; */
  /* int ret = tsk_treeseq_diversity(&ts, num_sample_sets, sample_set_sizes, */
  /*                                 sample_sets, 0, NULL, TSK_STAT_SITE, pi);
   */

  tsk_ld_calc_t ld_calc;
  tsk_ld_calc_init(&ld_calc, &ts);

  double r2;
  int ret = tsk_ld_calc_get_r2(&ld_calc, 0, 2, &r2);
  if (ret != 0) {
    puts("error!");
    exit(1);
  }
}
