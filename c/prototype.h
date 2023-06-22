#include "tskit.h"

typedef uint32_t tsk_bit_array_t;
#define BIT_ARRAY_CHUNK ((tsk_bit_array_t) 5)
#define BIT_ARRAY_NUM_BITS ((tsk_bit_array_t) 1 << BIT_ARRAY_CHUNK)

#define GET_2D_ROW(array, row_len, row) (array + (((size_t) (row_len)) * (size_t) row))

void union_bit_array(const tsk_bit_array_t *a, const tsk_bit_array_t *b,
    tsk_bit_array_t *out, const tsk_size_t len);

void subtract_bit_arrays(
    tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len);

void add_bit_arrays(tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len);

/* void add_bits_to_bit_array( */
/*     tsk_bit_array_t *a, const tsk_id_t *bits, const tsk_size_t len); */

/* void add_bit_to_bit_array(tsk_bit_array_t *a, const tsk_id_t bit); */
void add_bit_to_bit_array(tsk_bit_array_t *a, const tsk_bit_array_t bit);

void count_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, tsk_size_t *c);

void print_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, int newline);

/* static int get_allele_samples(const tsk_site_t *site, const tsk_size_t
 * num_sample_chunks, */
/*     const tsk_bit_array_t *state, tsk_bit_array_t *allele_samples); */

int get_mutation_samples(const tsk_treeseq_t *self, const tsk_size_t tree_index,
    const tsk_size_t num_sample_chunks, const tsk_id_t *right_child,
    const tsk_id_t *left_sib, const tsk_id_t *parent, tsk_size_t *out_offset,
    tsk_size_t *mut_offset, tsk_bit_array_t **mut_allele_samples);

int two_locus_stat(tsk_treeseq_t *self);
