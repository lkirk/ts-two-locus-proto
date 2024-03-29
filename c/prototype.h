#pragma once
#include "tskit.h"
#include "summary_functions.h"

typedef uint32_t tsk_bit_array_t;
#define BIT_ARRAY_CHUNK ((tsk_bit_array_t) 5)
#define BIT_ARRAY_NUM_BITS ((tsk_bit_array_t) 1 << BIT_ARRAY_CHUNK)
#define BIT_ARRAY_NUM_CHUNKS(n)                                                         \
    (((n) >> BIT_ARRAY_CHUNK) + (((n) % BIT_ARRAY_NUM_BITS) ? 1 : 0))

#define GET_2D_ROW(array, row_len, row) (array + (((size_t) (row_len)) * (size_t) row))

typedef int(summary_func)(tsk_size_t, const double *, tsk_size_t, double *, void *);

int process_tree(summary_func *summary_function, const char *tree_filename,
    double **result, tsk_size_t *num_stat, bool print_weights);

void intersect_bit_array(const tsk_bit_array_t *a, const tsk_bit_array_t *b,
    tsk_bit_array_t *out, const tsk_size_t len);

void subtract_bit_arrays(
    tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len);

void add_bit_arrays(tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len);

void add_bit_to_bit_array(tsk_bit_array_t *a, const tsk_bit_array_t bit);

void count_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, tsk_size_t *c);

void print_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, int newline);

bool bit_in_array(const tsk_bit_array_t *a, const tsk_bit_array_t bit);

typedef int norm_func_t(tsk_size_t state_dim, const double *hap_weights, tsk_size_t n_a,
    tsk_size_t n_b, double *result, void *params);

int get_mutation_samples(const tsk_treeseq_t *self, const tsk_size_t tree_index,
    const tsk_size_t num_sample_chunks, const tsk_id_t *right_child,
    const tsk_id_t *left_sib, const tsk_id_t *parent, const tsk_size_t *site_offsets,
    tsk_bit_array_t *restrict *allele_samples, tsk_size_t *restrict *num_alleles);

void get_all_samples_bits(
    tsk_bit_array_t *all_samples, tsk_size_t n, tsk_size_t n_chunks);

int two_site_general_stat(const tsk_treeseq_t *self, tsk_size_t state_dim,
    const double *sample_weights, tsk_size_t result_dim, general_stat_func_t *f,
    norm_func_t *norm_f, void *f_params, tsk_size_t num_windows, const double *windows,
    tsk_flags_t options, double **result, tsk_size_t *num_stat, bool print_weights);

void sample_weights_to_bit_array(const double *weights, tsk_size_t num_sample_sets,
    tsk_size_t num_samples, tsk_size_t num_sample_chunks, double **total_weight,
    tsk_bit_array_t **sample_bits);

int get_all_mutation_samples(const tsk_treeseq_t *self,
    const tsk_size_t num_sample_chunks, const tsk_size_t *site_offsets,
    tsk_size_t *num_alleles, tsk_bit_array_t **allele_samples);

int compute_general_two_site_stat_result(tsk_size_t site_a, tsk_size_t site_a_offset,
    tsk_size_t site_b, tsk_size_t site_b_offset, tsk_size_t num_sample_chunks,
    const tsk_size_t *num_alleles, const tsk_bit_array_t *state, tsk_size_t state_dim,
    tsk_bit_array_t *sample_sets, tsk_size_t result_dim, general_stat_func_t *f,
    void *f_params, norm_func_t *norm_f, const double *total_weight, bool polarised,
    tsk_flags_t options, double *result, bool print_weights);
