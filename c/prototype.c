#include <stdio.h>
#include <stdlib.h>

#include "prototype.h"

void
union_bit_array(const tsk_bit_array_t *a, const tsk_bit_array_t *b, tsk_bit_array_t *out,
    const tsk_size_t len)
{
    for (tsk_size_t i = 0; i < len; i++) {
        out[i] = a[i] & b[i];
    }
}

void
subtract_bit_arrays(tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len)
{
    for (tsk_size_t i = 0; i < len; i++) {
        a[i] &= ~(b[i]);
    }
}

void
add_bit_arrays(tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len)
{
    for (tsk_size_t i = 0; i < len; i++) {
        a[i] |= b[i];
    }
}

void
add_bit_to_bit_array(tsk_bit_array_t *a, const tsk_bit_array_t bit)
{
    tsk_bit_array_t i = bit >> BIT_ARRAY_CHUNK;
    a[i] |= (tsk_bit_array_t) 1 << (bit - (BIT_ARRAY_NUM_BITS * i));
}

void
count_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, tsk_size_t *c)
{
    tsk_bit_array_t tmp;
    *c = 0;
    for (tsk_size_t i = 0; i < len; i++) {
        tmp = a[i];
        while (tmp) {
            tmp &= (tmp - 1);
            *c += 1;
        }
    }
}

#define PRINT_ARRAY_GENERIC(a, len, format, newline)                                    \
    do {                                                                                \
        printf("[ ");                                                                   \
        for (int i = 0; i < len; i++) {                                                 \
            printf(format, a[i]);                                                       \
            printf(", ");                                                               \
        };                                                                              \
        if (newline) {                                                                  \
            puts("]");                                                                  \
        } else {                                                                        \
            printf("]");                                                                \
        }                                                                               \
    } while (0)

// TODO: figure inspect pass by ref/value here
static int
get_allele_samples(const tsk_site_t *site, const tsk_size_t num_sample_chunks,
    const tsk_bit_array_t *state, tsk_bit_array_t *allele_samples,
    tsk_size_t *num_alleles)
{
    int ret = 0;
    tsk_mutation_t mutation, parent_mut;
    tsk_size_t mutation_index, allele, alt_allele_length;
    /* The allele table */
    tsk_size_t max_alleles = site->mutations_length + 1;
    /* printf("max_alleles: %lu, %lu\n", max_alleles, max_alleles * sizeof(const char
     * *)); */
    const char **alleles = tsk_malloc(max_alleles * sizeof(*alleles));
    tsk_size_t *allele_lengths = tsk_calloc(max_alleles, sizeof(*allele_lengths));
    const char *alt_allele;
    const tsk_bit_array_t *state_row;
    tsk_bit_array_t *allele_samples_row;
    tsk_bit_array_t *alt_allele_samples_row;

    if (alleles == NULL || allele_lengths == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }

    tsk_bug_assert(state != NULL);
    alleles[0] = site->ancestral_state;
    allele_lengths[0] = site->ancestral_state_length;
    *num_alleles = 1;

    for (mutation_index = 0; mutation_index < site->mutations_length; mutation_index++) {
        mutation = site->mutations[mutation_index];
        /* Compute the allele index for this derived state value. */
        allele = 0;
        while (allele < *num_alleles) {
            if (mutation.derived_state_length == allele_lengths[allele]
                && tsk_memcmp(
                       mutation.derived_state, alleles[allele], allele_lengths[allele])
                       == 0) {
                break;
            }
            allele++;
        }
        if (allele == *num_alleles) {
            tsk_bug_assert(allele < max_alleles);
            alleles[allele] = mutation.derived_state;
            allele_lengths[allele] = mutation.derived_state_length;
            (*num_alleles)++;
        }

        /* Add the mutation's samples to this allele */
        allele_samples_row = GET_2D_ROW(allele_samples, num_sample_chunks, allele);
        state_row = GET_2D_ROW(state, num_sample_chunks, mutation_index);
        add_bit_arrays(allele_samples_row, state_row, num_sample_chunks);

        /* Get the index for the alternate allele that we must substract from */
        alt_allele = site->ancestral_state;
        alt_allele_length = site->ancestral_state_length;
        if (mutation.parent != TSK_NULL) {
            parent_mut = site->mutations[mutation.parent - site->mutations[0].id];
            alt_allele = parent_mut.derived_state;
            alt_allele_length = parent_mut.derived_state_length;
        }
        allele = 0;
        while (allele < *num_alleles) {
            if (alt_allele_length == allele_lengths[allele]
                && tsk_memcmp(alt_allele, alleles[allele], allele_lengths[allele])
                       == 0) {
                break;
            }
            allele++;
        }
        tsk_bug_assert(allele < *num_alleles);

        alt_allele_samples_row = GET_2D_ROW(allele_samples, num_sample_chunks, allele);
        subtract_bit_arrays(
            alt_allele_samples_row, allele_samples_row, num_sample_chunks);
    }
out:
    tsk_safe_free(alleles);
    tsk_safe_free(allele_lengths);
    return ret;
}

int
get_mutation_samples(const tsk_treeseq_t *self, const tsk_size_t tree_index,
    const tsk_size_t num_sample_chunks, const tsk_id_t *right_child,
    const tsk_id_t *left_sib, const tsk_id_t *parent, tsk_size_t *out_offset,
    tsk_size_t *mut_offset, tsk_bit_array_t *restrict *mut_allele_samples,
    tsk_size_t *restrict *num_alleles)
{
    int ret = 0;

    const tsk_size_t num_edges = self->tables->edges.num_rows;
    const tsk_flags_t *restrict flags = self->tables->nodes.flags;
    const tsk_size_t num_sites = self->tree_sites_length[tree_index];
    const tsk_site_t *restrict sites = self->tree_sites[tree_index];
    const tsk_id_t *restrict mut_nodes = self->tables->mutations.node;
    const tsk_site_t *restrict site;

    tsk_id_t u, node, top_mut_node, *fc_row, *first_curr_node, *stack;
    tsk_bit_array_t *mut_samples, *mut_samples_row, *out_row;
    tsk_size_t num_mutations, state_offset;
    int stack_top;

    // Get the number of mutations and the top mutation node
    num_mutations = 0;
    top_mut_node = 0;
    for (tsk_size_t s = 0; s < num_sites; s++) {
        site = &sites[s];
        for (tsk_size_t m = 0; m < site->mutations_length; m++) {
            num_mutations++;
            if (site->mutations[m].node > top_mut_node) {
                top_mut_node = site->mutations[m].node;
            }
        }
    }

    mut_samples = tsk_calloc(num_mutations * num_sample_chunks, sizeof(*mut_samples));
    first_curr_node = tsk_malloc(num_mutations * 2 * sizeof(*first_curr_node));
    // Size of stack is derived from tsk_tree_get_size_bound
    stack = tsk_malloc((1 + self->num_samples + num_edges) * sizeof(*stack));

    if (stack == NULL || first_curr_node == NULL || mut_samples == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }

    // Sentinel value is -2 instead of -1
    tsk_memset(first_curr_node, 0xfe, num_mutations * 2 * sizeof(*first_curr_node));

    stack_top = 0;
    stack[stack_top] = parent[top_mut_node];

    // Traverse the current tree, tracking which samples are under each mutation
    while (stack_top >= 0) {
        node = stack[stack_top];
        stack_top--;
        for (tsk_size_t m = 0; m < num_mutations; m++) {
            // TODO: rename row to something better
            fc_row = GET_2D_ROW(first_curr_node, 2, m);
            mut_samples_row = GET_2D_ROW(mut_samples, num_sample_chunks, m);
            if (mut_nodes[m + *mut_offset] == node) {
                fc_row[0] = node;
                if (flags[node] & TSK_NODE_IS_SAMPLE) {
                    add_bit_to_bit_array(mut_samples_row, (tsk_bit_array_t) node);
                }
            }
            if (fc_row[0] != -2) { // TODO: consider another sentinel value
                if (fc_row[1] == parent[node] || fc_row[1] == left_sib[node]) {
                    fc_row[1] = node;
                    if (flags[node] & TSK_NODE_IS_SAMPLE) {
                        add_bit_to_bit_array(mut_samples_row, (tsk_bit_array_t) node);
                    }
                } else if (fc_row[0] == parent[node]) {
                    fc_row[1] = node;
                    if (flags[node] & TSK_NODE_IS_SAMPLE) {
                        add_bit_to_bit_array(mut_samples_row, (tsk_bit_array_t) node);
                    }
                }
            }
        }
        u = right_child[node];
        while (u != TSK_NULL) {
            stack_top++;
            stack[stack_top] = u;
            u = left_sib[u];
        }
    }

    state_offset = 0;
    // TODO: get rid of the global offsets
    for (tsk_size_t s = 0; s < num_sites; s++) {
        site = &self->tree_sites[tree_index][s];

        mut_samples_row = GET_2D_ROW(mut_samples, num_sample_chunks, state_offset);
        state_offset += self->site_mutations_length[site->id];

        out_row = GET_2D_ROW(*mut_allele_samples, num_sample_chunks, *out_offset);
        *out_offset += self->site_mutations_length[site->id] + 1;
        *mut_offset += self->site_mutations_length[site->id];
        get_allele_samples(site, num_sample_chunks, mut_samples_row, out_row,
            &(*num_alleles)[site->id]);
    }

out:
    tsk_safe_free(mut_samples);
    tsk_safe_free(first_curr_node);
    tsk_safe_free(stack);

    return ret;
}

// TODO: static inline, keeping as is for testing
void
get_all_samples_bits(tsk_bit_array_t *all_samples, tsk_size_t n, tsk_size_t n_chunks)
{
    const tsk_bit_array_t all = ~((tsk_bit_array_t) 0);
    const tsk_bit_array_t remainder_samples = n % BIT_ARRAY_NUM_BITS;
    all_samples[n_chunks - 1] = remainder_samples ? ~(all << remainder_samples) : all;
    for (tsk_size_t i = 0; i < n_chunks - 1; i++) {
        all_samples[i] = all;
    }
}

void
sample_weights_to_bit_array(const double *weights, tsk_size_t num_sample_sets,
    tsk_size_t num_samples, tsk_size_t num_sample_chunks, double **total_weight,
    tsk_bit_array_t **sample_bits)
{
    const double *weight_row;
    tsk_bit_array_t *bits_row;
    for (tsk_bit_array_t i = 0; i < num_samples; i++) {
        weight_row = GET_2D_ROW(weights, num_sample_sets, i);
        for (tsk_size_t k = 0; k < num_sample_sets; k++) {
            (*total_weight)[k] += weight_row[k];
            if (weight_row[k] != 0) {
                bits_row = GET_2D_ROW(*sample_bits, num_sample_chunks, k);
                add_bit_to_bit_array(bits_row, i);
            }
        }
    }
}

// TODO: should site_1 be an id_t?
static int
compute_general_two_stat_site_result(tsk_size_t site_1, tsk_size_t site_1_offset,
    tsk_size_t site_2, tsk_size_t site_2_offset, tsk_size_t num_sample_chunks,
    const tsk_size_t *num_alleles, const tsk_bit_array_t *state, tsk_size_t state_dim,
    tsk_bit_array_t *sample_sets, tsk_size_t result_dim, general_stat_func_t *f,
    void *f_params, double *total_weight, bool polarised, double *result)
{
    int ret = 0;
    const tsk_bit_array_t *A_samples;
    const tsk_bit_array_t *B_samples;
    tsk_size_t w_A, w_B, w_AB, w_Ab, w_aB = 0;
    tsk_bit_array_t *sample_set_row;
    // ss_ prefix means that the samples are specific to a sample set
    tsk_bit_array_t ss_A_samples[num_sample_chunks], ss_B_samples[num_sample_chunks],
        ss_AB_samples[num_sample_chunks], AB_samples[num_sample_chunks];

    for (tsk_size_t mut_1 = polarised ? 1 : 0; mut_1 < num_alleles[site_1]; mut_1++) {
        for (tsk_size_t mut_2 = polarised ? 1 : 0; mut_2 < num_alleles[site_2];
             mut_2++) {
            A_samples = GET_2D_ROW(state, num_sample_chunks, site_1_offset + mut_1);
            B_samples = GET_2D_ROW(state, num_sample_chunks, site_2_offset + mut_2);

            union_bit_array(A_samples, B_samples, AB_samples, num_sample_chunks);
            for (tsk_size_t k = 0; k < state_dim; k++) {
                sample_set_row = GET_2D_ROW(sample_sets, num_sample_chunks, k);
                union_bit_array(
                    A_samples, sample_set_row, ss_A_samples, num_sample_chunks);
                union_bit_array(
                    B_samples, sample_set_row, ss_B_samples, num_sample_chunks);
                union_bit_array(
                    AB_samples, sample_set_row, ss_AB_samples, num_sample_chunks);
                count_bit_array(AB_samples, num_sample_chunks, &w_AB);
                count_bit_array(A_samples, num_sample_chunks, &w_A);
                count_bit_array(B_samples, num_sample_chunks, &w_B);

                w_Ab = w_A - w_AB;
                w_aB = w_B - w_AB;
                printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", site_1, site_2, mut_1,
                    mut_2, w_AB, w_Ab, w_aB, (tsk_size_t) total_weight[k]);
            }
        }
    }
    return ret;
}

int
two_site_general_stat(const tsk_treeseq_t *self, tsk_size_t state_dim,
    const double *sample_weights, tsk_size_t result_dim, general_stat_func_t *f,
    void *f_params, tsk_size_t num_windows, const double *windows, tsk_flags_t options,
    double *result)
{
    int ret = 0;
    const tsk_size_t num_nodes = self->tables->nodes.num_rows;
    const tsk_size_t num_sites = self->tables->sites.num_rows;
    const tsk_id_t num_edges = (tsk_id_t) self->tables->edges.num_rows;
    const tsk_id_t *restrict edges_in = self->tables->indexes.edge_insertion_order;
    const tsk_id_t *restrict edges_out = self->tables->indexes.edge_removal_order;
    const double *restrict edge_left = self->tables->edges.left;
    const double *restrict edge_right = self->tables->edges.right;
    const tsk_id_t *restrict edge_parent = self->tables->edges.parent;
    const tsk_id_t *restrict edge_child = self->tables->edges.child;
    const double sequence_length = self->tables->sequence_length;
    const tsk_size_t total_alleles
        = self->tables->mutations.num_rows + self->tables->sites.num_rows;
    const tsk_size_t num_samples = self->num_samples;

    tsk_size_t num_sample_chunks = (num_samples >> BIT_ARRAY_CHUNK)
                                   + ((num_samples % BIT_ARRAY_NUM_BITS) ? 1 : 0);

    tsk_id_t *restrict parent = tsk_malloc(num_nodes * sizeof(*parent));
    tsk_id_t *restrict right_child = tsk_malloc(num_nodes * sizeof(*right_child));
    tsk_id_t *restrict left_sib = tsk_malloc(num_nodes * sizeof(*left_sib));
    tsk_id_t *restrict right_sib = tsk_malloc(num_nodes * sizeof(*right_sib));

    tsk_size_t *restrict site_offsets = tsk_malloc(num_sites * sizeof(*site_offsets));
    tsk_size_t *restrict num_alleles = tsk_malloc(num_sites * sizeof(*num_alleles));
    tsk_bit_array_t *restrict mut_allele_samples
        = tsk_calloc(total_alleles * num_sample_chunks, sizeof(*mut_allele_samples));
    tsk_bit_array_t *restrict sample_sets
        = tsk_calloc(num_sample_chunks * state_dim, sizeof(*sample_sets));

    tsk_memset(parent, 0xff, num_nodes * sizeof(*parent));
    tsk_memset(right_child, 0xff, num_nodes * sizeof(*right_child));
    tsk_memset(left_sib, 0xff, num_nodes * sizeof(*left_sib));
    tsk_memset(right_sib, 0xff, num_nodes * sizeof(*right_sib));

    tsk_bit_array_t *allele_samples_row;
    tsk_size_t tree_index, out_offset, mut_offset;
    double t_left, t_right;
    tsk_id_t tj, tk, h, u, v, c;

    // TODO: check malloc return status

    tsk_bit_array_t all_samples_bits[num_sample_chunks];
    get_all_samples_bits(all_samples_bits, num_samples, num_sample_chunks);

    // Initialize the mutation sample tracking with all samples in the ancestral allele
    // TODO: rework this slightly to use the site offsets array
    tsk_size_t num_muts_cumsum = 0;
    for (tsk_size_t t = 0; t < self->num_trees; t++) {
        for (tsk_size_t s = 0; s < self->tree_sites_length[t]; s++) {
            allele_samples_row
                = GET_2D_ROW(mut_allele_samples, num_sample_chunks, num_muts_cumsum);
            add_bit_arrays(allele_samples_row, all_samples_bits, num_sample_chunks);
            site_offsets[self->tree_sites[t][s].id] = num_muts_cumsum;
            num_muts_cumsum
                += self->site_mutations_length[self->tree_sites[t][s].id] + 1;
        }
    }

    tj = 0;
    tk = 0;
    t_left = 0;
    tree_index = 0;
    out_offset = 0;
    mut_offset = 0;
    while (tj < num_edges || t_left < sequence_length) {
        while (tk < num_edges && edge_right[edges_out[tk]] == t_left) {
            h = edges_out[tk];
            tk++;
            u = edge_child[h];
            v = edge_parent[h];

            if (left_sib[u] != TSK_NULL) {
                right_sib[left_sib[u]] = right_sib[u];
            }

            if (right_sib[u] == TSK_NULL) {
                right_child[v] = left_sib[u];
            } else {
                left_sib[right_sib[u]] = left_sib[u];
            }

            left_sib[u] = TSK_NULL;
            right_sib[u] = TSK_NULL;
            parent[u] = TSK_NULL;
        }
        while (tj < num_edges && edge_left[edges_in[tj]] == t_left) {
            h = edges_in[tj];
            tj++;
            u = edge_child[h];
            v = edge_parent[h];
            parent[u] = v;

            c = right_child[v];
            if (c == TSK_NULL) {
                left_sib[u] = TSK_NULL;
                right_sib[u] = TSK_NULL;
            } else {
                right_sib[c] = u;
                left_sib[u] = c;
                right_sib[u] = TSK_NULL;
            }

            right_child[v] = u;
        }
        t_right = sequence_length;
        if (tj < num_edges) {
            t_right = TSK_MIN(t_right, edge_left[edges_in[tj]]);
        }
        if (tk < num_edges) {
            t_right = TSK_MIN(t_right, edge_right[edges_out[tk]]);
        }

        get_mutation_samples(self, tree_index, num_sample_chunks, right_child, left_sib,
            parent, &out_offset, &mut_offset, &mut_allele_samples, &num_alleles);

        tree_index++;
        t_left = t_right;
    }
    if (parent != NULL) {
        free(parent);
        parent = NULL;
    }
    if (right_child != NULL) {
        free(right_child);
        right_child = NULL;
    }
    if (left_sib != NULL) {
        free(left_sib);
        left_sib = NULL;
    }
    if (right_sib != NULL) {
        free(right_sib);
        right_sib = NULL;
    }

    // TODO: Maybe this should be inside of a function, then we call the below routine?
    tsk_size_t inner = 0;
    tsk_size_t site_1_offset = 0;
    tsk_size_t site_2_offset = 0;

    // Number of combinations w/ replacement (sites)
    tsk_size_t num_stat = (num_sites * (1 + num_sites)) >> (tsk_size_t) 1;
    result_dim = 1; // TODO: num sample sets
    // TODO check malloc return status
    result = tsk_malloc(num_stat * result_dim * sizeof(*result));
    double *total_weight = calloc(state_dim, sizeof(*total_weight));
    tsk_bit_array_t *sample_bits
        = tsk_calloc(num_sample_chunks * state_dim, sizeof(*sample_bits));
    sample_weights_to_bit_array(sample_weights, state_dim, num_samples,
        num_sample_chunks, &total_weight, &sample_bits);
    bool polarised = false;

    for (tsk_size_t site_1 = 0; site_1 < num_sites; site_1++) {
        site_1_offset = site_offsets[site_1];
        for (tsk_size_t site_2 = inner; site_2 < num_sites; site_2++) {
            site_2_offset = site_offsets[site_2];

            compute_general_two_stat_site_result(site_1, site_1_offset, site_2,
                site_2_offset, num_sample_chunks, num_alleles, mut_allele_samples,
                state_dim, sample_sets, result_dim, f, f_params, total_weight, polarised,
                result);
        }
        inner++;
    }
out:
    if (parent != NULL) {
        free(parent);
    }
    if (right_child != NULL) {
        free(right_child);
    }
    if (left_sib != NULL) {
        free(left_sib);
    }
    if (right_sib != NULL) {
        free(right_sib);
    }
    if (mut_allele_samples != NULL) {
        free(mut_allele_samples);
    }
    if (site_offsets != NULL) {
        free(site_offsets);
    }
    if (num_alleles != NULL) {
        free(num_alleles);
    }
    return ret;
}
