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
    tsk_size_t *mut_offset, tsk_bit_array_t **mut_allele_samples,
    tsk_size_t **num_alleles)
{

    const tsk_size_t num_edges = self->tables->edges.num_rows;
    const tsk_flags_t *flags = self->tables->nodes.flags;
    const tsk_size_t num_sites = self->tree_sites_length[tree_index];
    const tsk_site_t *sites = self->tree_sites[tree_index];
    const tsk_id_t *mut_nodes = self->tables->mutations.node;
    const tsk_site_t *site;

    tsk_id_t u;
    int stack_top;
    tsk_bit_array_t *mut_samples_row;

    tsk_id_t top_mut_node = 0;
    tsk_size_t num_mutations = 0;

    for (tsk_size_t s = 0; s < num_sites; s++) {
        site = &sites[s];
        for (tsk_size_t m = 0; m < site->mutations_length; m++) {
            num_mutations++;
            if (site->mutations[m].node > top_mut_node) {
                top_mut_node = site->mutations[m].node;
            }
        }
    }

    tsk_bit_array_t *mut_samples
        = tsk_calloc(num_mutations * num_sample_chunks, sizeof(*mut_samples));

    // Stores the first and current node for a given mutation.
    tsk_id_t *first_curr_node = tsk_malloc(num_mutations * 2 * sizeof(*first_curr_node));
    // Sentinel value is -2 instead of -1
    tsk_memset(first_curr_node, 0xfe, num_mutations * 2 * sizeof(*first_curr_node));

    tsk_id_t *stack = tsk_malloc((1 + self->num_samples + num_edges) * sizeof(*stack));

    stack_top = 0;
    stack[stack_top] = parent[top_mut_node];

    tsk_id_t node;
    tsk_id_t *fc_row;
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

    tsk_bit_array_t *out_row;
    tsk_size_t state_offset = 0;
    // TODO: get rid of the global offsets
    for (tsk_size_t s = 0; s < num_sites; s++) {
        site = &self->tree_sites[tree_index][s];

        mut_samples_row = GET_2D_ROW(mut_samples, num_sample_chunks, state_offset);
        state_offset += self->site_mutations_length[site->id];

        /* printf("out_offset %lu\n", *out_offset); */

        out_row = GET_2D_ROW(*mut_allele_samples, num_sample_chunks, *out_offset);
        *out_offset += self->site_mutations_length[site->id] + 1;
        *mut_offset += self->site_mutations_length[site->id];
        get_allele_samples(site, num_sample_chunks, mut_samples_row, out_row,
            &(*num_alleles)[site->id]);
    }

    return 0;
}

/* int */
/* compute_general_two_stat_site_result() */
/* { */
/* } */

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

int
two_locus_stat(tsk_treeseq_t *self)
{
    const tsk_size_t num_nodes = self->tables->nodes.num_rows;
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

    tsk_id_t *restrict parent = tsk_malloc(num_nodes * sizeof(*parent));
    tsk_id_t *restrict right_child = tsk_malloc(num_nodes * sizeof(*right_child));
    tsk_id_t *restrict left_sib = tsk_malloc(num_nodes * sizeof(*left_sib));
    tsk_id_t *restrict right_sib = tsk_malloc(num_nodes * sizeof(*right_sib));

    tsk_memset(parent, 0xff, num_nodes * sizeof(*parent));
    tsk_memset(right_child, 0xff, num_nodes * sizeof(*right_child));
    tsk_memset(left_sib, 0xff, num_nodes * sizeof(*left_sib));
    tsk_memset(right_sib, 0xff, num_nodes * sizeof(*right_sib));

    tsk_size_t tree_index, num_sample_chunks;
    double t_left, t_right;
    tsk_id_t tj, tk, h, u, v, c;

    num_sample_chunks = (num_samples >> BIT_ARRAY_CHUNK)
                        + ((num_samples % BIT_ARRAY_NUM_BITS) ? 1 : 0);

    // TODO: this could be smaller if we concatenated mut_allele segments after
    //       computing the number of alleles
    // TODO: restrict?
    tsk_bit_array_t *mut_allele_samples
        = tsk_calloc(total_alleles * num_sample_chunks, sizeof(*mut_allele_samples));

    tsk_bit_array_t all_samples_bits[num_sample_chunks];
    get_all_samples_bits(all_samples_bits, num_samples, num_sample_chunks);

    tsk_size_t *site_offsets
        = tsk_malloc(self->tables->sites.num_rows * sizeof(*site_offsets));
    tsk_size_t *num_alleles
        = tsk_malloc(self->tables->sites.num_rows * sizeof(*site_offsets));
    tsk_bit_array_t *allele_samples_row;
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
    tsk_size_t out_offset = 0;
    tsk_size_t mut_offset = 0;

    tj = 0;
    tk = 0;
    t_left = 0;
    tree_index = 0;
    // TODO: check for TSK_NULLs
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

    // TODO: should free all unused arrays at this point.
    //       Maybe this should be inside of a function, then we call the below routine?
    tsk_size_t inner = 0;
    const tsk_bit_array_t *left_allele_samples;
    const tsk_bit_array_t *right_allele_samples;
    tsk_size_t left_offset = 0;
    tsk_size_t right_offset = 0;
    tsk_size_t w_A, w_B, w_AB, w_Ab, w_aB = 0;
    tsk_bit_array_t AB_samples[num_sample_chunks];

    for (tsk_size_t s_0 = 0; s_0 < self->tables->sites.num_rows; s_0++) {
        left_offset = site_offsets[s_0];
        for (tsk_size_t s_1 = inner; s_1 < self->tables->sites.num_rows; s_1++) {
            right_offset = site_offsets[s_1];
            for (tsk_size_t m_0 = 0; m_0 < num_alleles[s_0]; m_0++) {
                for (tsk_size_t m_1 = 0; m_1 < num_alleles[s_1]; m_1++) {
                    /* for (tsk_size_t m_0 = 0; m_0 < self->site_mutations_length[s_0] +
                     * 1; m_0++) { */
                    /*     for (tsk_size_t m_1 = 0; m_1 <
                     * self->site_mutations_length[s_1] + 1; */
                    /*          m_1++) { */
                    left_allele_samples = GET_2D_ROW(
                        mut_allele_samples, num_sample_chunks, left_offset + m_0);
                    right_allele_samples = GET_2D_ROW(
                        mut_allele_samples, num_sample_chunks, right_offset + m_1);

                    union_bit_array(left_allele_samples, right_allele_samples,
                        AB_samples, num_sample_chunks);
                    count_bit_array(AB_samples, num_sample_chunks, &w_AB);
                    count_bit_array(left_allele_samples, num_sample_chunks, &w_A);
                    count_bit_array(right_allele_samples, num_sample_chunks, &w_B);
                    w_Ab = w_A - w_AB;
                    w_aB = w_B - w_AB;
                    printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", s_0, s_1, m_0,
                        m_1, w_AB, w_Ab, w_aB, self->num_samples);
                }
            }
        }
        inner++;
    }

    return 0;
}
