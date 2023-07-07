#include <stdio.h>
#include <stdlib.h>

#include "prototype.h"

void
intersect_bit_array(const tsk_bit_array_t *a, const tsk_bit_array_t *b,
    tsk_bit_array_t *out, const tsk_size_t len)
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

void // Useful for debugging
print_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, int newline)
{
    tsk_bit_array_t val;
    printf("[ ");
    for (tsk_bit_array_t i = 0; i < len; i++) {
        for (tsk_bit_array_t j = 0; j < BIT_ARRAY_NUM_BITS; j++) {
            val = j + (BIT_ARRAY_NUM_BITS * i);
            if (bit_in_array(a, val)) {
                printf("%u ", val);
            }
        }
    }
    if (newline) {
        puts("]");
    } else {
        printf("]");
    }
}

bool
bit_in_array(const tsk_bit_array_t *a, const tsk_bit_array_t bit)
{
    tsk_bit_array_t i = bit >> BIT_ARRAY_CHUNK;
    return a[i] & ((tsk_bit_array_t) 1 << (bit - (BIT_ARRAY_NUM_BITS * i)));
}

static bool
node_in_array(tsk_bit_array_t *a, const tsk_id_t node)
{
    if (node == TSK_NULL) {
        return false;
    }
    return bit_in_array(a, (tsk_bit_array_t) node);
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
// TODO: the restricts here, not sure. Are the row declarations aliases? I think so.
{
    int ret = 0;

    const tsk_size_t num_edges = self->tables->edges.num_rows;
    const tsk_flags_t *restrict flags = self->tables->nodes.flags;
    const tsk_size_t num_sites = self->tree_sites_length[tree_index];
    const tsk_site_t *restrict sites = self->tree_sites[tree_index];
    const tsk_id_t *restrict mut_nodes = self->tables->mutations.node;
    const tsk_site_t *restrict site;

    tsk_id_t u, node, top_mut_node, *stack;
    tsk_bit_array_t *mut_samples, *mut_samples_row, *node_paths, *path, *out_row;
    tsk_size_t num_mutations, state_offset, num_nodes, num_node_chunks;
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
    num_nodes = tsk_treeseq_get_num_nodes(self);
    num_node_chunks = BIT_ARRAY_NUM_CHUNKS(num_nodes);

    node_paths = tsk_calloc(num_node_chunks * num_mutations, sizeof(*node_paths));

    // Size of stack is derived from tsk_tree_get_size_bound
    stack = tsk_malloc((1 + self->num_samples + num_edges) * sizeof(*stack));

    if (stack == NULL || node_paths == NULL || mut_samples == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }

    stack_top = 0;
    stack[stack_top] = parent[top_mut_node];

    // Traverse the current tree, tracking which samples are under each mutation
    while (stack_top >= 0) {
        node = stack[stack_top];
        stack_top--;
        for (tsk_size_t m = 0; m < num_mutations; m++) {
            path = GET_2D_ROW(node_paths, num_node_chunks, m);
            mut_samples_row = GET_2D_ROW(mut_samples, num_sample_chunks, m);
            if (mut_nodes[m + *mut_offset] == node) {
                add_bit_to_bit_array(path, (tsk_bit_array_t) node);
                if (flags[node] & TSK_NODE_IS_SAMPLE) {
                    add_bit_to_bit_array(mut_samples_row, (tsk_bit_array_t) node);
                }
            }
            if (node_in_array(path, parent[node])) {
                add_bit_to_bit_array(path, (tsk_bit_array_t) node);
                if (flags[node] & TSK_NODE_IS_SAMPLE) {
                    add_bit_to_bit_array(mut_samples_row, (tsk_bit_array_t) node);
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
    // TODO: don't we only want to do this from the current site offset onwards??
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
    tsk_safe_free(node_paths);
    tsk_safe_free(stack);

    return ret;
}

static int
norm_hap_weighted(tsk_size_t state_dim, const double *hap_weights,
    tsk_size_t TSK_UNUSED(n_a), tsk_size_t TSK_UNUSED(n_b), double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    const double *weight_row;
    double n;
    for (tsk_size_t k = 0; k < state_dim; k++) {
        weight_row = GET_2D_ROW(hap_weights, 3, k);
        n = (double) args.sample_set_sizes[k];
        // TODO: what to do when n = 0
        result[k] = weight_row[0] / n;
    }
    return 0;
}

static int
norm_af_weighted(tsk_size_t state_dim, const double *hap_weights,
    tsk_size_t TSK_UNUSED(n_a), tsk_size_t TSK_UNUSED(n_b), double *result, void *params)
{
    sample_count_stat_params_t args = *(sample_count_stat_params_t *) params;
    const double *weight_row;
    double p_A;
    double p_B;
    double n;
    for (tsk_size_t k = 0; k < state_dim; k++) {
        weight_row = GET_2D_ROW(hap_weights, 3, k);
        n = (double) args.sample_set_sizes[k];
        // TODO: what to do when n = 0
        p_A = (weight_row[0] + weight_row[1]) / n;
        p_B = (weight_row[0] + weight_row[2]) / n;
        result[k] = p_A * p_B;
    }
    return 0;
}

static int
norm_total_weighted(tsk_size_t state_dim, const double *TSK_UNUSED(hap_weights),
    tsk_size_t n_a, tsk_size_t n_b, double *result, void *TSK_UNUSED(params))
{
    for (tsk_size_t k = 0; k < state_dim; k++) {
        result[k] = 1 / (double) (n_a * n_b);
    }
    return 0;
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

// TODO: should site_1/2 be an id_t?
int
compute_general_two_site_stat_result(tsk_size_t site_1, tsk_size_t site_1_offset,
    tsk_size_t site_2, tsk_size_t site_2_offset, tsk_size_t num_sample_chunks,
    const tsk_size_t *num_alleles, const tsk_bit_array_t *state, tsk_size_t state_dim,
    tsk_bit_array_t *sample_sets, tsk_size_t result_dim, general_stat_func_t *f,
    void *f_params, norm_func_t *norm_func, const double *total_weight, bool polarised,
    tsk_flags_t TSK_UNUSED(options), double *result, bool print_weights)
{
    int ret = 0;
    const tsk_bit_array_t *A_samples, *B_samples;
    tsk_size_t w_A = 0, w_B = 0, w_AB = 0;
    tsk_bit_array_t *ss_row; // ss_ prefix refers to a sample set
    tsk_bit_array_t ss_A_samples[num_sample_chunks], ss_B_samples[num_sample_chunks],
        ss_AB_samples[num_sample_chunks], AB_samples[num_sample_chunks];
    uint8_t polarised_val = polarised ? 1 : 0;

    // TODO: naming of hap_weights?
    double hap_weights[3 * state_dim];
    double hap_norm[state_dim];
    double *hap_weight_row;
    tsk_size_t row_len = num_alleles[site_2] * state_dim;
    // TODO: is this stack allocation dangerous??
    double result_tmp[row_len * num_alleles[site_1]];
    double *result_tmp_row;

    for (tsk_size_t mut_1 = polarised_val; mut_1 < num_alleles[site_1]; mut_1++) {
        result_tmp_row = GET_2D_ROW(result_tmp, row_len, mut_1);
        for (tsk_size_t mut_2 = polarised_val; mut_2 < num_alleles[site_2]; mut_2++) {
            A_samples = GET_2D_ROW(state, num_sample_chunks, site_1_offset + mut_1);
            B_samples = GET_2D_ROW(state, num_sample_chunks, site_2_offset + mut_2);
            intersect_bit_array(A_samples, B_samples, AB_samples, num_sample_chunks);
            for (tsk_size_t k = 0; k < state_dim; k++) {
                ss_row = GET_2D_ROW(sample_sets, num_sample_chunks, k);
                hap_weight_row = GET_2D_ROW(hap_weights, 3, k);

                intersect_bit_array(A_samples, ss_row, ss_A_samples, num_sample_chunks);
                intersect_bit_array(B_samples, ss_row, ss_B_samples, num_sample_chunks);
                intersect_bit_array(
                    AB_samples, ss_row, ss_AB_samples, num_sample_chunks);

                count_bit_array(ss_AB_samples, num_sample_chunks, &w_AB);
                count_bit_array(ss_A_samples, num_sample_chunks, &w_A);
                count_bit_array(ss_B_samples, num_sample_chunks, &w_B);

                hap_weight_row[0] = (double) w_AB;
                hap_weight_row[1] = (double) (w_A - w_AB); // w_Ab
                hap_weight_row[2] = (double) (w_B - w_AB); // w_aB

                if (print_weights) {
                    printf("%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", site_1, site_2,
                        mut_1, mut_2, w_AB, (tsk_size_t) hap_weight_row[1],
                        (tsk_size_t) hap_weight_row[2], (tsk_size_t) total_weight[k]);
                }
            }
            ret = f(state_dim, hap_weights, result_dim, result_tmp_row, f_params);
            if (ret != 0) {
                goto out;
            }
            ret = norm_func(state_dim, hap_weights, num_alleles[site_1] - polarised_val,
                num_alleles[site_1] - polarised_val, hap_norm, f_params);
            if (ret != 0) {
                goto out;
            }
            for (tsk_size_t k = 0; k < state_dim; k++) {
                *result += result_tmp_row[k] * hap_norm[k];
            }
            result_tmp_row += state_dim; // Advance to the next column
        }
    }

out:
    return ret;
}

int
get_all_mutation_samples(const tsk_treeseq_t *self, const tsk_size_t num_sample_chunks,
    tsk_size_t *num_alleles, tsk_bit_array_t **mut_allele_samples)
{
    int ret = 0;
    const tsk_size_t num_nodes = self->tables->nodes.num_rows;
    const tsk_id_t num_edges = (tsk_id_t) self->tables->edges.num_rows;
    const tsk_id_t *restrict edges_in = self->tables->indexes.edge_insertion_order;
    const tsk_id_t *restrict edges_out = self->tables->indexes.edge_removal_order;
    const double *restrict edge_left = self->tables->edges.left;
    const double *restrict edge_right = self->tables->edges.right;
    const tsk_id_t *restrict edge_parent = self->tables->edges.parent;
    const tsk_id_t *restrict edge_child = self->tables->edges.child;
    const double sequence_length = self->tables->sequence_length;

    tsk_id_t *restrict parent = tsk_malloc(num_nodes * sizeof(*parent));
    tsk_id_t *restrict right_child = tsk_malloc(num_nodes * sizeof(*right_child));
    tsk_id_t *restrict left_sib = tsk_malloc(num_nodes * sizeof(*left_sib));
    tsk_id_t *restrict right_sib = tsk_malloc(num_nodes * sizeof(*right_sib));

    if (parent == NULL || right_child == NULL || left_sib == NULL || right_sib == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }

    tsk_memset(parent, 0xff, num_nodes * sizeof(*parent));
    tsk_memset(right_child, 0xff, num_nodes * sizeof(*right_child));
    tsk_memset(left_sib, 0xff, num_nodes * sizeof(*left_sib));
    tsk_memset(right_sib, 0xff, num_nodes * sizeof(*right_sib));

    tsk_size_t tree_index, out_offset, mut_offset;
    double t_left, t_right;
    tsk_id_t tj, tk, h, u, v, c;

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
            parent, &out_offset, &mut_offset, mut_allele_samples, &num_alleles);

        tree_index++;
        t_left = t_right;
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
    return ret;
}

int
two_site_general_stat(const tsk_treeseq_t *self, tsk_size_t state_dim,
    const double *sample_weights, tsk_size_t result_dim, general_stat_func_t *f,
    void *f_params, tsk_size_t TSK_UNUSED(num_windows),
    const double *TSK_UNUSED(windows), tsk_flags_t options, double **result,
    bool print_weights)
{
    int ret = 0;
    const tsk_size_t num_sites = self->tables->sites.num_rows;
    const tsk_size_t num_samples = self->num_samples;
    const tsk_size_t max_alleles = self->tables->mutations.num_rows + num_sites;
    tsk_size_t num_sample_chunks = BIT_ARRAY_NUM_CHUNKS(num_samples);
    tsk_size_t *restrict site_offsets = tsk_malloc(num_sites * sizeof(*site_offsets));
    tsk_size_t *restrict num_alleles = tsk_malloc(num_sites * sizeof(*num_alleles));
    tsk_bit_array_t *restrict sample_sets
        = tsk_calloc(num_sample_chunks * state_dim, sizeof(*sample_sets));
    tsk_bit_array_t *mut_allele_samples
        = tsk_calloc(max_alleles * num_sample_chunks, sizeof(*mut_allele_samples));
    tsk_bit_array_t all_samples_bits[num_sample_chunks];

    if (site_offsets == NULL || num_alleles == NULL || sample_sets == NULL
        || mut_allele_samples == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }

    get_all_samples_bits(all_samples_bits, num_samples, num_sample_chunks);

    // Initialize the mutation sample tracking with all samples in the ancestral
    // allele
    // TODO: rework this slightly to use the site offsets array
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

    get_all_mutation_samples(self, num_sample_chunks, num_alleles, &mut_allele_samples);

    // Number of pairs w/ replacement (sites)
    tsk_size_t num_stat = (num_sites * (1 + num_sites)) >> (tsk_size_t) 1;

    *result = tsk_malloc(num_stat * result_dim * sizeof(*result));
    double *total_weight = tsk_calloc(state_dim, sizeof(*total_weight));
    tsk_bit_array_t *sample_bits
        = tsk_calloc(num_sample_chunks * state_dim, sizeof(*sample_bits));

    if (total_weight == NULL || sample_bits == NULL || result == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }

    sample_weights_to_bit_array(sample_weights, state_dim, num_samples,
        num_sample_chunks, &total_weight, &sample_bits);
    bool polarised = false;

    if (options & TSK_STAT_POLARISED) {
        polarised = true;
    }

    norm_func_t *norm_func;
    if (options & TSK_HAP_WEIGHTED) {
        norm_func = norm_hap_weighted;
    } else if (options & TSK_AF_WEIGHTED) {
        norm_func = norm_af_weighted;
    } else if (options & TSK_TOTAL_WEIGHTED) {
        norm_func = norm_total_weighted;
    } else {
        // TODO: maybe we should have a "bad normalisation strategy" error??
        ret = TSK_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    tsk_size_t inner = 0, site_1_offset = 0, site_2_offset = 0, result_offset = 0;
    for (tsk_size_t site_1 = 0; site_1 < num_sites; site_1++) {
        site_1_offset = site_offsets[site_1];
        for (tsk_size_t site_2 = inner; site_2 < num_sites; site_2++) {
            site_2_offset = site_offsets[site_2];
            ret = compute_general_two_site_stat_result(site_1, site_1_offset, site_2,
                site_2_offset, num_sample_chunks, num_alleles, mut_allele_samples,
                state_dim, sample_bits, result_dim, f, f_params, norm_func, total_weight,
                polarised, options, &((*result)[result_offset]), print_weights);
            if (ret != 0) {
                goto out;
            }
            result_offset++;
        }
        inner++;
    }

    printf("result ");
    for (tsk_size_t i = 0; i < num_stat; i++) {
        printf("%.18f", (*result)[i]);
        if (i != num_stat - 1) {
            printf(",");
        } else {
            printf("\n");
        }
    }
out:
    if (mut_allele_samples != NULL) {
        free(mut_allele_samples);
    }
    if (site_offsets != NULL) {
        free(site_offsets);
    }
    if (num_alleles != NULL) {
        free(num_alleles);
    }
    if (sample_sets != NULL) {
        free(sample_sets);
    }
    tsk_safe_free(sample_bits);
    tsk_safe_free(total_weight);
    return ret;
}

static tsk_flags_t
pick_norm_strategy(summary_func *func)
{
    if (func == D) {
        return TSK_TOTAL_WEIGHTED;
    } else if (func == D2) {
        return TSK_TOTAL_WEIGHTED;
    } else if (func == r2) {
        return TSK_HAP_WEIGHTED;
    } else if (func == D_prime) {
        return TSK_HAP_WEIGHTED;
    } else if (func == r) {
        return TSK_TOTAL_WEIGHTED;
    } else if (func == Dz) {
        return TSK_TOTAL_WEIGHTED;
    } else if (func == pi2) {
        return TSK_TOTAL_WEIGHTED;
    }
    return 0;
}

static bool
pick_polarisation(summary_func *func)
{
    if (func == D) {
        return true;
    } else if (func == D2) {
        return false;
    } else if (func == r2) {
        return false;
    } else if (func == D_prime) {
        return true;
    } else if (func == r) {
        return true;
    } else if (func == Dz) {
        return false;
    } else if (func == pi2) {
        return false;
    }
    return false;
}

int
process_tree(
    summary_func *summary_function, const char *tree_filename, bool print_weights)
{
    tsk_flags_t options = 0;
    options |= pick_norm_strategy(summary_function);
    if (pick_polarisation(summary_function)) {
        options |= TSK_STAT_POLARISED;
    }

    tsk_treeseq_t ts;
    tsk_treeseq_load(&ts, tree_filename, 0);

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

    int ret
        = two_site_general_stat(&ts, num_sample_sets, sample_weights, num_sample_sets,
            summary_function, &args, 0, NULL, options, &result, print_weights);

    if (ret != 0) {
        puts(tsk_strerror(ret));
        exit(1);
    }

    tsk_treeseq_free(&ts);
    tsk_safe_free(result);
    tsk_safe_free(sample_weights);

    return ret;
}
