#include <stdio.h>
#include <stdlib.h>

#include "tskit.h"
#include "tskit/core.h"
#include "tskit/tables.h"
#include "tskit/trees.h"

typedef uint32_t tsk_bit_array_t;
#define BIT_ARRAY_CHUNK ((tsk_bit_array_t) 5)
#define BIT_ARRAY_NUM_BITS ((tsk_bit_array_t) 1 << BIT_ARRAY_CHUNK)

#define GET_2D_ROW(array, row_len, row) (array + (((size_t) (row_len)) * (size_t) row))
/* static int */
/* cmp_uint32(const void *a, const void *b) */
/* { */
/*     const uint32_t *ia = (const uint32_t *) a; */
/*     const uint32_t *ib = (const uint32_t *) b; */
/*     return (*ia > *ib) - (*ia < *ib); */
/* } */

void
union_bit_array(const tsk_bit_array_t *a, const tsk_bit_array_t *b, tsk_bit_array_t *out,
    const tsk_size_t len)
{
    for (int i = 0; i < len; i++) {
        out[i] = a[i] & b[i];
    }
}

void
subtract_bit_arrays(tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len)
{
    for (int i = 0; i < len; i++) {
        a[i] &= ~(b[i]);
    }
}

void
add_bit_arrays(tsk_bit_array_t *a, const tsk_bit_array_t *b, const tsk_size_t len)
{
    for (int i = 0; i < len; i++) {
        a[i] |= b[i];
    }
}

void
add_bits_to_bit_array(tsk_bit_array_t *a, const tsk_id_t *bits, const tsk_size_t len)
{
    for (int i = 0; i < len; i++) {
        a[bits[i] >> BIT_ARRAY_CHUNK] |= 1 << (bits[i] - (bits[i] >> BIT_ARRAY_CHUNK));
    }
}

void
add_bit_to_bit_array(tsk_bit_array_t *a, const tsk_id_t bit)
{
    a[bit >> BIT_ARRAY_CHUNK] |= (tsk_bit_array_t) 1 << (bit - (bit >> BIT_ARRAY_CHUNK));
}

void
count_bit_array(tsk_bit_array_t *a, const tsk_size_t len, tsk_size_t *c)
{
    for (int i = 0; i < len; i++) {
        while (a[i]) {
            a[i] &= (a[i] - 1);
            c += 1;
        }
    }
}

void
print_bit_array(const tsk_bit_array_t *a, const tsk_size_t len, int newline)
{
    printf("[ ");
    for (int i = 0; i < len; i++) {
        printf("%i ", a[i]);
    }
    if (newline) {
        puts("]");
    } else {
        printf("]");
    }
}

static int
get_allele_samples(const tsk_site_t *site, const tsk_size_t num_sample_chunks,
    const tsk_bit_array_t *state, tsk_bit_array_t *allele_samples)
{
    int ret = 0;
    tsk_mutation_t mutation, parent_mut;
    tsk_size_t mutation_index, allele, num_alleles, alt_allele_length;
    /* The allele table */
    tsk_size_t max_alleles = site->mutations_length + 1;
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
    num_alleles = 1;

    print_bit_array(state, max_alleles, 1);
    for (mutation_index = 0; mutation_index < site->mutations_length; mutation_index++) {
        mutation = site->mutations[mutation_index];
        /* Compute the allele index for this derived state value. */
        allele = 0;
        while (allele < num_alleles) {
            if (mutation.derived_state_length == allele_lengths[allele]
                && tsk_memcmp(
                       mutation.derived_state, alleles[allele], allele_lengths[allele])
                       == 0) {
                break;
            }
            allele++;
        }
        if (allele == num_alleles) {
            tsk_bug_assert(allele < max_alleles);
            alleles[allele] = mutation.derived_state;
            allele_lengths[allele] = mutation.derived_state_length;
            num_alleles++;
        }

        /* Add the mutation's samples to this allele */
        allele_samples_row = GET_2D_ROW(allele_samples, num_sample_chunks, allele);
        state_row = GET_2D_ROW(state, num_sample_chunks, mutation_index);
        add_bit_arrays(allele_samples_row, state_row, num_sample_chunks);
        printf(" add %lu ", allele);
        print_bit_array(allele_samples, max_alleles, 1);

        /* Get the index for the alternate allele that we must substract from */
        alt_allele = site->ancestral_state;
        alt_allele_length = site->ancestral_state_length;
        if (mutation.parent != TSK_NULL) {
            parent_mut = site->mutations[mutation.parent - site->mutations[0].id];
            alt_allele = parent_mut.derived_state;
            alt_allele_length = parent_mut.derived_state_length;
        }
        allele = 0;
        while (allele < num_alleles) {
            if (alt_allele_length == allele_lengths[allele]
                && tsk_memcmp(alt_allele, alleles[allele], allele_lengths[allele])
                       == 0) {
                break;
            }
            allele++;
        }
        tsk_bug_assert(allele < num_alleles);

        alt_allele_samples_row = GET_2D_ROW(allele_samples, num_sample_chunks, allele);
        subtract_bit_arrays(
            alt_allele_samples_row, allele_samples_row, num_sample_chunks);
        printf(" sub %lu ", allele);
        print_bit_array(allele_samples, max_alleles, 1);
    }
out:
    tsk_safe_free(alleles);
    tsk_safe_free(allele_lengths);
    return ret;
}

int
two_locus_stat(tsk_treeseq_t *self)
{
    tsk_size_t num_nodes = self->tables->nodes.num_rows;
    const tsk_id_t num_edges = (tsk_id_t) self->tables->edges.num_rows;
    const tsk_id_t *restrict edges_in = self->tables->indexes.edge_insertion_order;
    const tsk_id_t *restrict edges_out = self->tables->indexes.edge_removal_order;
    const double *restrict edge_left = self->tables->edges.left;
    const double *restrict edge_right = self->tables->edges.right;
    const tsk_id_t *restrict edge_parent = self->tables->edges.parent;
    const tsk_id_t *restrict edge_child = self->tables->edges.child;
    const double sequence_length = self->tables->sequence_length;
    const tsk_flags_t *flags = self->tables->nodes.flags;

    tsk_id_t *restrict parent = tsk_malloc(num_nodes * sizeof(*parent));
    tsk_id_t *restrict left_child = tsk_malloc(num_nodes * sizeof(*left_child));
    tsk_id_t *restrict right_child = tsk_malloc(num_nodes * sizeof(*right_child));
    tsk_id_t *restrict left_sib = tsk_malloc(num_nodes * sizeof(*left_sib));
    tsk_id_t *restrict right_sib = tsk_malloc(num_nodes * sizeof(*right_sib));

    tsk_memset(parent, 0xff, num_nodes * sizeof(*parent));
    tsk_memset(left_child, 0xff, num_nodes * sizeof(*left_child));
    tsk_memset(right_child, 0xff, num_nodes * sizeof(*right_child));
    tsk_memset(left_sib, 0xff, num_nodes * sizeof(*left_sib));
    tsk_memset(right_sib, 0xff, num_nodes * sizeof(*right_sib));

    tsk_size_t tree_index;
    double t_left, t_right;
    tsk_id_t tj, tk, h, u, v, c;

    int stack_top;
    tsk_id_t *stack = tsk_malloc((1 + self->num_samples + num_edges) * sizeof(*stack));

    tsk_site_t site;
    tsk_size_t total_mutations = 0;
    tsk_size_t total_alleles = 0;
    for (int t = 0; t < self->num_trees; t++) {
        for (int s = 0; s < self->tree_sites_length[t]; s++) {
            total_alleles++; // one for the ancestral state
            site = self->tree_sites[t][s];
            for (int m = 0; m < self->site_mutations_length[site.id]; m++) {
                total_mutations++;
                total_alleles++;
            }
        }
    }
    printf("total_mutations %lu\n", total_mutations);
    tsk_size_t num_sample_chunks = (self->num_samples >> BIT_ARRAY_CHUNK)
                                           + (self->num_samples % BIT_ARRAY_NUM_BITS)
                                       ? 1
                                       : 0;

    // TODO: this could be smaller if we concatenated mut_allele segments after
    //       computing the number of alleles
    tsk_bit_array_t *mut_allele_samples
        = tsk_calloc(total_alleles * num_sample_chunks, sizeof(*mut_allele_samples));

    int tmp = num_sample_chunks - 1;
    tsk_bit_array_t all_samples_bits[num_sample_chunks];

    all_samples_bits[tmp]
        = ~(~((tsk_bit_array_t) 0) // TODO: are these casts necessary?
            << (self->num_samples % ((tsk_bit_array_t) 1 << BIT_ARRAY_CHUNK)));
    tmp--;
    while (tmp >= 0) {
        all_samples_bits[tmp] = ~((tsk_bit_array_t) 0);
        tmp--;
    }

    tsk_bit_array_t *allele_samples_row;
    int num_muts_cumsum = 0;
    for (int t = 0; t < self->num_trees; t++) {
        for (int s = 0; s < self->tree_sites_length[t]; s++) {
            allele_samples_row
                = GET_2D_ROW(mut_allele_samples, num_sample_chunks, num_muts_cumsum);
            add_bit_arrays(allele_samples_row, all_samples_bits, num_sample_chunks);
            num_muts_cumsum
                += self->site_mutations_length[self->tree_sites[t][s].id] + 1;
        }
    }
    /* for (int t = 0; t < self->num_trees; t++) { */
    /*     allele_samples_row = GET_2D_ROW( */
    /*         mut_allele_samples, num_sample_chunks, self->site_mutations_length[s]); */
    /*     for (int s = 0; s < self->tree_sites_length[t]; s++) { */
    /*     } */
    /* } */

    tsk_size_t out_offset = 0;

    tj = 0;
    tk = 0;
    t_left = 0;
    tree_index = 0;
    while (tj < num_edges || t_left < sequence_length) {
        while (tk < num_edges && edge_right[edges_out[tk]] == t_left) {
            h = edges_out[tk];
            tk++;
            u = edge_child[h];
            v = edge_parent[h];
            while (v != TSK_NULL) {
                if (left_sib[u] == TSK_NULL) {
                    left_child[v] = right_sib[u];
                } else {
                    right_sib[left_sib[u]] = right_sib[u];
                }
                if (right_sib[u] == TSK_NULL) {
                    right_child[v] = left_sib[u];
                } else {
                    left_sib[right_sib[u]] = left_sib[u];
                }
                left_sib[u] = -1;
                right_sib[u] = -1;

                v = parent[v];
            }
            parent[u] = -1;
        }
        while (tj < num_edges && edge_left[edges_in[tj]] == t_left) {
            h = edges_in[tj];
            tj++;
            u = edge_child[h];
            v = edge_parent[h];
            parent[u] = v;

            while (v != TSK_NULL) {
                c = right_child[v];
                if (c == TSK_NULL) {
                    left_child[v] = u;
                    left_sib[u] = -1;
                    right_sib[u] = -1;
                } else {
                    right_sib[c] = u;
                    left_sib[u] = c;
                    right_sib[u] = -1;
                }
                right_child[v] = u;

                v = parent[v];
            }
        }
        t_right = sequence_length;
        if (tj < num_edges) {
            t_right = TSK_MIN(t_right, edge_left[edges_in[tj]]);
        }
        if (tk < num_edges) {
            t_right = TSK_MIN(t_right, edge_right[edges_out[tk]]);
        }

        tsk_id_t top_mut_node = 0;
        tsk_size_t num_sites = self->tree_sites_length[tree_index];
        tsk_site_t *sites = self->tree_sites[tree_index];
        const tsk_site_t *site;
        tsk_size_t num_mutations = 0;
        tsk_id_t *mut_nodes = tsk_malloc(num_mutations * sizeof(*mut_nodes));

        for (int i = 0; i < num_sites; i++) {
            site = &sites[i];
            for (int j = 0; j < site->mutations_length; j++) {
                mut_nodes[num_mutations] = site->mutations[j].node;
                num_mutations++;
                if (site->mutations[j].node > top_mut_node) {
                    top_mut_node = site->mutations[j].node;
                }
            }
        }
        /* printf("top mut node is: %u \n", top_mut_node); */

        tsk_bit_array_t *mut_samples
            = tsk_calloc(num_mutations * num_sample_chunks, sizeof(*mut_samples));

        tsk_bit_array_t *mut_samples_row
            = GET_2D_ROW(mut_samples, num_mutations * num_sample_chunks, tree_index);
        tsk_bit_array_t *mut_samples_sub_row;

        tsk_id_t *first_curr_node
            = tsk_malloc(num_mutations * 2 * sizeof(*first_curr_node));
        tsk_memset(first_curr_node, 0xff, num_mutations * 2 * sizeof(*first_curr_node));

        // TODO: how to do this w/ memset
        // TODO: consider another sentinel value
        for (int i = 0; i < (num_mutations * 2); i++) {
            first_curr_node[i] *= 2;
        }

        stack_top = 0;
        stack[stack_top] = top_mut_node;

        tsk_id_t node;
        tsk_id_t *row;
        while (stack_top >= 0) {
            node = stack[stack_top];
            stack_top--;
            /* printf("N- %d\n", node); */
            for (int i = 0; i < num_mutations; i++) {
                // TODO: rename row to something better
                row = GET_2D_ROW(first_curr_node, 2, i);
                mut_samples_sub_row = GET_2D_ROW(mut_samples_row, num_sample_chunks, i);
                if (mut_nodes[i] == node) {
                    row[0] = node;
                    if (flags[node] & TSK_NODE_IS_SAMPLE) {
                        add_bit_to_bit_array(mut_samples_sub_row, node);
                    }
                }
                if (row[0] != -2) { // TODO: consider another sentinel value
                    if (row[1] == parent[node] || row[1] == left_sib[node]) {
                        row[1] = node;
                        if (flags[node] & TSK_NODE_IS_SAMPLE) {
                            add_bit_to_bit_array(mut_samples_sub_row, node);
                        }
                    } else if (row[0] == parent[node]) {
                        row[1] = node;
                        if (flags[node] & TSK_NODE_IS_SAMPLE) {
                            add_bit_to_bit_array(mut_samples_sub_row, node);
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
        for (int s = 0; s < num_sites; s++) {
            site = &self->tree_sites[tree_index][s];
            mut_samples_sub_row
                = GET_2D_ROW(mut_samples_row, num_sample_chunks, state_offset);
            state_offset += self->site_mutations_length[site->id];

            printf("out_offset %lu\n", out_offset);
            out_row = GET_2D_ROW(mut_allele_samples, num_sample_chunks, out_offset);
            out_offset += self->site_mutations_length[site->id] + 1;
            get_allele_samples(site, num_sample_chunks, mut_samples_sub_row, out_row);
        }

        printf("tree %lu \n", tree_index);
        printf("[ ");
        for (int i = 0; i < num_mutations * num_sample_chunks; i++) {
            printf("%d ", mut_samples_row[i]);
        }
        puts("]");

        tree_index++;
        t_left = t_right;
    }
    puts("final");
    print_bit_array(mut_allele_samples, total_alleles * num_sample_chunks, 1);
    return 0;
}

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
}
