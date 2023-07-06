#include <CUnit/Basic.h>

#include "CUnit/CUnit.h"
#include "testlib.h"
#include "prototype.h"

static void
test_sample_counting(void)
{
    tsk_size_t count;
    tsk_bit_array_t a[3] = { 4294967295, 4294967295, 67108863 };
    count_bit_array(a, 3, &count);
    CU_ASSERT_EQUAL(count, 90);
}

static void
test_add_to_bit_array(void)
{
    tsk_bit_array_t a[4] = { 0, 0, 0, 0 };
    for (tsk_bit_array_t s = 0; s < 90; s++) {
        add_bit_to_bit_array(a, s);
    }

    tsk_bit_array_t truth_a[4] = { 4294967295, 4294967295, 67108863, 0 };
    for (int i = 0; i < 3; i++) {
        CU_ASSERT_EQUAL(a[i], truth_a[i]);
    }

    tsk_bit_array_t b[4] = { 0, 0, 0, 0 };
    for (tsk_bit_array_t s = 80; s < 100; s++) {
        add_bit_to_bit_array(b, s);
    }

    tsk_bit_array_t truth_b[4] = { 0, 0, 4294901760, 15 };
    for (int i = 0; i < 4; i++) {
        CU_ASSERT_EQUAL(b[i], truth_b[i]);
    }
}

static void
test_bit_array_intersect(void)
{
    tsk_bit_array_t a[4] = { 4294967295, 4294967295, 67108863, 0 };
    tsk_bit_array_t b[4] = { 0, 0, 4294901760, 15 };
    tsk_bit_array_t truth_intersect[4] = { 0, 0, 67043328, 0 };

    tsk_bit_array_t result[4];
    intersect_bit_array(a, b, result, 4);

    puts("");
    for (int i = 0; i < 4; i++) {
        CU_ASSERT_EQUAL(result[i], truth_intersect[i]);
    }
}

static void
test_get_all_samples_bits(void)
{
    tsk_size_t num_samples;
    tsk_size_t num_sample_chunks;

    // number of samples does not fill the entire array
    num_samples = 124;
    num_sample_chunks = BIT_ARRAY_NUM_CHUNKS(num_samples);
    CU_ASSERT_EQUAL(num_sample_chunks, 4);
    tsk_bit_array_t result1[num_sample_chunks];
    get_all_samples_bits(result1, num_samples, num_sample_chunks);

    // NB: this test is valid for BIT_ARRAY_NUM_BITS = 32
    tsk_bit_array_t truth1[4] = { 4294967295, 4294967295, 4294967295, 268435455 };
    for (tsk_size_t i = 0; i < num_sample_chunks; i++) {
        CU_ASSERT_EQUAL(result1[i], truth1[i]);
    }

    // number of samples is an exact multiple of the BIT_ARRAY_NUM_BITS
    num_samples = 160;
    num_sample_chunks = BIT_ARRAY_NUM_CHUNKS(num_samples);
    CU_ASSERT_EQUAL(num_sample_chunks, 5);
    tsk_bit_array_t result2[num_sample_chunks];
    get_all_samples_bits(result2, num_samples, num_sample_chunks);

    // NB: this test is valid for BIT_ARRAY_NUM_BITS = 32
    tsk_bit_array_t truth2[5]
        = { 4294967295, 4294967295, 4294967295, 4294967295, 4294967295 };
    for (tsk_size_t i = 0; i < num_sample_chunks; i++) {
        CU_ASSERT_EQUAL(result2[i], truth2[i]);
    }
}

static void
test_sample_weights_to_bit_array(void)
{
    double weights[32] = { 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0,
        1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    tsk_size_t num_sample_sets = 4;
    tsk_size_t num_samples = 8;
    tsk_size_t num_sample_chunks = 1;
    tsk_bit_array_t *sample_bits
        = tsk_calloc(num_sample_chunks * num_sample_sets, sizeof(*sample_bits));
    double *total_weight = tsk_calloc(num_sample_sets, sizeof(*total_weight));
    tsk_bit_array_t bits_truth[4] = { 7, 56, 112, 6 };
    double total_truth[4] = { 3, 3, 3, 2 };
    sample_weights_to_bit_array(weights, num_sample_sets, num_samples, num_sample_chunks,
        &total_weight, &sample_bits);
    for (tsk_size_t i = 0; i < num_sample_sets * num_sample_chunks; i++) {
        CU_ASSERT_EQUAL(sample_bits[i], bits_truth[i]);
    }
    for (tsk_size_t i = 0; i < num_sample_sets; i++) {
        CU_ASSERT_EQUAL(total_weight[i], total_truth[i]);
    }
}

static void
test_bit_in_array(void)
{
    tsk_bit_array_t arr[2] = { 202, 2147483648 };
    CU_ASSERT_TRUE(bit_in_array(arr, 1));
    CU_ASSERT_TRUE(bit_in_array(arr, 3));
    CU_ASSERT_TRUE(bit_in_array(arr, 6));
    CU_ASSERT_TRUE(bit_in_array(arr, 7));
    CU_ASSERT_TRUE(bit_in_array(arr, 63));
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_sample_counting", test_sample_counting },
        { "test_add_to_bit_array", test_add_to_bit_array },
        { "test_bit_array_intersect", test_bit_array_intersect },
        { "test_get_all_samples_bits", test_get_all_samples_bits },
        { "test_sample_weights_to_bit_array", test_sample_weights_to_bit_array },
        { "test_bit_in_array", test_bit_in_array },
        { NULL, NULL },
    };
    return test_main(tests, argc, argv);
}
