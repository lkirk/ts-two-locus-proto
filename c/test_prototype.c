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
test_bit_array_union(void)
{
    tsk_bit_array_t a[4] = { 4294967295, 4294967295, 67108863, 0 };
    tsk_bit_array_t b[4] = { 0, 0, 4294901760, 15 };
    tsk_bit_array_t truth_union[4] = { 0, 0, 67043328, 0 };

    tsk_bit_array_t result[4];
    union_bit_array(a, b, result, 4);

    puts("");
    for (int i = 0; i < 4; i++) {
        printf("%u\n", result[i]);
        CU_ASSERT_EQUAL(result[i], truth_union[i]);
    }
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_sample_counting", test_sample_counting },
        { "test_add_to_bit_array", test_add_to_bit_array },
        { "test_bit_array_union", test_bit_array_union },
        { NULL, NULL },
    };
    return test_main(tests, argc, argv);
}
