#include <CUnit/Basic.h>

#include "CUnit/CUnit.h"
#include "testlib.h"
#include "prototype.h"

static void
test_sample_counting(void)
{
    tsk_size_t count;
    tsk_bit_array_t b[3] = { 4294967295, 4294967295, 67108863 };
    count_bit_array(b, 3, &count);
    CU_ASSERT_EQUAL(count, 90);
}

static void
test_bit_arrays(void)
{
    tsk_bit_array_t b[3] = { 0, 0, 0 };
    for (tsk_bit_array_t s = 0; s < 90; s++) {
        add_bit_to_bit_array(b, s);
    }

    tsk_bit_array_t truth[3] = { 4294967295, 4294967295, 67108863 };
    for (int i = 0; i < 3; i++) {
        CU_ASSERT_EQUAL(b[i], truth[i]);
    }
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_sample_counting", test_sample_counting },
        { "test_bit_arrays", test_bit_arrays },
        { NULL, NULL },
    };
    return test_main(tests, argc, argv);
}
