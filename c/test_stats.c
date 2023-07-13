#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>

#include <CUnit/Basic.h>

#include "test_stats_truth.h"
#include "testlib.h"
#include "prototype.h"

const char *test_cases[]
    = { "case1.tree", "case2.tree", /* "case3.tree",  */ "case4.tree", "case5.tree",
          "case6.tree", "case7.tree", "case11.tree", NULL };

// Caution, the below string/path code is not exactly robust/portable/safe! but it works
// for this limited purpose
static void
get_tree_dir(char *path)
{
    char *res;
    res = realpath("../trees", path);
    if (!res) {
        // hack, but gets the tests working
        res = realpath("../../trees", path);
    }
    if (!res) {
        char *err_str = strerror(errno);
        printf("ERROR: %s: %s\n", path, err_str);
        exit(1);
    }
}

static void
get_path_to_tree(const char *name, char *out)
{
    char tree_dir[PATH_MAX];

    get_tree_dir(tree_dir);
    strncpy(out, tree_dir, PATH_MAX);
    strncat(out, "/", PATH_MAX);
    strncat(out, name, PATH_MAX);

    if (access(out, F_OK) != 0) {
        printf("ERROR: %s: File does not exist\n", out);
        exit(1);
    }
}

static void
assert_results_equal_truth(double *truth, double *result, tsk_size_t n)
{
    double tolerance = 1e-17;
    for (tsk_size_t i = 0; i < n; i++) {
        /* printf("%f %f %f\n", truth[i] - result[i], truth[i], result[i]); */
        CU_ASSERT_DOUBLE_EQUAL(truth[i], result[i], tolerance);
    }
}

static void
test_D(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(D, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(D_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

static void
test_D2(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(D2, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(D2_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

static void
test_r2(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(r2, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(r2_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

static void
test_D_prime(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(D_prime, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(D_prime_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

static void
test_r(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(r, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(r_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

static void
test_Dz(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(Dz, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(Dz_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

static void
test_pi2(void)
{
    char tree_path[PATH_MAX];
    double *result;
    tsk_size_t num_stat;
    for (int i = 0; test_cases[i] != NULL; i++) {
        get_path_to_tree(test_cases[i], tree_path);
        process_tree(pi2, tree_path, &result, &num_stat, false);
        assert_results_equal_truth(pi2_truth[i], result, num_stat);
        tsk_safe_free(result);
    }
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_D", test_D },
        { "test_D2", test_D2 },
        { "test_r2", test_r2 },
        { "test_D_prime", test_D_prime },
        { "test_r", test_r },
        { "test_Dz", test_Dz },
        { "test_pi2", test_pi2 },

        { NULL, NULL },
    };
    return test_main(tests, argc, argv);
}
