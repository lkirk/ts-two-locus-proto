double D_truth[11][1024] = {
    { 0.2222222222222222, -0.1111111111111111, 0.2222222222222222 },
    { 0.046875, 0.0625, 0.0625, 0.25, 0.25, 0.25 },
    /* { 0.1875, -0.078125, 0.1875, 0.03125, 0.0, -0.010416666666666671, 0.0, 0.1875, */
    /*     -0.03125, 0.05859375, -0.078125, -0.0234375, 0.0, 0.0078125, 0.0, -0.078125,
     */
    /*     0.0234375, 0.1875, 0.03125, 0.0, -0.010416666666666671, 0.0, 0.1875, -0.03125,
     */
    /*     0.109375, 0.0, 0.005208333333333334, 0.0, 0.03125, 0.015625, 0.0, 0.0, 0.0,
       0.0, */
    /*     0.0, 0.012152777777777773, 0.0, -0.010416666666666671, 0.008680555555555552,
       0.0, */
    /*     0.0, 0.0, 0.1875, -0.03125, 0.026041666666666668 }, */
    { 0.0625, 0.08333333333333333, 0.08333333333333333, 0.041666666666666685,
        0.2222222222222222, 0.1111111111111111, 0.05555555555555558, 0.2222222222222222,
        0.11111111111111116, 0.13888888888888884 },
    { 0.1875, 0.1875, 0.1875, 0.1875, 0.1875, 0.03125, -0.09375, -0.1875, 0.1875,
        0.03125, 0.1875, -0.1875, 0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125,
        -0.1875, -0.1875, -0.1875, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875,
        0.1875, 0.1875, 0.1875, 0.1875, 0.03125, -0.09375, -0.1875, 0.1875, 0.03125,
        0.1875, -0.1875, 0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125, -0.1875,
        -0.1875, -0.1875, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.1875,
        0.1875, 0.1875, 0.03125, -0.09375, -0.1875, 0.1875, 0.03125, 0.1875, -0.1875,
        0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125, -0.1875, -0.1875, -0.1875,
        -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.1875, 0.1875, 0.03125,
        -0.09375, -0.1875, 0.1875, 0.03125, 0.1875, -0.1875, 0.1875, 0.1875, -0.1875,
        -0.1875, 0.1875, 0.03125, -0.1875, -0.1875, -0.1875, -0.1875, 0.1875, -0.1875,
        0.125, 0.03125, 0.1875, 0.1875, 0.03125, -0.09375, -0.1875, 0.1875, 0.03125,
        0.1875, -0.1875, 0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125, -0.1875,
        -0.1875, -0.1875, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.109375,
        -0.015625, -0.03125, 0.03125, -0.015625, 0.03125, -0.03125, 0.03125, 0.03125,
        -0.03125, -0.03125, 0.03125, -0.015625, -0.03125, -0.03125, -0.03125, -0.03125,
        0.03125, -0.03125, -0.0625, -0.015625, 0.03125, 0.109375, 0.09375, -0.09375,
        -0.015625, -0.09375, 0.09375, -0.09375, -0.09375, 0.09375, 0.09375, -0.09375,
        -0.015625, 0.09375, 0.09375, 0.09375, 0.09375, -0.09375, 0.09375, -0.0625,
        -0.015625, -0.09375, 0.1875, -0.1875, -0.03125, -0.1875, 0.1875, -0.1875,
        -0.1875, 0.1875, 0.1875, -0.1875, -0.03125, 0.1875, 0.1875, 0.1875, 0.1875,
        -0.1875, 0.1875, -0.125, -0.03125, -0.1875, 0.1875, 0.03125, 0.1875, -0.1875,
        0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125, -0.1875, -0.1875, -0.1875,
        -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.109375, 0.03125, -0.03125,
        0.03125, 0.03125, -0.03125, -0.03125, 0.03125, -0.015625, -0.03125, -0.03125,
        -0.03125, -0.03125, 0.03125, -0.03125, -0.0625, -0.015625, 0.03125, 0.1875,
        -0.1875, 0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125, -0.1875, -0.1875,
        -0.1875, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.1875, -0.1875,
        -0.1875, 0.1875, 0.1875, -0.1875, -0.03125, 0.1875, 0.1875, 0.1875, 0.1875,
        -0.1875, 0.1875, -0.125, -0.03125, -0.1875, 0.1875, 0.1875, -0.1875, -0.1875,
        0.1875, 0.03125, -0.1875, -0.1875, -0.1875, -0.1875, 0.1875, -0.1875, 0.125,
        0.03125, 0.1875, 0.1875, -0.1875, -0.1875, 0.1875, 0.03125, -0.1875, -0.1875,
        -0.1875, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.1875, 0.1875,
        -0.1875, -0.03125, 0.1875, 0.1875, 0.1875, 0.1875, -0.1875, 0.1875, -0.125,
        -0.03125, -0.1875, 0.1875, -0.1875, -0.03125, 0.1875, 0.1875, 0.1875, 0.1875,
        -0.1875, 0.1875, -0.125, -0.03125, -0.1875, 0.1875, 0.03125, -0.1875, -0.1875,
        -0.1875, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.109375, -0.03125,
        -0.03125, -0.03125, -0.03125, 0.03125, -0.03125, 0.0625, -0.015625, 0.03125,
        0.1875, 0.1875, 0.1875, 0.1875, -0.1875, 0.1875, -0.125, -0.03125, -0.1875,
        0.1875, 0.1875, 0.1875, -0.1875, 0.1875, -0.125, -0.03125, -0.1875, 0.1875,
        0.1875, -0.1875, 0.1875, -0.125, -0.03125, -0.1875, 0.1875, -0.1875, 0.1875,
        -0.125, -0.03125, -0.1875, 0.1875, -0.1875, 0.125, 0.03125, 0.1875, 0.1875,
        -0.125, -0.03125, -0.1875, 0.25, 0.0625, 0.125, 0.109375, 0.03125, 0.1875 },
    { 0.25, 0.25, 0.25 },
    { 0.043209876543209874, -0.018518518518518517, 0.05555555555555555 },
    { 0.05555555555555555, 0.0, 0.05555555555555555 },
    { 0.25, 0.0, 0.25 },
    { 0.25, -0.25, 0.25 },
    { 0.05555555555555555, 0.037037037037037035, 0.043209876543209874 },
};

double D2_truth[11][1024] = {
    { 0.04938271604938271, 0.012345679012345678, 0.04938271604938271 },
    { 0.022569444444444444, 0.03125, 0.03125, 0.0625, 0.0625, 0.0625 },
    /* { 0.03515625, 0.012369791666666666, 0.03515625, 0.0009765625, 0.0, 0.01513671875,
     */
    /*     0.0, 0.03515625, 0.01220703125, 0.01833767361111111, 0.012369791666666664, */
    /*     0.0011393229166666667, 0.0, 0.004964192708333333, 0.0, 0.012369791666666664,
     */
    /*     0.006429036458333332, 0.03515625, 0.0009765625, 0.0, 0.01513671875, 0.0, */
    /*     0.03515625, 0.01220703125, 0.011962890625, 0.0, 0.0008544921875, 0.0, */
    /*     0.0009765625, 0.0030517578125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01031494140625, 0.0,
     */
    /*     0.01513671875, 0.00579833984375, 0.0, 0.0, 0.0, 0.03515625, 0.01220703125, */
    /*     0.01177978515625 }, */
    { 0.022976680384087792, 0.026748971193415634, 0.014403292181069956,
        0.00360082304526749, 0.04938271604938271, 0.012345679012345678,
        0.003086419753086418, 0.04938271604938271, 0.01234567901234568,
        0.019290123456790122 },
    { 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625,
        0.0087890625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0009765625, 0.0087890625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0009765625, 0.0087890625, 0.03515625, 0.03515625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.0009765625, 0.0087890625, 0.03515625, 0.03515625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.03515625,
        0.0009765625, 0.0087890625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.011962890625, 0.000244140625,
        0.0009765625, 0.0009765625, 0.000244140625, 0.0009765625, 0.0009765625,
        0.0009765625, 0.0009765625, 0.0009765625, 0.0009765625, 0.0009765625,
        0.000244140625, 0.0009765625, 0.0009765625, 0.0009765625, 0.0009765625,
        0.0009765625, 0.0009765625, 0.00390625, 0.000244140625, 0.0009765625,
        0.011962890625, 0.0087890625, 0.0087890625, 0.000244140625, 0.0087890625,
        0.0087890625, 0.0087890625, 0.0087890625, 0.0087890625, 0.0087890625,
        0.0087890625, 0.000244140625, 0.0087890625, 0.0087890625, 0.0087890625,
        0.0087890625, 0.0087890625, 0.0087890625, 0.00390625, 0.000244140625,
        0.0087890625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625,
        0.0009765625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625,
        0.0009765625, 0.03515625, 0.011962890625, 0.0009765625, 0.0009765625,
        0.0009765625, 0.0009765625, 0.0009765625, 0.0009765625, 0.0009765625,
        0.000244140625, 0.0009765625, 0.0009765625, 0.0009765625, 0.0009765625,
        0.0009765625, 0.0009765625, 0.00390625, 0.000244140625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625,
        0.0009765625, 0.03515625, 0.03515625, 0.0009765625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625, 0.0009765625,
        0.03515625, 0.011962890625, 0.0009765625, 0.0009765625, 0.0009765625,
        0.0009765625, 0.0009765625, 0.0009765625, 0.00390625, 0.000244140625,
        0.0009765625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.015625, 0.0009765625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.015625, 0.0009765625,
        0.03515625, 0.03515625, 0.03515625, 0.015625, 0.0009765625, 0.03515625,
        0.03515625, 0.015625, 0.0009765625, 0.03515625, 0.0625, 0.00390625, 0.015625,
        0.011962890625, 0.0009765625, 0.03515625 },
    { 0.0625, 0.0625, 0.0625 },
    { 0.023844603634269844, 0.02384460363426984, 0.02384460363426984 },
    { 0.024691358024691357, 0.0, 0.024691358024691357 },
    { 0.0625, 0.0, 0.0625 },
    { 0.0625, 0.0625, 0.0625 },
    { 0.02384460363426984, 0.009212687767786075, 0.023844603634269844 },
};

double r2_truth[11][1024] = { { 0.9999999999999996, 0.2499999999999999,
                                  0.9999999999999996 },
    { 1.0, 0.6666666666666666, 0.6666666666666666, 1.0, 1.0, 1.0 },
    /* { 1.0, 0.42857142857142855, 1.0, 0.047619047619047616, 0.0, 0.6031746031746031,
       0.0, */
    /*     1.0, 0.35873015873015873, 1.0, 0.42857142857142855, 0.06802721088435375, 0.0,
     */
    /*     0.18984126984126984, 0.0, 0.42857142857142855, 0.22984126984126987, 1.0, */
    /*     0.047619047619047616, 0.0, 0.6031746031746031, 0.0, 1.0, 0.35873015873015873,
     */
    /*     1.0, 0.0, 0.06802721088435375, 0.0, 0.047619047619047616, 0.15374149659863945,
     */
    /*     0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.6031746031746031, 0.3415873015873016,
       0.0, */
    /*     0.0, 0.0, 1.0, 0.35873015873015873, 1.0 }, */
    { 1.0, 0.5999999999999999, 0.3499999999999999, 0.13999999999999993,
        0.9999999999999996, 0.2499999999999999, 0.09999999999999984, 0.9999999999999996,
        0.4, 0.9999999999999998 },
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616, 0.42857142857142855, 1.0, 1.0,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0,
        1.0, 1.0, 1.0, 0.047619047619047616, 0.42857142857142855, 1.0, 1.0,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0,
        1.0, 1.0, 0.047619047619047616, 0.42857142857142855, 1.0, 1.0,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0,
        1.0, 0.047619047619047616, 0.42857142857142855, 1.0, 1.0, 0.047619047619047616,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0, 0.047619047619047616,
        0.42857142857142855, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 0.02040816326530612, 0.047619047619047616,
        0.047619047619047616, 0.02040816326530612, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.02040816326530612, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.14285714285714285, 0.02040816326530612,
        0.047619047619047616, 1.0, 0.4285714285714285, 0.4285714285714285,
        0.02040816326530612, 0.4285714285714285, 0.4285714285714285, 0.4285714285714285,
        0.4285714285714285, 0.4285714285714285, 0.4285714285714285, 0.4285714285714285,
        0.02040816326530612, 0.4285714285714285, 0.4285714285714285, 0.4285714285714285,
        0.4285714285714285, 0.4285714285714285, 0.4285714285714285, 0.14285714285714285,
        0.02040816326530612, 0.4285714285714285, 1.0, 1.0, 0.047619047619047616, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.3333333333333333, 0.047619047619047616, 1.0, 1.0, 0.047619047619047616, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.3333333333333333, 0.047619047619047616, 1.0, 1.0, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.02040816326530612, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.14285714285714285, 0.02040816326530612,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.3333333333333333, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 1.0, 0.047619047619047616, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.047619047619047616, 0.047619047619047616,
        0.047619047619047616, 0.14285714285714285, 0.02040816326530612,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 1.0, 1.0, 0.3333333333333333,
        0.047619047619047616, 1.0, 1.0, 1.0, 0.3333333333333333, 0.047619047619047616,
        1.0, 1.0, 0.3333333333333333, 0.047619047619047616, 1.0, 1.0,
        0.14285714285714285, 0.3333333333333333, 1.0, 0.047619047619047616, 1.0 },
    { 1.0, 1.0, 1.0 }, { 0.9999999999999998, 0.9999999999999998, 0.9999999999999998 },
    { 0.9999999999999996, 0.0, 0.9999999999999996 }, { 1.0, 0.0, 1.0 },
    { 1.0, 1.0, 1.0 }, { 0.9999999999999998, 0.23079365079365075, 0.9999999999999998 } };

double D_prime_truth[11][1024] = {
    { 0.33333333333333326, 0.0, 0.33333333333333326 },
    { 0.75, 0.5, 0.5, 0.5, 0.5, 0.5 },
    /* { 0.75, -0.041666666666666664, 0.75, 0.125, nan, nan, 0.0, 0.75, 0.375, 0.375, */
    /*     -0.041666666666666664, 0.0, nan, nan, 0.0, -0.041666666666666664, */
    /*     0.20833333333333331, 0.75, 0.125, nan, nan, 0.0, 0.75, 0.375, 0.125, nan, nan,
     */
    /*     0.0, 0.125, 0.125, nan, nan, 0.0, nan, nan, nan, 0.0, nan, nan, 0.0, 0.0, 0.0,
     */
    /*     0.75, 0.375, 0.625 }, */
    { 0.4999999999999999, 0.33333333333333326, 0.4999999999999999, 0.5,
        0.33333333333333326, 0.33333333333333326, 0.3333333333333333, 0.6666666666666665,
        0.6666666666666666, 0.8333333333333331 },
    { 0.75, 0.75, 0.75, 0.75, 0.75, 0.125, 0.0, 0.0, 0.75, 0.125, 0.75, 0.0, 0.75, 0.75,
        0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.75,
        0.75, 0.75, 0.75, 0.125, 0.0, 0.0, 0.75, 0.125, 0.75, 0.0, 0.75, 0.75, 0.0, 0.0,
        0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.75, 0.75, 0.75,
        0.125, 0.0, 0.0, 0.75, 0.125, 0.75, 0.0, 0.75, 0.75, 0.0, 0.0, 0.75, 0.125, 0.0,
        0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.75, 0.75, 0.125, 0.0, 0.0, 0.75,
        0.125, 0.75, 0.0, 0.75, 0.75, 0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75,
        0.0, 0.5, 0.125, 0.75, 0.75, 0.125, 0.0, 0.0, 0.75, 0.125, 0.75, 0.0, 0.75, 0.75,
        0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.125,
        0.0, 0.0, 0.125, 0.0, 0.125, 0.0, 0.125, 0.125, 0.0, 0.0, 0.125, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.125, 0.125, 0.125, 0.0, 0.0, 0.0, 0.125, 0.0,
        0.0, 0.125, 0.125, 0.0, 0.0, 0.125, 0.125, 0.125, 0.125, 0.0, 0.125, 0.0, 0.0,
        0.0, 0.25, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.25, 0.25, 0.25,
        0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.75, 0.125, 0.75, 0.0, 0.75, 0.75, 0.0, 0.0,
        0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.125, 0.125, 0.0,
        0.125, 0.125, 0.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.0,
        0.125, 0.75, 0.0, 0.75, 0.75, 0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75,
        0.0, 0.5, 0.125, 0.75, 0.25, 0.0, 0.0, 0.25, 0.25, 0.0, 0.0, 0.25, 0.25, 0.25,
        0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.75, 0.75, 0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0,
        0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.75, 0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0,
        0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.25, 0.25, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25,
        0.0, 0.25, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25, 0.0, 0.25, 0.0,
        0.0, 0.0, 0.75, 0.125, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75, 0.125,
        0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.125, 0.0, 0.125, 0.25, 0.25, 0.25, 0.25, 0.0,
        0.25, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.25, 0.25, 0.0,
        0.25, 0.0, 0.0, 0.0, 0.25, 0.0, 0.25, 0.0, 0.0, 0.0, 0.75, 0.0, 0.5, 0.125, 0.75,
        0.25, 0.0, 0.0, 0.0, 0.5, 0.125, 0.5, 0.125, 0.125, 0.75 },
    { 0.5, 0.5, 0.5 },
    { 0.7777777777777777, 0.4444444444444444, 0.6666666666666666 },
    { 0.6666666666666665, 0.0, 0.6666666666666665 },
    { 0.5, 0.0, 0.5 },
    { 0.5, 0.0, 0.5 },
    { 0.6666666666666666, 0.3333333333333333, 0.7777777777777777 },
};

double r_truth[11][1024] = {
    { 0.9999999999999999, -0.49999999999999994, 0.9999999999999999 },
    { 0.21132486540518708, 0.21132486540518708, 0.21132486540518708, 1.0, 1.0, 1.0 },
    /* { 1.0, -0.4939935020206553, 1.0, 0.2182178902359924, 0.0, -0.08488133583335669,
       0.0, */
    /*     1.0, -0.1494829254768914, 0.3908910548820038, -0.4939935020206553, */
    /*     -0.18053751654656763, 0.0, 0.05849429484450966, 0.0, -0.4939935020206553, */
    /*     0.1236712356742527, 1.0, 0.2182178902359924, 0.0, -0.08488133583335669,
       0.0, 1.0, */
    /*     -0.1494829254768914, 1.0, 0.0, 0.024850710549522523, 0.0, 0.2182178902359924,
     */
    /*     0.09785954587161398, 0.0, 0.0, 0.0, 0.0, 0.0, 0.056587557222237794, 0.0, */
    /*     -0.08488133583335669, 0.049230711257251356, 0.0, 0.0, 0.0, 1.0, */
    /*     -0.1494829254768914, 0.16227353026548486 }, */
    { 0.341886116991581, 0.341886116991581, 0.4081138830084189, 0.25811388300841903,
        0.9999999999999999, 0.49999999999999994, 0.316227766016838, 0.9999999999999999,
        0.632455532033676, 0.9999999999999998 },
    { 1.0, 1.0, 1.0, 1.0, 1.0, 0.2182178902359924, -0.6546536707079772, -1.0, 1.0,
        0.2182178902359924, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.2182178902359924,
        -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 0.5773502691896258, 0.2182178902359924, 1.0,
        1.0, 1.0, 1.0, 1.0, 0.2182178902359924, -0.6546536707079772, -1.0, 1.0,
        0.2182178902359924, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.2182178902359924,
        -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 0.5773502691896258, 0.2182178902359924, 1.0,
        1.0, 1.0, 1.0, 0.2182178902359924, -0.6546536707079772, -1.0, 1.0,
        0.2182178902359924, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.2182178902359924,
        -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 0.5773502691896258, 0.2182178902359924, 1.0,
        1.0, 1.0, 0.2182178902359924, -0.6546536707079772, -1.0, 1.0, 0.2182178902359924,
        1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.2182178902359924, -1.0, -1.0, -1.0, -1.0,
        1.0, -1.0, 0.5773502691896258, 0.2182178902359924, 1.0, 1.0, 0.2182178902359924,
        -0.6546536707079772, -1.0, 1.0, 0.2182178902359924, 1.0, -1.0, 1.0, 1.0, -1.0,
        -1.0, 1.0, 0.2182178902359924, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0,
        0.5773502691896258, 0.2182178902359924, 1.0, 1.0, -0.14285714285714285,
        -0.2182178902359924, 0.2182178902359924, -0.14285714285714285,
        0.2182178902359924, -0.2182178902359924, 0.2182178902359924, 0.2182178902359924,
        -0.2182178902359924, -0.2182178902359924, 0.2182178902359924,
        -0.14285714285714285, -0.2182178902359924, -0.2182178902359924,
        -0.2182178902359924, -0.2182178902359924, 0.2182178902359924,
        -0.2182178902359924, -0.3779644730092272, -0.14285714285714285,
        0.2182178902359924, 1.0, 0.6546536707079772, -0.6546536707079772,
        -0.14285714285714285, -0.6546536707079772, 0.6546536707079772,
        -0.6546536707079772, -0.6546536707079772, 0.6546536707079772, 0.6546536707079772,
        -0.6546536707079772, -0.14285714285714285, 0.6546536707079772,
        0.6546536707079772, 0.6546536707079772, 0.6546536707079772, -0.6546536707079772,
        0.6546536707079772, -0.3779644730092272, -0.14285714285714285,
        -0.6546536707079772, 1.0, -1.0, -0.2182178902359924, -1.0, 1.0, -1.0, -1.0, 1.0,
        1.0, -1.0, -0.2182178902359924, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0,
        -0.5773502691896258, -0.2182178902359924, -1.0, 1.0, 0.2182178902359924, 1.0,
        -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.2182178902359924, -1.0, -1.0, -1.0, -1.0, 1.0,
        -1.0, 0.5773502691896258, 0.2182178902359924, 1.0, 1.0, 0.2182178902359924,
        -0.2182178902359924, 0.2182178902359924, 0.2182178902359924, -0.2182178902359924,
        -0.2182178902359924, 0.2182178902359924, -0.14285714285714285,
        -0.2182178902359924, -0.2182178902359924, -0.2182178902359924,
        -0.2182178902359924, 0.2182178902359924, -0.2182178902359924,
        -0.3779644730092272, -0.14285714285714285, 0.2182178902359924, 1.0, -1.0, 1.0,
        1.0, -1.0, -1.0, 1.0, 0.2182178902359924, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0,
        0.5773502691896258, 0.2182178902359924, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0,
        -0.2182178902359924, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -0.5773502691896258,
        -0.2182178902359924, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 0.2182178902359924, -1.0,
        -1.0, -1.0, -1.0, 1.0, -1.0, 0.5773502691896258, 0.2182178902359924, 1.0, 1.0,
        -1.0, -1.0, 1.0, 0.2182178902359924, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0,
        0.5773502691896258, 0.2182178902359924, 1.0, 1.0, 1.0, -1.0, -0.2182178902359924,
        1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -0.5773502691896258, -0.2182178902359924, -1.0,
        1.0, -1.0, -0.2182178902359924, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0,
        -0.5773502691896258, -0.2182178902359924, -1.0, 1.0, 0.2182178902359924, -1.0,
        -1.0, -1.0, -1.0, 1.0, -1.0, 0.5773502691896258, 0.2182178902359924, 1.0, 1.0,
        -0.2182178902359924, -0.2182178902359924, -0.2182178902359924,
        -0.2182178902359924, 0.2182178902359924, -0.2182178902359924, 0.3779644730092272,
        -0.14285714285714285, 0.2182178902359924, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0,
        -0.5773502691896258, -0.2182178902359924, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0,
        -0.5773502691896258, -0.2182178902359924, -1.0, 1.0, 1.0, -1.0, 1.0,
        -0.5773502691896258, -0.2182178902359924, -1.0, 1.0, -1.0, 1.0,
        -0.5773502691896258, -0.2182178902359924, -1.0, 1.0, -1.0, 0.5773502691896258,
        0.2182178902359924, 1.0, 1.0, -0.5773502691896258, -0.2182178902359924, -1.0,
        1.0, 0.3779644730092272, 0.5773502691896258, 1.0, 0.2182178902359924, 1.0 },
    { 1.0, 1.0, 1.0 },
    { 0.18377223398316206, -0.12212786219416509, 0.2609542781331212 },
    { 0.24999999999999997, 0.0, 0.24999999999999997 },
    { 1.0, 0.0, 1.0 },
    { 1.0, -1.0, 1.0 },
    { 0.2609542781331212, 0.15896939941870186, 0.18377223398316206 },
};

double Dz_truth[11][1024] = {
    { 0.024691358024691357, -0.012345679012345678, 0.024691358024691357 },
    { 0.006944444444444444, 0.0, 0.0, 0.0, 0.0, 0.0 },
    /* { 0.046875, 0.0234375, 0.046875, -0.01171875, 0.0, 0.013671875, 0.0, 0.046875, */
    /*     0.001953125, 0.018663194444444444, 0.0234375, -0.009765625, 0.0, */
    /*     0.005533854166666666, 0.0, 0.0234375, 0.001627604166666666, 0.046875, */
    /*     -0.01171875, 0.0, 0.013671875, 0.0, 0.046875, 0.001953125, 0.0615234375, 0.0,
     */
    /*     -0.00732421875, 0.0, -0.01171875, 0.00146484375, 0.0, 0.0, 0.0, 0.0, 0.0, */
    /*     0.010498046875, 0.0, 0.013671875, -0.000732421875, 0.0, 0.0, 0.0, 0.046875, */
    /*     0.001953125, 0.001708984375 }, */
    { 0.006858710562414267, 0.004115226337448559, -0.00823045267489712,
        -0.00823045267489712, 0.024691358024691357, -0.012345679012345678,
        -0.012345679012345671, 0.024691358024691353, 0.024691358024691357,
        0.061728395061728406 },
    { 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.03515625,
        0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, -0.01171875, 0.03515625, 0.046875, 0.046875, -0.01171875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875,
        -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0,
        -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.03515625,
        0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875, 0.046875,
        -0.01171875, 0.03515625, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875,
        0.046875, -0.01171875, 0.03515625, 0.046875, 0.046875, -0.01171875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, -0.01171875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875,
        0.046875, 0.0615234375, -0.0087890625, -0.01171875, -0.01171875, -0.0087890625,
        -0.01171875, -0.01171875, -0.01171875, -0.01171875, -0.01171875, -0.01171875,
        -0.01171875, -0.0087890625, -0.01171875, -0.01171875, -0.01171875, -0.01171875,
        -0.01171875, -0.01171875, 0.0, -0.0087890625, -0.01171875, 0.0615234375,
        0.03515625, 0.03515625, -0.0087890625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, -0.0087890625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0, -0.0087890625,
        0.03515625, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875,
        -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.0, -0.01171875, 0.046875, 0.0615234375, -0.01171875, -0.01171875,
        -0.01171875, -0.01171875, -0.01171875, -0.01171875, -0.01171875, -0.0087890625,
        -0.01171875, -0.01171875, -0.01171875, -0.01171875, -0.01171875, -0.01171875,
        0.0, -0.0087890625, -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875, 0.046875,
        0.046875, -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.0, -0.01171875, 0.046875, 0.046875, 0.046875, -0.01171875, 0.046875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875,
        0.046875, -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875,
        0.046875, 0.0, -0.01171875, 0.046875, 0.0615234375, -0.01171875, -0.01171875,
        -0.01171875, -0.01171875, -0.01171875, -0.01171875, 0.0, -0.0087890625,
        -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0,
        -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0,
        -0.01171875, 0.046875, 0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875,
        0.046875, 0.046875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.046875,
        0.046875, 0.0, -0.01171875, 0.046875, 0.046875, 0.0, -0.01171875, 0.046875, 0.0,
        0.0, 0.0, 0.0615234375, -0.01171875, 0.046875 },
    { 0.0, 0.0, 0.0 },
    { 0.0033870175616860566, 0.003387017561686057, 0.003387017561686057 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.003387017561686057, -0.00027096140493488457, 0.0033870175616860566 },
};

double pi2_truth[11][1024] = {
    { 0.04938271604938272, 0.04938271604938272, 0.04938271604938272 },
    { 0.043402777777777776, 0.05208333333333333, 0.05208333333333333, 0.0625, 0.0625,
        0.0625 },
    /* { 0.03515625, 0.033203125, 0.03515625, 0.0205078125, 0.0, 0.02490234375, 0.0, */
    /*     0.03515625, 0.03369140625, 0.03135850694444444, 0.033203125, */
    /*     0.019368489583333336, 0.0, 0.023518880208333332, 0.0, 0.033203125, */
    /*     0.031819661458333336, 0.03515625, 0.0205078125, 0.0, 0.02490234375, 0.0, */
    /*     0.03515625, 0.03369140625, 0.011962890625, 0.0, 0.0145263671875, 0.0, */
    /*     0.0205078125, 0.0196533203125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01763916015625, 0.0,
     */
    /*     0.02490234375, 0.02386474609375, 0.0, 0.0, 0.0, 0.03515625, 0.03369140625, */
    /*     0.03228759765625 }, */
    { 0.041495198902606306, 0.04526748971193416, 0.04526748971193416,
        0.028292181069958854, 0.04938271604938272, 0.04938271604938272,
        0.030864197530864203, 0.04938271604938272, 0.030864197530864203,
        0.01929012345679012 },
    { 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0205078125,
        0.0205078125, 0.03515625, 0.03515625, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0205078125, 0.0205078125, 0.03515625, 0.03515625, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0205078125, 0.0205078125, 0.03515625, 0.03515625, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.0205078125, 0.0205078125, 0.03515625, 0.03515625, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.03515625,
        0.0205078125, 0.0205078125, 0.03515625, 0.03515625, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.011962890625, 0.011962890625,
        0.0205078125, 0.0205078125, 0.011962890625, 0.0205078125, 0.0205078125,
        0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125,
        0.011962890625, 0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125,
        0.0205078125, 0.0205078125, 0.02734375, 0.011962890625, 0.0205078125,
        0.011962890625, 0.0205078125, 0.0205078125, 0.011962890625, 0.0205078125,
        0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125,
        0.0205078125, 0.011962890625, 0.0205078125, 0.0205078125, 0.0205078125,
        0.0205078125, 0.0205078125, 0.0205078125, 0.02734375, 0.011962890625,
        0.0205078125, 0.03515625, 0.03515625, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875,
        0.0205078125, 0.03515625, 0.03515625, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875,
        0.0205078125, 0.03515625, 0.011962890625, 0.0205078125, 0.0205078125,
        0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125,
        0.011962890625, 0.0205078125, 0.0205078125, 0.0205078125, 0.0205078125,
        0.0205078125, 0.0205078125, 0.02734375, 0.011962890625, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875,
        0.0205078125, 0.03515625, 0.03515625, 0.0205078125, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875, 0.0205078125,
        0.03515625, 0.011962890625, 0.0205078125, 0.0205078125, 0.0205078125,
        0.0205078125, 0.0205078125, 0.0205078125, 0.02734375, 0.011962890625,
        0.0205078125, 0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.03515625, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.046875, 0.0205078125, 0.03515625,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.03515625, 0.046875, 0.0205078125,
        0.03515625, 0.03515625, 0.03515625, 0.046875, 0.0205078125, 0.03515625,
        0.03515625, 0.046875, 0.0205078125, 0.03515625, 0.0625, 0.02734375, 0.046875,
        0.011962890625, 0.0205078125, 0.03515625 },
    { 0.0625, 0.0625, 0.0625 },
    { 0.04579247743399549, 0.04579247743399549, 0.0457924774339955 },
    { 0.04938271604938272, 0.04938271604938272, 0.04938271604938272 },
    { 0.0625, 0.0625, 0.0625 },
    { 0.0625, 0.0625, 0.0625 },
    { 0.0457924774339955, 0.04579247743399549, 0.04579247743399549 },
};