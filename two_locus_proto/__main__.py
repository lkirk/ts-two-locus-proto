import argparse
import sys

import tskit

from two_locus_proto.site import two_site_general_stat
from two_locus_proto.summary_functions import *


NORM_METHOD = {
    D: "total",
    D_prime: "hap_weighted",
    D2: "total",
    Dz: "total",
    pi2: "total",
    r: "total",
    r2: "hap_weighted",
}

POLARIZATION = {
    D: True,
    D_prime: True,
    D2: False,
    Dz: False,
    pi2: False,
    r: True,
    r2: False,
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, help="tree sequence file to parse")
    parser.add_argument("--summary-func", required=True, help="summary func to apply to tree sequence file")
    parser.add_argument("--debug", default=False, action="store_true", help="debug mode")
    parser.add_argument("--print-weights", default=False, action="store_true", help="print weights to stdout")
    return parser.parse_args()


def main():
    args = parse_args()
    ts = tskit.load(args.tree)
    summary_func = globals()[args.summary_func]
    two_site_general_stat(
        ts=ts,
        summary_func=summary_func,
        norm_method=NORM_METHOD[summary_func],
        polarized=POLARIZATION[summary_func],
        debug=args.debug,
        print_weights=args.print_weights,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
