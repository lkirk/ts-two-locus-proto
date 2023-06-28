import argparse
import sys

import tskit

from two_locus_proto.site import two_site_general_stat
from two_locus_proto.summary_functions import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', required=True, help='tree sequence file to parse')
    parser.add_argument('--summary-func', required=True, help='summary func to apply to tree sequence file')
    parser.add_argument('--polarized', default=False, action='store_true', help='polarize results')
    return parser.parse_args()


def main():
    args = parse_args()
    ts = tskit.load(args.tree)
    summary_func = globals()[args.summary_func]
    two_site_general_stat(ts, summary_func, 'total', args.polarized)
    return 0


if __name__ == '__main__':
    sys.exit(main())
