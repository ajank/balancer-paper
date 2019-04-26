#!/usr/bin/env python

# Script by Sascha Meiers

from __future__ import print_function
import argparse
import pysam
import sys
VERSION = 0.1

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Add flag to all entries in a bam file. Output written to stdout')
    parser.add_argument('bam', help='Input BAM')
    parser.add_argument('flag', type=int, help='SAM flag to be added')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(VERSION))
    args = parser.parse_args()

    with pysam.AlignmentFile(args.bam,'rb') as bam:
        header = bam.header
        if not 'PG' in header:
            header['PG'] = []
        header['PG'].append(dict(ID = 'addFlag_'+str(args.flag), 
                                 PN = 'addFlag.py', 
                                 VN = str(VERSION),
                                 CL = " ".join(sys.argv)))

        with pysam.AlignmentFile("-", "wb", header=header) as out:
            for x in bam:
                x.flag |= args.flag
                out.write(x)