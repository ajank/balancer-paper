#! /usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import print_function
import argparse
import os
import re
import sys
from subprocess import check_output, call
from tempfile import NamedTemporaryFile
import pysam
import gzip
from datetime import date
from sortedcontainers import SortedListWithKey
VERSION = 0.3

# external tools:
samtools = '/g/korbel/meiers/tools/bin/samtools'
tabix    = '/g/software/bin/tabix-1.2.1'

def iterVCF(file, mode='no-multi-allelic-variants'):
    """this is a very rudimentary vcf parser - use with care"""
    try:
        f=gzip.GzipFile(fileobj=file) # a first line is ignored (should be header)
        f.next() # reads a few bytes until it detects whether it's gzipped or not
    except IOError:
        f=file # fallback
        f.next() # read to next line
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        s = line.split()
        if mode == 'no-multi-allelic-variants' and (',' in s[3] or ',' in s[4]):
            continue
        info = dict() if len(s)<8 else \
               dict([tuple(kv.split('=')[:2]) if '=' in kv else (kv,True) for kv in s[7].split(';')])
        yield (s[0],          #CHROM  
               int(s[1]) - 1, #POS --> make 1-based coordinate 0-based
               s[2],          #ID
               s[3],          #REF
               s[4],          #ALT     
               float(s[5]),   #QUAL
               s[6],          #FILTER
               info)          #INFO


if __name__ == "__main__":
    # arguments
    parser = argparse.ArgumentParser(description='Separates reads based on SNVs '
            'into REF and ALT reads. Annotation is performed in the CO field of the '
            'bam file. SNPs with base quality < Q are ignored. If a read shows an '
            'allele different from the two expected ones, UNKNOWN is stated.')
    parser.add_argument('-v', '--variants', dest='f_vcf',  metavar='VCF.gz',
            type=argparse.FileType('rb'), required=True,
            help='Tabix indexed VCF.gz file containing SNVs/indels.')
    parser.add_argument('-b', '--bam', metavar='BAM', dest='f_bam',
            required=True, help='BAM file with reads to be separated.')
    parser.add_argument('-o', '--out', metavar='BAM', dest='f_out',
            required=True, help='Output BAM file.')
    parser.add_argument('-q', '--qual', metavar='Q', dest='i_qual', choices=range(0,255),
            required=False, default=10, help='Minimum base quality. [0-255], default: 10).')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(VERSION))
    args = parser.parse_args()

    # Look up table (sorted list) for VCF entries
    print ("Parsing VCF...", file=sys.stderr)
    vcfList = SortedListWithKey(iterVCF(args.f_vcf),
                                key=lambda x: (x[0],x[1]))

    global_count=0
    print ("Parsing BAM...", file=sys.stderr)
    with pysam.AlignmentFile(args.f_bam,'rb') as bam:
        
        header = bam.header
        if not 'PG' in header:
            header['PG'] = []
        header['PG'].append(dict(ID = 'separator_' + str(date.today()),
                                 PN = 'separator.py',
                                 VN = str(VERSION),
                                 CL = " ".join(sys.argv)))

        with pysam.AlignmentFile(args.f_out, 'wb', header=header) as bamOut:
            
            # go through reads sorted by read name
            rdname = None
            annotations = []
            readBuffer = []
            counter = dict(ref = 0, alt = 0, unknown = 0, bam = 0)
            for read in bam:
                
                # after all alignments of a read, write them to file.
                if rdname != read.query_name:
                    for r in readBuffer:
                        for a in annotations:
                            r.set_tag('CO', a, replace=False)
                        bamOut.write(r)
                    annotations = []
                    readBuffer = []
                rdname = read.query_name

                counter["bam"] +=1
                if counter["bam"] % 1000000 == 0:
                    print("  processed {}M bam entries".format(counter["bam"]/1000000), file=sys.stderr)

                if not read.is_unmapped:
                    chrom = bam.getrname(read.reference_id)
                    _min, _max = read.reference_start, read.reference_end

                    for snp in vcfList.irange((chrom, _min), (chrom, _max)):
                        # fetch [q]uery position and [r]eference position. It's over-complex, 
                        # but with this method I am at least sure I get exactly the position I want
                        aligned_positions = SortedListWithKey(read.get_aligned_pairs(matches_only=True),
                                                              key=lambda x: x[1])
                        try: 
                            qp,rp = list(aligned_positions.irange_key(snp[1],snp[1]))[0]
                        except: # pos not aligned.
                            continue

                        rdflag = ('SA' if read.flag&2048 else '') + ('2' if read.flag&128 else'1')

                        if read.query_qualities[qp] < args.i_qual:
                            continue

                        if len(snp[3]) > len(snp[4]): # deletion
                            print ("Ignoring non-SNV variants!", snp, file=sys.stderr)
                            #if read.query_sequence[qp:qp+len(snp[3])] == snp[3]:
                            #    annotations.append('REF,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                            #elif read.query_sequence[qp:qp+len(snp[4])] == snp[4]:
                            #    annotations.append('ALT,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                        elif len(snp[4]) > len(snp[3]): # insertion
                            print ("Ignoring non-SNV variants!", snp, file=sys.stderr)
                            #if read.query_sequence[qp:qp+len(snp[4])] == snp[4]:
                            #    annotations.append('ALT,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                            #elif read.query_sequence[qp:qp+len(snp[3])] == snp[3]:
                            #    annotations.append('REF,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                        else:
                            if read.query_sequence[qp:qp+len(snp[3])] == snp[3]:
                                annotations.append('REF,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                                counter["ref"] +=1
                            elif read.query_sequence[qp:qp+len(snp[4])] == snp[4]:
                                annotations.append('ALT,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                                counter["alt"] +=1
                            else:
                                annotations.append('UNKNOWN,read{},{},{},{},{}'.format(rdflag,snp[1],snp[3],snp[4],read.query_qualities[qp]))
                                counter["unknown"] +=1

                # also output unmapped reads!
                readBuffer.append(read)

            # after last read
            for r in readBuffer:
                for a in annotations:
                    r.set_tag('CO', a, replace=False)
                bamOut.write(read)
                annotations = []
                readBuffer = []

            print("Finished.", file=sys.stderr)
            for x in counter:
                print (x,counter[x],sep="\t",file=sys.stderr)
