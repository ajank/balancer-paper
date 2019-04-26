#! /usr/bin/env python2.7
# -*- coding: latin-1 -*-

# Script by Sascha Meiers

from __future__ import print_function
import argparse
import os
import re
import sys
from subprocess import check_output, call
from tempfile import NamedTemporaryFile
import pysam
import gzip
from sortedcontainers import SortedListWithKey


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
    parser = argparse.ArgumentParser(description='simplified version of read separation. '
            'Output into 1 file, only reads actually under the SNP are annotated. '
            'Annotations in VRG and BAL')
    parser.add_argument('-v', '--variants', dest='f_vcf',  metavar='VCF.gz',
            type=argparse.FileType('rb'), required=True,
            help='Tabix indexed VCF.gz file containing SNVs/indels.')
    parser.add_argument('-b', '--bam', metavar='BAM', dest='f_bam',
            required=True, help='BAM file with reads to be separated.')
    parser.add_argument('-o', '--out', metavar='BAM', dest='f_out',
            required=True, help='Output BAM file.')
    args = parser.parse_args()

    # Look up table (sorted list) for VCF entries
    print ("Parsing VCF ...")
    vcfList = SortedListWithKey(iterVCF(args.f_vcf),
                                key=lambda x: (x[0],x[1]))

    global_count=0
    print ("Parsing BAM ...")
    with pysam.AlignmentFile(args.f_bam,'rb') as bam:
        with pysam.AlignmentFile(args.f_out, 'wb', template=bam) as bamOut:
            
            # go through reads, order does not matter
            for read in bam:
                rdname = read.query_name
                a = []

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
                    
                        # Annotation: allele, ref pos, query pos, read 1/2
                        # Please only inupt SNPs, I don't know whether this works on indels, too. Likely it does not.
                        if read.query_sequence[qp:qp+len(snp[3])] == snp[3]:
                            a.append('VRG' + '|' + str(rp) + '|' + str(qp) + '|' + snp[3] + '|' +  snp[4])
                        if read.query_sequence[qp:qp+len(snp[4])] == snp[4]:
                            a.append('BAL' + '|' + str(rp) + '|' + str(qp) + '|' + snp[3] + '|' +  snp[4])
                    if a:
                        read.set_tag('CO', ",".join(a), replace=False)
                
                # also output unmapped reads!
                bamOut.write(read)


