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
from sortedcontainers import SortedListWithKey


# external tools:
samtools = '/g/korbel/meiers/tools/bin/samtools'
tabix    = '/g/software/bin/tabix-1.2.1'

def iterVCF(file):
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

        # Ignore multi-allelic variants
        if (',' in s[3] or ',' in s[4]):
            continue

        # Ignore non-SNPs
        if len(s[3]) != 1 or len(s[4]) != 1:
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

def write_reads(rdBuf, anno, f_dict):
    # determine which file to write to
    if not anno:
        bamOut = f_dict["amb"]
    elif len(set([x[0] for x in anno])) > 1:
        bamOut = f_dict["err"]
    elif anno[0][0] == "REF":
        bamOut = f_dict["ref"]
    else:
        bamOut = f_dict["alt"]

    anno_string = ",".join(["snp:{}".format(x[1]) for x in anno])

    # write down reads
    for r in rdBuf:
        r.set_tag('CO', anno_string, replace=False)
        bamOut.write(r)


if __name__ == "__main__":
    # arguments
    parser = argparse.ArgumentParser(description='Separate reads based on SNVs/small '
            'indels into REF and ALT reads. Annotation is performed in the CO field of the '
            'bam file. If neither REF nor ALT apply, UNKNOWN is stated.')
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

    if args.f_out.endswith(".bam"):
        args.f_out = args.f_out[:-4]

    global_count=0
    print ("Parsing BAM ...")
    with pysam.AlignmentFile(args.f_bam,'rb') as bam:
        with pysam.AlignmentFile(args.f_out + ".alt.bam", 'wb', template=bam) as bam_alt, \
             pysam.AlignmentFile(args.f_out + ".ref.bam", 'wb', template=bam) as bam_ref, \
             pysam.AlignmentFile(args.f_out + ".err.bam", 'wb', template=bam) as bam_err, \
             pysam.AlignmentFile(args.f_out + ".amb.bam", 'wb', template=bam) as bam_amb:
            
            f_out = dict(alt=bam_alt, ref=bam_ref, err=bam_err, amb=bam_amb)
            # go through reads sorted by read name
            rdname = None
            annotations = []
            readBuffer = []
            for read in bam:
                
                # after all alignments of a read, write them to file.
                if rdname != read.query_name:
                    write_reads(readBuffer, annotations, f_out)
                    annotations = []
                    readBuffer = []
                rdname = read.query_name

                # Ignore unmapped and suppl. reads
                if not read.is_unmapped and not read.flag&2048:

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
                        if read.query_sequence[qp:qp+len(snp[3])] == snp[3]:
                            annotations.append( ('REF', rp, qp, '2' if read.flag&128 else'1'))
                        if read.query_sequence[qp:qp+len(snp[4])] == snp[4]:
                            annotations.append( ('ALT', rp, qp, '2' if read.flag&128 else'1'))
                
                # also output unmapped reads!
                readBuffer.append(read)

            # after last read
            write_reads(readBuffer, annotations, f_out)
            annotations = []
            readBuffer = []

