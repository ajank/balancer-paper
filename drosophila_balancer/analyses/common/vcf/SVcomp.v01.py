from __future__ import print_function
import pybedtools as pybed
import re
import vcf
import gzip
import sys
import argparse
import itertools as it
import os


def get_intersection(A, B, **args):
	'''Given 2 BED files ('wb' format) that carry their integer key in 
	column 4, extract the matching keys'''

	# Prevent error on empty bed files
	if not A or not B:
		return []

	# Bedtools intersect
	isc = A.intersect(B, **args)
	
	# Read out keys
	matches = []
	for l in isc:
		s = str(l).strip().split()
		matches.append( (int(s[3]), int(s[7])))
	# clean up open file
	os.remove(isc.fn)
	del isc

	return matches


if __name__ == "__main__":


	def format_input_file(s):
		# Return name and VCF file object.
		if s == '-':
			return ("<stdin>", vcf.Reader(sys.stdin, 'r'))
		else:
			try:
				return (os.path.basename(s), 
					    vcf.Reader(open(s), 'rb', compressed=True) if s.endswith('.gz') else vcf.Reader(open(s), 'r', compressed=False))
			except SyntaxError:
				raise
			except:
				msg = "Not a valid VCF file. '{0}'.".format(s)
				raise argparse.ArgumentTypeError(msg)



	parser = argparse.ArgumentParser(description='compare two or more SV VCF files')
	parser.add_argument('sv_files', metavar='FILE', type=format_input_file, nargs='+',
			help='input files containing SV calls [vcf, vcf.gz]')
	parser.add_argument('-v', '--version', action='version',                    
			version='%(prog)s (version 0.1)')
	parser.add_argument('-r', '--recpr', type=float, metavar='F',
	        help='require reciprocal overlap of this fraction. F in (0,1]')
	args = parser.parse_args()


	print ("Loading files into memory...")
	data = []

	for name,f in args.sv_files:
		vcfdata = dict(records=[])
		vcfdata['name'] = name

		# Read VCF
		vcf_bed = []
		for i,rec in enumerate(f):

			# get the samples (only once)
			if i==0:
				vcfdata['samples'] = set([call.sample for call in rec.samples])

			# Ignore non-SV lines
			if type(rec.ALT[0]) != type(vcf.model._SV(None)) or 'END' not in rec.INFO or type(rec.INFO['END']) != type(1):
				print ("Ignoring non-SV entry in VCF {}. Could also be a header with VCF header.".format(name), file=sys.stderr)
				continue
			# store the records in memory
			vcfdata['records'].append(rec)
			
			# make a bed file
			vcf_bed.append( (rec.CHROM, rec.POS, rec.INFO['END'], i) )

		if not 'samples' in vcfdata:
			vcfdata['samples'] = set()

		vcfdata['bed'] = pybed.BedTool(vcf_bed)
		data.append(vcfdata)

	# Prepare intesection parameters
	isc_args= dict(wb=True)
	if args.recpr:
		if args.recpr>1 or args.recpr<0:
			print("-r must be within (0,1]", sys.stderr)
			sys.exit(1)
		isc_args['r']=True
		isc_args['f']=args.recpr
		print("Required reciprocal overlap of {}".format(args.recpr))
	print()


	print ("Single VCF files:")
	for A in data:
		print ("  -", A['name'])
		print (" "*5, "{} calls in total".format(len(A['records'])))
		m = get_intersection(A['bed'], A['bed'], **isc_args)
		print (" "*5, "{} pairwise overlaps among those".format(len([x for x in m if x[0]!=x[1]])))


	print()
	print ("Pairwise comparisons:")
	for A,B in it.combinations(data, 2):
		print ("  -", "{} vs. {}".format(A['name'], B['name']))
		m = get_intersection(A['bed'], B['bed'], **isc_args)
		print (" "*5, "{} overlap".format(len(m)))

		# Looking at the sample overlap
		cmn_smp = A['samples'].intersection(B['samples'])
		print (" "*5, "They share {} samples: {}{}".format(
				len(cmn_smp), 
				", ".join(list(cmn_smp)[:5]),
				",..." if len(cmn_smp)>5 else ""))

		# Comparing genotypes:
		c_ge1 = c_all = 0
		for a,b in m:
			recA = A['records'][a]
			recB = B['records'][b]
			cmp_GT = [recA.genotype(s)['GT'] == recB.genotype(s)['GT'] for s in cmn_smp]
			if any(cmp_GT):
				c_ge1 += 1
				if all(cmp_GT):
					c_all += 1
		print (" "*5, "{} GTs agree in at least one sample".format(c_ge1))
		print (" "*5, "{} GTs agree in all samples".format(c_ge1))
			

