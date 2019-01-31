#!/usr/bin/env python
## Martijn Derks
## Pipeline to gwt coverage information from SV output.

import commands
from collections import OrderedDict
import sys
import gzip
import vcf
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser( description='Add coverage info to VCF file')
parser.add_argument("-v", "--vcf_file", help="VCF_file (compressed) with only indels", nargs=1)

class Coverage():

	def coverage_info(self):
		self.cov_dic = {}
		cov_table_file = open("Table_coverage.txt")
		for sample in cov_table_file:
			samplename, cov = sample.strip().split(" ")
			self.cov_dic[samplename]=float(cov.replace("X",""))
			
	def add_coverage_info_VCF(self, vcf_file):
		bampath="/lustre/nobackup/WUR/ABGC/shared/great_tit/BAM_30birds/"
		vcf_reader = vcf.Reader(filename = vcf_file, strict_whitespace=True)
		vcf_writer = vcf.Writer(open("Filtered."+vcf_file, 'w'), vcf_reader)
		samples = vcf_reader.samples
		for record in vcf_reader:
			if float(record.QUAL) < 20:
				continue
			type = record.INFO['SVTYPE']
			if type == "BND":
				vcf_writer.write_record(record)
				continue
			chr=record.CHROM
			start = int(record.POS)
			end = int(record.INFO['END'])
			length = end-start
			for i in range(len(record.samples)):
				sample = record.samples[i]
				samplename=sample.sample.split("-")[-1]
				region_cov = commands.getoutput("samtools depth -a -r "+str(chr)+":"+str(start)+"-"+str(end)+" "+bampath+"chr1A."+samplename+".dedup.real.bam |  awk '{sum+=$3} END { print (sum+1)/NR}'")
				if "fatal" in region_cov:
					region_cov=0.01
				
				ratio = (float(region_cov)/self.cov_dic[samplename])
				sample.add_field("CV",[float(region_cov), self.cov_dic[samplename], round(ratio,3)])
				record.samples[i] = sample
			vcf_writer.write_record(record)

if __name__ == '__main__':
	args = parser.parse_args()
	vcf_input = args.vcf_file[0]
	C=Coverage()
	C.coverage_info()
	C.add_coverage_info_VCF(vcf_input)				
