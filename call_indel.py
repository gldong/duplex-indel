#!/usr/bin/env python

# Performs double-stranded and single-stranded Indel calling from Tn5 duplex sequencing data. 
# NOTE: 
# This code is inspired by and adapted from https://github.com/lh3/lianti (21a15c8) which was
# originally written by Heng Li to perform duplex SNV calling. We thank Heng Li for allowing us 
# to adapt this code for indel calling in this repo.

import argparse
import os
import re
import sys
import pysam
from dataclasses import dataclass, field
import math


def parse_args():
	parser = argparse.ArgumentParser(
		description="Call somatic Indels and output an annotated VCF.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	# General options
	parser.add_argument('-b', type=int, default=1, help='Number of bulk samples')
	parser.add_argument('-hs', dest='fn_hap', help='File listing haploid samples')
	parser.add_argument('-H', action='store_true', help='Mark all samples as haploid')
	parser.add_argument('-q', type=int, default=50, help='Min mapping quality')
	parser.add_argument('-e', dest='fn_excl', help='File listing samples to exclude')
	parser.add_argument('-v', dest='fn_var', help='VCF file with sites to exclude')
	parser.add_argument('-r', dest='fn_rep', help='Replicate sample mapping file')
	parser.add_argument('-F', action='store_true', help='Print mutations filtered by -w and -v')
	parser.add_argument('-u', action='store_true', help='Process autosomes only')
	parser.add_argument('-o', required=True, help='Output VCF file')

	# Duplex-specific filters
	parser.add_argument('-a', type=int, default=5, help='Min ALT read depth in total to call a double-stranded mutation')
	parser.add_argument('-s', type=int, default=2, help='Min ALT read depth per strand to call a double-stranded mutation')
	parser.add_argument('-S', type=int, default=4, help='Min strand depth at candidate single-stranded sites')
	parser.add_argument('-B', type=float, default=0.2, help='Min ALT allele balance')
	parser.add_argument('-j', type=int, default=2, help='Min allele depth to call joint mutations')
	parser.add_argument('-J', type=int, default=1, help='Min allele depth on both strands for joint mutations')
	parser.add_argument('-l', type=int, default=1, help='Max conflicting duplex reads')
	parser.add_argument('-L', type=int, default=10, help='Min distance towards read end')
	parser.add_argument('-w', type=int, default=100, help='Window size to filter clustered mutations')
	parser.add_argument('-R', type=int, default=0, help='Max REF read depth at mutation site')
	parser.add_argument('-T', type=int, choices=[1, 2, 3], default=3, help='Filtering stringency level (1=most strict, 3=least)')

	# Bulk-specific filters
	parser.add_argument('-D', type=int, default=20, help='Min bulk read depth')
	parser.add_argument('-A', type=int, default=8, help='Min bulk ALT read depth to call a het')
	parser.add_argument('-m', type=int, default=0, help='Max bulk ALT read depth to call a mutation')
	parser.add_argument('-P', action='store_true', help='The bulk is haploid')

	# Positional argument: input VCF file
	parser.add_argument('vcf', help='Joint VCF input file')

	return parser.parse_args()


def read_list(fn):
	h = {}
	if not fn or fn.strip() == "":
		return h
	with open(fn, 'r') as f:
		for line in f:
			if line.strip() == "":
				continue
			key = line.strip().split('\t')[0]
			h[key] = 1
	return h


@dataclass
class SampleMeta:
	name: str
	ploidy: int
	col: int
	ado: list  # [ref, alt] allele dropout
	fn: int
	mut: int
	dmg: int
	dmg_fp: int
	dmg_fn: list  # [ref, alt]
	calls: list


def write_output_vcf_header(vcf_in, output_path, output_sample_name):
	# Create a new empty header
	header_out = pysam.VariantHeader()
	# Copy fileformat, contig and reference lines
	for rec in vcf_in.header.records:
		if rec.key in ("fileformat", "contig", "reference"):
			header_out.add_line(str(rec))
	# Add INFO and FORMAT fields
	header_out.add_line('##INFO=<ID=BDP,Number=1,Type=String,Description="Bulk depths of called alleles">')
	header_out.add_line('##INFO=<ID=VT,Number=1,Type=String,Description="Variant type: DS/SS/JOINT">')
	header_out.add_line('##INFO=<ID=NBC,Number=1,Type=Integer,Description="Number of BC bams with variant called">')
	header_out.add_line('##FORMAT=<ID=BC,Number=1,Type=String,Description="BC bam names with variant called">')
	header_out.add_line('##FORMAT=<ID=ADF,Number=1,Type=String,Description="Allelic depths on the forward strand per BC bam">')
	header_out.add_line('##FORMAT=<ID=ADR,Number=1,Type=String,Description="Allelic depths on the reverse strand per BC bam">')
	header_out.add_line('##FORMAT=<ID=NTN,Number=1,Type=String,Description="Number of unique Tn5 sites per BC bam">')
	header_out.add_line('##FORMAT=<ID=TN,Number=1,Type=String,Description="Left|Right end positions of Tn5 sites per BC bam">')
	header_out.add_line('##FORMAT=<ID=MGL,Number=1,Type=String,Description="Left end position of merged window per BC bam">')
	header_out.add_line('##FORMAT=<ID=MGR,Number=1,Type=String,Description="Right end position of merged window per BC bam">')
	# Add command line
	cmdline = " ".join(sys.argv)
	header_out.add_line(f'##cmd="{cmdline}"')
	# Add sample name
	header_out.add_sample(output_sample_name)
	# Open output file and write header
	vcf_out = pysam.VariantFile(output_path, "w", header=header_out)
	return vcf_out


def parse_vcf_header(args, sample_excl, sample_hap, rep_str):
	vcf_in = pysam.VariantFile(args.vcf)
	samples = list(vcf_in.header.samples)

	col2sample = {} # dict mapping column number to sample
	sample_meta = []

	for idx, sample in enumerate(samples[args.b:]):  # skip bulk samples
		s1 = sample
		s2 = sample.replace('.bam', '')

		if s1 in sample_excl or s2 in sample_excl:
			continue
		if rep_str and (s1 in rep_str or s2 in rep_str):
			continue

		ploidy = 1 if args.H or s1 in sample_hap or s2 in sample_hap else 2
		name = s2
		col = idx + args.b  # actual column index in VCF

		sample_meta.append(SampleMeta(
			name=name,
			ploidy=ploidy,
			col=col,
			ado=[0, 0],
			fn=0,
			mut=0,
			dmg=0,
			dmg_fp=0,
			dmg_fn=[0, 0],
			calls=[]
		))

	# Build col2sample mapping
	for i, cm in enumerate(sample_meta):
		col2sample[cm.col] = i

	# Build replicate ID map
	rep_id = {i: i for i in range(len(samples))}  # default: identity

	if rep_str:
		sample2id = {s: i for i, s in enumerate(samples)}
		for i, s in enumerate(samples):
			s2 = s.replace('.bam', '')
			rep_key = rep_str.get(s, rep_str.get(s2))
			if rep_key and rep_key in sample2id:
				rep_id[i] = sample2id[rep_key]

	return vcf_in, samples, col2sample, sample_meta, rep_id


def load_var_filter(fn_var):
	var_map = set()
	if not fn_var:
		return var_map
	with open(fn_var) as f:
		for line in f:
			t = line.strip().split('\t')
			if t[0][0] == "#":
				continue
			var_map.add(t[0]+":"+t[1])
	return var_map


def load_rep_map(fn_rep):
	rep_str = {}
	if not fn_rep:
		return rep_str
	with open(fn_rep) as f:
		for line in f:
			t = line.strip().split()
			if not t:
				continue
			canonical = t[0]
			for alias in t[1:]:
				rep_str[alias] = canonical
	return rep_str


def is_autosome(ctg):
	return re.match(r'^(chr)?([0-9]+)$', ctg) is not None


def parse_ad(sample_data):
	adf = sample_data['ADF']
	adr = sample_data['ADR']
	if adf is None or adr is None:
		raise RuntimeError("missing ADF or ADR in FORMAT")
	adf = [int(x) for x in adf]
	adr = [int(x) for x in adr]
	if len(adf) != len(adr):
		raise RuntimeError("Inconsistent VCF (ADF/ADR length mismatch)")
	ad = [adf[i] + adr[i] for i in range(len(adf))]
	dp_ref = ad[0] if ad else 0
	dp_alt = sum(ad[1:]) if len(ad) > 1 else 0
	return adf, adr, ad, dp_ref, dp_alt


def filter_germline_adj(x, germline_adj):
	# Filter ds/ss calls adjacent to germline indels
	for s,e in germline_adj:
		if x['pos'] > s and x['pos'] < e:
			return True
	return False


def print_calls(vcf_out, output_sample_name, x, sample_meta, show_flt, min_joint, min_joint_strand):
	bulk_ad = [0,0]
	bulk_alt = []
	sample_hit_joint = []
	sample_hit_ds = []
	sample_hit_ss = []
	sample_hit_flt = []

	# bulk info
	for b in x['bulks']:
		bulk_ad[0] += b['ad'][0] if len(b['ad'])>0 else 0
		bulk_ad[1] += b['ad'][1] if len(b['ad'])>1 else 0
		bulk_alt.append(str(b['ad'][1] if len(b['ad'])>1 else 0))
	if bulk_ad[1] != 0:
		bulk_ad[1] = "|".join(bulk_alt)

	# sample info
	for i, c in enumerate(x['BC_samples']):

		# double-stranded call (ds)
		if (not c['flt']) and c['alt'] and (not x['flt']):
			sample_meta[i].mut += 1
			sample_hit_ds.append({
				"BC": sample_meta[i].name,
				"ADF": ",".join(map(str,c['adf'])),
				"ADR": ",".join(map(str,c['adr'])),
				"NTN": str(c['count_ltpos']),
				"TN": str(c['uniq_ltpos']),
				"MGL": str(c['min_mg_start']) if math.isfinite(c['min_mg_start']) else ".",
				"MGR": str(c['max_mg_end']) if c['max_mg_end']>0 else "."
			})

		# single-stranded call (ss)
		if (not c['flt']) and c['dmg'] and (not x['flt']):
			sample_meta[i].dmg += 1
			sample_hit_ss.append({
				"BC": sample_meta[i].name,
				"ADF": ",".join(map(str,c['adf'])),
				"ADR": ",".join(map(str,c['adr'])),
				"NTN": str(c['count_ltpos']),
				"TN": str(c['uniq_ltpos']),
				"MGL": str(c['min_mg_start']) if math.isfinite(c['min_mg_start']) else ".",
				"MGR": str(c['max_mg_end']) if c['max_mg_end']>0 else "."
			})

		# joint call (joint)
		if (not x['flt']) and x['n_joint_alt'] >= 2:
			if c['ad'][1] >= min_joint and c['adf'][1] >= min_joint_strand and c['adr'][1] >= min_joint_strand:
				sample_hit_joint.append({
					"BC": sample_meta[i].name,
					"ADF": ",".join(map(str,c['adf'])),
					"ADR": ",".join(map(str,c['adr'])),
					"NTN": str(c['count_ltpos']),
					"TN": str(c['uniq_ltpos']),
					"MGL": str(c['min_mg_start']) if math.isfinite(c['min_mg_start']) else ".",
					"MGR": str(c['max_mg_end']) if c['max_mg_end']>0 else "."
				})

		# filtered call
		if (not c['flt']) and x['flt']:
			sample_hit_flt.append({
				"BC": sample_meta[i].name,
				"ADF": ",".join(map(str,c['adf'])),
				"ADR": ",".join(map(str,c['adr'])),
				"NTN": str(c['count_ltpos']),
				"TN": str(c['uniq_ltpos']),
				"MGL": str(c['min_mg_start']) if math.isfinite(c['min_mg_start']) else ".",
				"MGR": str(c['max_mg_end']) if c['max_mg_end']>0 else "."
			})

	# VCF output
	# rec = vcf_out.new_record(
	# 	contig=x['ctg'],
	# 	start=x['pos'] - 1, #pysam expects 0-based
	# 	stop=start + len(x['ref']),
	# 	alleles=(x['ref'], x['alt']),
	# 	id='.',
	# 	qual=None,
	# 	filter='PASS' if not x['flt'] else 'LowQual'
	# )
	# # Set INFO
	# rec.info['BDP'] = ",".join(map(str, bulk_ad))
	# vt = "."
	# nbc = 0

	def write_with(vt, sample_hits):
		rec = vcf_out.new_record(
			contig=x['ctg'],
			start=x['pos'] - 1, #pysam expects 0-based
			stop=x['pos'] - 1 + len(x['ref']),
			alleles=(x['ref'], x['alt']),
			id='.',
			qual=None,
			filter='PASS' if not x['flt'] else 'LowQual'
		)
		# Set INFO
		rec.info['BDP'] = ",".join(map(str, bulk_ad))
		rec.info['VT']  = vt
		rec.info['NBC'] = len(sample_hits)
		rec.samples[output_sample_name]["BC"]  = ";".join(h["BC"]  for h in sample_hits)
		rec.samples[output_sample_name]["ADF"] = ";".join(h["ADF"] for h in sample_hits)
		rec.samples[output_sample_name]["ADR"] = ";".join(h["ADR"] for h in sample_hits)
		rec.samples[output_sample_name]["NTN"] = ";".join(h["NTN"] for h in sample_hits)
		rec.samples[output_sample_name]["TN"]  = ";".join(h["TN"]  for h in sample_hits)
		rec.samples[output_sample_name]["MGL"] = ";".join(h["MGL"] for h in sample_hits)
		rec.samples[output_sample_name]["MGR"] = ";".join(h["MGR"] for h in sample_hits)
		vcf_out.write(rec)

	if sample_hit_ds: write_with("DS", sample_hit_ds)
	if sample_hit_ss: write_with("SS", sample_hit_ss)
	if sample_hit_joint: write_with("JOINT", sample_hit_joint)


	# if sample_hit_ds: 
	# 	vt = "DS"
	# 	nbc = len(sample_hit_ds)
	# 	rec.info['VT'] = vt
	# 	rec.info['NBC'] = nbc
	# 	# Set FORMAT
	# 	rec.samples[output_sample_name]["BC"] = ";".join([hit[0] for hit in sample_hit_ds])
	# 	rec.samples[output_sample_name]["ADF"] = ";".join([hit[1] for hit in sample_hit_ds])
	# 	rec.samples[output_sample_name]["ADR"] = ";".join([hit[2] for hit in sample_hit_ds])
	# 	rec.samples[output_sample_name]["NTN"] = ";".join([hit[3] for hit in sample_hit_ds])
	# 	rec.samples[output_sample_name]["TN"] = ";".join([hit[4] for hit in sample_hit_ds])
	# 	rec.samples[output_sample_name]["MGL"] = ";".join([hit[5] for hit in sample_hit_ds])
	# 	rec.samples[output_sample_name]["MGR"] = ";".join([hit[6] for hit in sample_hit_ds])
	# 	# Write record
	# 	vcf_out.write(rec)
	# if sample_hit_ss: 
	# 	vt = "SS"
	# 	nbc = len(sample_hit_ss)
	# 	rec.info['VT'] = vt
	# 	rec.info['NBC'] = nbc
	# 	# Set FORMAT
	# 	rec.samples[output_sample_name]["BC"] = ";".join([hit[0] for hit in sample_hit_ss])
	# 	rec.samples[output_sample_name]["ADF"] = ";".join([hit[1] for hit in sample_hit_ss])
	# 	rec.samples[output_sample_name]["ADR"] = ";".join([hit[2] for hit in sample_hit_ss])
	# 	rec.samples[output_sample_name]["NTN"] = ";".join([hit[3] for hit in sample_hit_ss])
	# 	rec.samples[output_sample_name]["TN"] = ";".join([hit[4] for hit in sample_hit_ss])
	# 	rec.samples[output_sample_name]["MGL"] = ";".join([hit[5] for hit in sample_hit_ss])
	# 	rec.samples[output_sample_name]["MGR"] = ";".join([hit[6] for hit in sample_hit_ss])
	# 	# Write record
	# 	vcf_out.write(rec)
	# if sample_hit_joint: 
	# 	vt = "JOINT"
	# 	nbc = len(sample_hit_joint)
	# 	rec.info['VT'] = vt
	# 	rec.info['NBC'] = nbc
	# 	# Set FORMAT
	# 	rec.samples[output_sample_name]["BC"] = ";".join([hit[0] for hit in sample_hit_joint])
	# 	rec.samples[output_sample_name]["ADF"] = ";".join([hit[1] for hit in sample_hit_joint])
	# 	rec.samples[output_sample_name]["ADR"] = ";".join([hit[2] for hit in sample_hit_joint])
	# 	rec.samples[output_sample_name]["NTN"] = ";".join([hit[3] for hit in sample_hit_joint])
	# 	rec.samples[output_sample_name]["TN"] = ";".join([hit[4] for hit in sample_hit_joint])
	# 	rec.samples[output_sample_name]["MGL"] = ";".join([hit[5] for hit in sample_hit_joint])
	# 	rec.samples[output_sample_name]["MGR"] = ";".join([hit[6] for hit in sample_hit_joint])
	# 	# Write record
	# 	vcf_out.write(rec)
	# rec.info['VT'] = vt
	# rec.info['NBC'] = nbc
	# # Set FORMAT
	# rec.samples[output_sample_name]["BC"] = ";".join([hit[0] for hit in sample_hit])
	# rec.samples[output_sample_name]["ADF"] = ";".join([hit[1] for hit in sample_hit])
	# rec.samples[output_sample_name]["ADR"] = ";".join([hit[2] for hit in sample_hit])
	# rec.samples[output_sample_name]["NTN"] = ";".join([hit[3] for hit in sample_hit])
	# rec.samples[output_sample_name]["TN"] = ";".join([hit[4] for hit in sample_hit])
	# rec.samples[output_sample_name]["MGL"] = ";".join([hit[5] for hit in sample_hit])
	# rec.samples[output_sample_name]["MGR"] = ";".join([hit[6] for hit in sample_hit])
	# # Write record
	# vcf_out.write(rec)


def main():
	args = parse_args()
	if args.s * 2 > args.a:
		sys.exit("Error: 2 * {-s} should not be larger than {-a}")

	print("VCF input:", args.vcf)
	print("Calling thresholds [a, s]:", args.a, args.s)
	print("Command line:", " ".join(sys.argv))

	sample_excl = read_list(args.fn_excl)
	sample_hap = read_list(args.fn_hap)
	rep_str = load_rep_map(args.fn_rep)
	var_map = load_var_filter(args.fn_var)

	vcf_in, samples, col2sample, sample_meta, rep_id = parse_vcf_header(args, sample_excl, sample_hap, rep_str)
	output_sample_name = os.path.basename(args.vcf).replace(".vcf.gz", "")

	# General options
	n_bulk = args.b
	is_hap_sample = args.H
	min_mapq = args.q
	output_path = args.o

	# Duplex-specific filters
	min_dp_alt = args.a
	min_dp_alt_strand = args.s
	min_dp_dmg_strand = args.S
	min_ab = args.B
	min_joint = args.j
	min_joint_strand = args.J
	max_lt = args.l
	min_end_len = args.L
	flt_win = args.w
	max_dp_ref = args.R
	flt_lv = args.T

	# Bulk-specific filters
	min_dp_bulk = args.D
	min_het_dp_bulk = args.A
	max_alt_dp_bulk = args.m
	is_hap_bulk = args.P
	min_het_ab_bulk = 0.3
	
	last_bulk = []
	last = []
	last_ctg = ''
	germline_indel_adjacent = []

	n_het_bulk = 0
	n_hom_bulk = 0
	n_het_bulk_detected = 0

	# Open output and write header 
	vcf_out = write_output_vcf_header(vcf_in, output_path, output_sample_name)

	#-----------------------------
	# Iterate records in VCF
	#-----------------------------
	for rec in vcf_in:
		ctg = rec.chrom
		pos = rec.pos
		if args.u and not is_autosome(ctg): # auto only
			continue

		# Skip sites with low mapQ
		flt_bulks = False
		flt_ds = False
		if 'AMQ' in rec.info:
			amq = rec.info['AMQ'] # pysam returns tuple
			if amq is not None: 
				valid_vals = [x for x in amq if x is not None and str(x).isdigit()]
				if valid_vals and any(int(x) < min_mapq for x in valid_vals):
					flt_bulks = True
					flt_ds = True

		#--------------------------
		# Parse genotype arrays
		#--------------------------
		BC_samples = []
		bulks = []
		
		for i, sname in enumerate(samples):
			# Map replicate columns
			rep_col = rep_id.get(i, i)
			if i >= n_bulk:
				sample_id = col2sample.get(rep_col)
				if sample_id is None:
					continue

			sample_data = rec.samples[sname]

			# ADF/ADR and DP
			adf, adr, ad, dp_ref, dp_alt = parse_ad(sample_data)
			dp = sum(ad)

			if i < n_bulk:
				bulks.append({'dp':dp, 'ad':ad, 'adf':adf, 'adr':adr})
			else: # for duplex data
				flt = False
				flt_dmg = False

				# LTDROP
				lt = 0
				ltdrop_v = sample_data.get('LTDROP')
				if ltdrop_v is not None:
					try:
						lt = int(ltdrop_v) if isinstance(ltdrop_v, (int,str)) else int(ltdrop_v[0])
					except Exception:
						lt = 0

				# Haploid sample with both alleles
				if sample_meta[sample_id].ploidy == 1 and dp_alt > 0 and dp_ref > 0:
					flt = True
				if lt > max_lt:
					flt = True

				# ALEN filter (min length to end of read)
				alen = sample_data.get('ALEN')
				if alen:
					for v in alen[1:]:
						if v in (None, '.', ''):
							continue
						if float(v) < float(min_end_len):
							flt = True
							flt_dmg = True

				# LTPOS: collect unique left and right positions of read pairs
				ltpos = sample_data.get('LTPOS') or []
				obj_count = {}
				for key in ltpos:
					if key in (None, "", ".", "|"):
						continue
					obj_count[key] = obj_count.get(key, 0) + 1
				uniq_ltpos = ",".join(sorted(obj_count.keys()))
				count_ltpos = len(obj_count)

				# MGPOS: left and right sites of merged window
				mgpos = sample_data.get('MGPOS') or []
				min_mg_start = math.inf
				max_mg_end = 0
				if mgpos:
					for key in mgpos:
						if key in (None, "", ".", "|"):
							continue
						parts = key.split("|")
						if len(parts) == 2:
							mg_start = int(parts[0])
							mg_end = int(parts[1])
							if mg_start < min_mg_start:
								min_mg_start = mg_start
							if mg_end > max_mg_end:
								max_mg_end = mg_end

					# compare indel extents to merged window per -T level
					is_deletion = True
					ref = rec.ref
					alt = rec.alts[0] if rec.alts else ""
					if len(ref)==1 and len(alt)>1:
						is_deletion = False
					indel_start = pos + 1 # indel_start: first inserted/deleted base
					if is_deletion:
						indel_end = indel_start + len(ref) - 2 # indel_end: last deleted base
					else:
						indel_end = indel_start + len(alt) - 2 # indel_end: last inserted base if aligned as soft clip

					if flt_lv == 1: # exact match btw indel and merged window
						if indel_start == min_mg_start and indel_end == max_mg_end:
							flt = True; flt_dmg = True
					elif flt_lv == 2: # allow up to 2bp difference btw indel and merged window
						if (min_mg_start-2) <= indel_start <= (min_mg_start+2) and (max_mg_end-2) <= indel_end <= (max_mg_end+2):
							flt = True; flt_dmg = True
					elif flt_lv == 3: # any overlap btw indel and merged window
						if indel_start <= max_mg_end and (indel_start >= min_mg_start or indel_end >= min_mg_start):
							flt = True; flt_dmg = True

				# Merge multi-appearances of same sample
				if len(BC_samples) <= sample_id or BC_samples[sample_id] is None:
					while len(BC_samples) <= sample_id:
						BC_samples.append(None)
					BC_samples[sample_id] = {
						'flt': flt,
						'dp': dp,
						'ad': ad,
						'adf': adf,
						'adr': adr,
						'lt': lt,
						'flt_dmg': flt_dmg,
						'count_ltpos': count_ltpos,
						'uniq_ltpos': uniq_ltpos,
						'min_mg_start': min_mg_start,
						'max_mg_end': max_mg_end
					}
				else:
					c = BC_samples[sample_id]
					if flt: c['flt'] = True
					if flt_dmg: c['flt_dmg'] = True
					if c['lt'] < lt: c['lt'] = lt
					# use min of ADF/ADR across replicates
					c['dp'] = 0
					for j in range(len(ad)):
						if c['adf'][j] > adf[j]: c['adf'][j] = adf[j]
						if c['adr'][j] > adr[j]: c['adr'][j] = adr[j]
						c['ad'][j] = c['adf'][j] + c['adr'][j]
						c['dp'] += c['ad'][j]
					c['count_ltpos'] = count_ltpos
					c['uniq_ltpos'] = uniq_ltpos
					c['min_mg_start'] = min_mg_start
					c['max_mg_end'] = max_mg_end

		# Only consider biallelic indel sites for calling
		alt_list = list(rec.alts) if rec.alts else []
		if len(alt_list) != 1 or ((len(alt_list[0]) == 1 and len(rec.ref) == 1)):
			flt_bulks = True
			flt_ds = True

		#----------------------------
		# Test het in the bulk(s)
		#----------------------------
		all_het = True
		all_hom = True
		all_good_alt = True
		for b in bulks:
			b['het'] = False
			b['hom'] = False
			# het: both strands support both alleles sufficiently, and AB >= min_het_ab_bulk (0.3)
			if (b['adf'][0] > 0 and b['adf'][1] > 0 and b['adr'][0] > 0 and b['adr'][1] > 0 and
				b['ad'][0] >= min_het_dp_bulk and b['ad'][1] >= min_het_dp_bulk):
				if b['ad'][0] >= b['dp'] * min_het_ab_bulk and b['ad'][1] >= b['dp'] * min_het_ab_bulk:
					b['het'] = True
			if (not b['het']) and b['ad'][1] >= min_het_dp_bulk and b['adf'][1] > 0 and b['adr'][1] > 0 and b['ad'][0] <= max_alt_dp_bulk:
				b['hom'] = True
			if b['ad'][1] < min_het_dp_bulk:
				all_good_alt = False
			if not b['het']:
				all_het = False
			if not b['hom']:
				all_hom = False
			if b['dp'] < min_dp_bulk:
				flt_bulks = True

		# output differences in bulk to error log
		if n_bulk > 1 and not is_hap_bulk:
			n_bulk_ref = 0
			n_bulk_alt = 0
			for b in bulks:
				if b['ad'][1] == 0:
					n_bulk_ref += 1
				elif b['ad'][1] >= min_dp_alt and b['adf'][1] >= min_dp_alt_strand and b['adr'][1] >= min_dp_alt_strand:
					n_bulk_alt += 1
			if n_bulk_ref > 0 and n_bulk_alt > 0:
				while last_bulk and (last_bulk[0]['ctg'] != ctg or last_bulk[0]['pos'] + flt_win < pos):
					x = last_bulk.pop(0)
					if not x['flt']:
						print('BV', *x['data'], sep='\t')
				ad_pairs = [f"{b['adf'][1]}:{b['adr'][1]}" for b in bulks]
				flt_this = flt_bulks
				if var_map and f"{ctg}:{pos}" in var_map:
					flt_this = True
				for j in range(len(last_bulk)):
					flt_this = True
					last_bulk[j]['flt'] = True
				last_bulk.append({
					'flt': flt_this,
					'ctg': ctg,
					'pos': pos,
					'data': [ctg, pos, rec.ref, alt_list[0] if alt_list else '.', *ad_pairs]
				})

		#----------------------
		# Count ADO
		#----------------------
		# BC_sample allele depth < min_joint (2 by default)
		if is_hap_bulk and all_hom and not flt_bulks: # for hap bulk, only count alt allele dropped
			n_hom_bulk += 1
			for j, c in enumerate(BC_samples):
				if c is None: continue
				if c['flt'] or c['ad'][1] < min_joint:
					sample_meta[j].ado[1] += 1
		if (not is_hap_bulk) and all_het and not flt_bulks:
			n_het_bulk += 1
			for j, c in enumerate(BC_samples):
				if c is None: continue
				if c['flt'] or c['ad'][0] < min_joint:
					sample_meta[j].ado[0] += 1 # ref allele dropped
				if c['flt'] or c['ad'][1] < min_joint:
					sample_meta[j].ado[1] += 1 # alt allele dropped

		#------------------------------------------
		# Test if ALT is callable and count FN
		#------------------------------------------
		n_joint_alt = 0
		alt_detected = False
		for i, c in enumerate(BC_samples):
			if c is None:
				continue
			# If sample is haploid and it has ref alleles, c.flt will be true
			c['alt'] = (not c['flt'] and
						c['ad'][1] >= min_dp_alt and
						c['adf'][1] >= min_dp_alt_strand and
						c['adr'][1] >= min_dp_alt_strand and
						c['ad'][1] >= c['dp'] * min_ab and
						c['ad'][0] <= max_dp_ref)
			c['joint_alt'] = (not c['flt'] and
							  c['ad'][1] >= min_joint and
							  c['adf'][1] >= min_joint_strand and
							  c['adr'][1] >= min_joint_strand)
			if c['joint_alt']:
				n_joint_alt += 1
			# Count FN when bulk is ALT but sample is not
			if (not flt_bulks) and (not c['alt']) and ((is_hap_bulk and all_hom) or ((not is_hap_bulk) and all_het)):
				sample_meta[i].fn += 1
			if (not flt_bulks) and c['alt']:
				alt_detected = True
			# Call single-stranded variants (DNA damage)
			c['dmg'] = (not c['flt_dmg'] and
						c['ad'][1] >= min_dp_dmg_strand and c['ad'][0] >= min_dp_dmg_strand and
						(c['adf'][1] * c['adr'][1] == 0) and
						(c['adf'][0] * c['adf'][1] == 0) and
						(c['adr'][0] * c['adr'][1] == 0))
			if all_het and not flt_bulks:
				if not ((not c['flt_dmg']) and c['adf'][0] >= min_dp_dmg_strand and c['adr'][0] >= min_dp_dmg_strand):
					sample_meta[i].dmg_fn[0] += 1 # bulk het, but sample doesn't have enough REF support, count dmg FN for REF
				if not ((not c['flt_dmg']) and c['adf'][1] >= min_dp_dmg_strand and c['adr'][1] >= min_dp_dmg_strand):
					sample_meta[i].dmg_fn[1] += 1 # bulk het, but sample doesn't have enough ALT support, count dmg FN for ALT
				if sample_meta[i].ploidy > 1 and c['dmg']:
					sample_meta[i].dmg_fp += 1 # bulk het, but sample dmg (ALT only on one strand), count dmg FP (no dmg_fp of this kind for haploid sample)

		# Count detected germline het in bulk
		if (not flt_bulks) and alt_detected and ((is_hap_bulk and all_hom) or ((not is_hap_bulk) and all_het)):
			n_het_bulk_detected += 1

		# Skip germline variants for calling, allowing somatic indel next to germline SNP; if germline indel (het or ALT-hom), add to list for filtering somatic calls
		if all_good_alt:
			# if germline indel, add adjacent region to list for filtering
			is_indel = (len(rec.ref) != 1 or len(alt_list[0]) != 1) if alt_list else False
			if (not flt_bulks) and (all_hom or all_het) and is_indel:
				# adjacent region: flanking bp = max(5, twice indel length)
				flank_default = 5
				indel_length = abs(len(rec.ref) - len(alt_list[0]))
				span = max(flank_default, 2*indel_length)
				germline_indel_adjacent.append((pos - span, pos + span))
			# skip site for somatic calling
			continue

		# Require at least one bulk with good RefHom
		n_bulk_ref = sum(1 for b in bulks if b['ad'][1] <= max_alt_dp_bulk)
		if n_bulk_ref == 0:
			flt_ds = True # flag the infavorable scenario: no bulks with good RefHom; this site may be used for window filtering later

		#--------------------------------------
		# Apply filters when contig changes
		#--------------------------------------
		if last_ctg != ctg and last:
			# apply adjacent germline indel filter
			for k in range(len(last)):
				if filter_germline_adj(last[k], germline_indel_adjacent):
					last[k]['flt'] = True
			# apply window filter: start from the first variant in the contig, if it's not filtered and its distance to the next variant is larger than flt_win, print variant
			while len(last) > 1:
				x = last.pop(0)
				if x['pos'] + flt_win > last[0]['pos']:
					x['flt'] = True
					last[0]['flt'] = True
				if args.F or (not x['flt']):
					print_calls(vcf_out, output_sample_name, x, sample_meta, args.F, min_joint, min_joint_strand)
			# print last variant
			x = last.pop(0)
			if args.F or (not x['flt']):
				print_calls(vcf_out, output_sample_name, x, sample_meta, args.F, min_joint, min_joint_strand)
			# reset germline list for the new contig
			germline_indel_adjacent.clear()

		flt_this = flt_ds or flt_bulks
		if var_map and f"{ctg}:{pos}" in var_map:
			flt_this = True

		last.append({
			'flt': flt_this,
			'n_joint_alt': n_joint_alt,
			'ctg': ctg,
			'pos': pos,
			'bulks': bulks,
			'BC_samples': [c if c is not None else {
				'flt': True, 'dp':0, 'ad':[0,0], 'adf':[0,0], 'adr':[0,0],
				'lt':0, 'flt_dmg':True, 'count_ltpos':0, 'uniq_ltpos':'',
				'min_mg_start': math.inf, 'max_mg_end':0,
				'alt': False, 'joint_alt': False, 'dmg': False
			} for c in BC_samples],
			'ref': rec.ref,
			'alt': alt_list[0] if alt_list else '.'
		})
		last_ctg = ctg

	# print BV to log
	while last_bulk:
		x = last_bulk.pop(0)
		if not x['flt']:
			print('BV', *x['data'], sep='\t')

	# Print last contig calls
	if last:
		for k in range(len(last)):
			if filter_germline_adj(last[k], germline_indel_adjacent):
				last[k]['flt'] = True
		while len(last) > 1:
			x = last.pop(0)
			if x['pos'] + flt_win > last[0]['pos']:
				x['flt'] = True
				last[0]['flt'] = True
			if args.F or (not x['flt']):
				print_calls(vcf_out, output_sample_name, x, sample_meta, args.F, min_joint, min_joint_strand)
		x = last.pop(0)
		if args.F or (not x['flt']):
			print_calls(vcf_out, output_sample_name, x, sample_meta, args.F, min_joint, min_joint_strand)
		germline_indel_adjacent.clear()

	#-------------------------
	# Output final stats
	#-------------------------
	mut = []
	fnr = []
	corr_mut = []
	dmg = []
	corr_dmg = []
	ado = []
	fnr_dmg = []
	fpr_dmg = []

	# avoid division by zero
	n_het_bulk_safe = max(n_het_bulk, 1)
	n_hom_bulk_safe = max(n_hom_bulk, 1)

	for i, c in enumerate(sample_meta):
		if is_hap_bulk:
			ado_i = c.ado[1] / n_hom_bulk_safe
			fnr_i = c.fn / n_hom_bulk_safe
		else:
			if c.ploidy == 1:
				ado_i = 2.0 * c.ado[1] / n_het_bulk_safe - 1.0
				fnr_i = 2.0 * c.fn / n_het_bulk_safe - 1.0
			else:
				ado_i = c.ado[1] / n_het_bulk_safe
				fnr_i = c.fn / n_het_bulk_safe
		mut.append(c.mut)
		fnr.append(f"{fnr_i:.4f}")
		# safe division
		denom = 1.0 - float(fnr_i)
		if denom <= 0:
			denom = 1e-12
		corr_mut.append(f"{(c.mut / denom):.2f}")

		dmg.append(c.dmg)
		if not is_hap_bulk:
			fpr_dmg_i = c.dmg_fp / n_het_bulk_safe
			fnr_dmg_i = (2.0 * c.dmg_fn[1] / n_het_bulk_safe - 1.0) if c.ploidy == 1 else (c.dmg_fn[1] / n_het_bulk_safe)
			# safe division
			denom2 = 1.0 - float(fnr_i)
			if denom2 <= 0:
				denom2 = 1e-12
			corr_dmg_i = (c.dmg - float(corr_mut[-1]) * fpr_dmg_i) / denom2
			fpr_dmg.append(f"{fpr_dmg_i:.4f}")
			fnr_dmg.append(f"{fnr_dmg_i:.4f}")
			corr_dmg.append(f"{corr_dmg_i:.2f}")

	print('NN', *[str(x) for x in mut], sep='\t')
	print('NR', *fnr, sep='\t')
	print('NC', *corr_mut, sep='\t')
	print('DN', *[str(x) for x in dmg], sep='\t')
	if not is_hap_bulk:
		print('DP', *fpr_dmg, sep='\t')
		print('DR', *fnr_dmg, sep='\t')
		print('DC', *corr_dmg, sep='\t')

	# Sensitivity
	print('Sensitivity_binary', n_het_bulk_detected, n_het_bulk, f"{(n_het_bulk_detected/max(1,n_het_bulk)):.4f}", sep='\t')
	sum_fnr = sum(float(x) for x in fnr) if fnr else 0.0
	print('Sensitivity_FNR', len(fnr), f"{sum_fnr:.4f}", f"{(len(fnr) - sum_fnr):.4f}", sep='\t')

	# Close files
	vcf_out.close()
	vcf_in.close()


if __name__ == '__main__':
	main()




