
import argparse
import HTSeq
import os
import pyfaidx
import regex
from statsmodels.distributions.empirical_distribution import ECDF
import sys
import numpy as np
from utility import reverseComplement
import pysam
import pandas as pd
import traceback
import logging
logger = logging.getLogger('root')
logger.propagate = False

# reference_start is 0-index, reference_end is 1-index
def is_outie(x,y):
	if x.is_reverse:
		if not y.is_reverse:
			if x.reference_start < y.reference_start:
				return True,x.reference_end - y.reference_start
	else:
		if y.is_reverse:
			if y.reference_start < x.reference_start:
				return True,y.reference_end - x.reference_start
	return False,0


# make sure BWA version >= 0.7.7
def get_conversion_pos(r,overlap_length=7):
	if overlap_length<0:
		return []
	out = []
	actual_start = 0
	# expected pos is len(overlap) + 2, 1-index, we gave it 2bp flank
	expected_min = overlap_length
	expected_max = overlap_length+4
	try:
		myList = r.get_aligned_pairs(with_seq=True)
	except Exception as e:
		logger.error(traceback.format_exc())
		return []
	# qualities = r.query_qualities
	if r.is_reverse:
		myList = myList[::-1]
	q_seq = r.query_sequence
	for q,ref_pos,base in myList:
		if q == None:
			continue
		actual_start+=1
		# if qualities[q]<=10:
		# 	continue
		if ref_pos == None:
			continue
		if expected_min<=actual_start <=expected_max:
			if q_seq[q] == "G" and base.upper() == "A":
				out.append(ref_pos)
			if q_seq[q] == "C" and base.upper() == "T":
				out.append(ref_pos)
		if actual_start>expected_max:
			break
	return out

def get_read_strand(read):
	if read.is_reverse:
		return "-"
	return "+"

def get_read_start(read):
	if read.is_reverse:
		return read.reference_end
	return read.reference_start+1

def tabulate_start_positions_BE(bam=None,label=None,output_dir=None, # inputs
								min_overlap=5,max_overlap=15,mapq_threshold=0, **kwargs): # other filters

	# print (bam)
	# output files
	SE_read_stats_file = f"{output_dir}/{label}.SE.read_stats.csv"
	SE_read_stats_list_of_dict = [] #  # read, chr, pos, converted_pos_list
	PE_read_stats_file = f"{output_dir}/{label}.PE.read_stats.csv"
	PE_read_stats_list_of_dict = [] #  # read, chr, pos1, pos2, is_outie, overlap_bp, converted_pos_list

	# expected conversion position is overlap_bp+2

	# housekeeping
	logger.info("Reading bam to a dictionary")
	reads_dict = bam_to_dict(pysam.AlignmentFile(bam,"rb"),mapq_threshold)
	logger.info("Finished bam to dictionary")
	ga_coverage_paired = HTSeq.GenomicArray("auto", stranded=False,typecode="d") # same as old ga_coverage
	ga_coverage_single = HTSeq.GenomicArray("auto", stranded=False,typecode="d") # subset of ga_noise, save this var in case we want to rescue some reads
	ga_converted = HTSeq.GenomicArray("auto", stranded=False,typecode="d")
	ga_overlap = HTSeq.GenomicArray("auto", stranded=False,typecode="d")
	ga_noise = HTSeq.GenomicArray("auto", stranded=False,typecode="d")


	# main
	read_count = 0
	SE_read_count = 0
	PE_read_count = 0
	for read in reads_dict:
		read_count += 1
		if not read_count % 100000:
			print(read_count/float(1000000), end=" ", file=sys.stderr)
		read1_list = reads_dict[read][0]
		read2_list = reads_dict[read][1]

		# parse single-end reads
		if len(read1_list) == 0 or len(read2_list) == 0:
			if len(read1_list) == 0:
				ending = "/1"
			if len(read2_list) == 0:
				ending = "/2"
			SE_read_count+=1
			read_list = read1_list+read2_list
			for r in read_list:
				read_start = get_read_start(r)
				# ga_coverage_single[HTSeq.GenomicPosition(r.reference_name, read_start)]+=1
				# converted_pos_list = get_conversion_pos(r)
				converted_pos_list = []
				SE_read_stats_list_of_dict.append({"read":read+ending,
												   "chr":r.reference_name,
												   "pos":read_start,
												   "strand":get_read_strand(r),
												   "converted_pos_list":converted_pos_list})
			# print (f"{read} (SE) is noise.")
			# ga_noise[HTSeq.GenomicPosition(r.reference_name, read_start)]+=1
			ga_noise[HTSeq.GenomicInterval(r.reference_name, r.reference_start,r.reference_end)]+=1

			continue
		PE_read_count += 1
		noise_flag = False # only count noise once
		# enumerate all possible pairs
		for i in read1_list:
			read1_start = get_read_start(i)
			read1_strand = get_read_strand(i)
			for j in read2_list:
				read2_start = get_read_start(j)
				read2_strand = get_read_strand(j)
				if i.reference_name != j.reference_name:
					noise_flag=True
					SE_read_stats_list_of_dict.append({"read":read+"/1",
													   "chr":i.reference_name,
													   "pos":read1_start,
													   "strand":read1_strand,
													   "converted_pos_list":[]})
					SE_read_stats_list_of_dict.append({"read":read+"/2",
													   "chr":j.reference_name,
													   "pos":read2_start,
													   "strand":read2_strand,
													   "converted_pos_list":[]})
					continue
				flag,overlap_bp = is_outie(i,j)
				if flag:
					R1_converted_list = get_conversion_pos(i,overlap_bp)
					R2_converted_list = get_conversion_pos(j,overlap_bp)
					converted_pos_list = R1_converted_list+R2_converted_list
					if min_overlap<=overlap_bp<=max_overlap:
						ga_coverage_paired[HTSeq.GenomicPosition(i.reference_name, read1_start)]+=1
						ga_coverage_paired[HTSeq.GenomicPosition(i.reference_name, read2_start)]+=1
						for pos in converted_pos_list:
							ga_converted[HTSeq.GenomicPosition(i.reference_name, pos)]+=1
					else:
						noise_flag=True
				else:
					converted_pos_list = []
					noise_flag=True
				# read, chr, pos1, pos2, is_outie, overlap_bp, converted_pos_list
				PE_read_stats_list_of_dict.append({"read":read,
												   "chr":i.reference_name,
												   "pos1":i.reference_start,
												   "strand1":read1_strand,
												   "pos2":j.reference_start,
												   "strand2":read2_strand,
												   "is_outie":flag,
												   "overlap_bp":overlap_bp,
												   "converted_pos_list":converted_pos_list})
		if noise_flag:
			# ga_noise[HTSeq.GenomicPosition(i.reference_name, read1_start)]+=1
			# ga_noise[HTSeq.GenomicPosition(i.reference_name, read2_start)]+=1
			ga_noise[HTSeq.GenomicInterval(i.reference_name, i.reference_start,i.reference_end)]+=1
			ga_noise[HTSeq.GenomicInterval(j.reference_name, j.reference_start,j.reference_end)]+=1
			# print (f"{read} (PE) is noise.")
	pd.DataFrame.from_dict(SE_read_stats_list_of_dict).to_csv(SE_read_stats_file,index=False)
	pd.DataFrame.from_dict(PE_read_stats_list_of_dict).to_csv(PE_read_stats_file,index=False)


	return ga_coverage_paired, ga_coverage_single, ga_converted, ga_overlap, ga_noise, PE_read_count


""" Find genomic windows (coordinate positions)
"""
def find_windows(ga_windows, window_size):
	# Initialize comparison position
	last = HTSeq.GenomicInterval("0", 0, 0)
	# Iterate through window GenomicArray and consolidate windows that are within 3 bp, up to a maximum of 10 bp.
	for iv, value in ga_windows.steps():
		if value:
			if iv.chrom != last.chrom or iv.start - last.end > window_size or iv.end - last.start > 10:
				last = iv
			else:
				consolidated_interval = HTSeq.GenomicInterval(iv.chrom, last.start, iv.end)
				ga_windows[consolidated_interval] = 1
				last = consolidated_interval

	return ga_windows  # Return consolidated GenomicArray

""" Find actual sequences of potential off-target sites
"""
def output_alignments(narrow_ga, narrow_ga_converted,narrow_ga_noise, ga_windows, reference_genome, target_sequence, target_name, target_cells,
					  bam_filename, mismatch_threshold, ga_pval, search_radius, out,nuclease_ga_overlap=None,control_ga_overlap=None):

	# dictionary to store the matched reads
	matched_dict = {}   
	# dictionary to add read_count for each pair chromosome:start_position among matched reads
	reads_dict = {}
	# dictionary to store window_start. For duplicated matched off-target.
	window_min = {}
	# dictionary to store window_end. For duplicated matched off-target.
	window_max = {}

	# dictionary to store the unmatched reads
	unmatched_dict = {}
	
	for iv, value in ga_windows.steps():
		if value:
			window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - search_radius, iv.end + search_radius)

			offtarget_sequence_no_bulge, mismatches, offtarget_sequence_length, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
			realigned_target, \
			bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end = \
				alignSequences(target_sequence, window_sequence, max_score=mismatch_threshold)

			# get genomic coordinates of sequences
			mm_start, mm_end, b_start, b_end = '', '', '', ''
			if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '+':
				mm_start = iv.start - search_radius + int(start_no_bulge)
				mm_end = iv.start - search_radius + int(end_no_bulge)
			if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '-':
				mm_start = iv.end + search_radius - int(end_no_bulge)
				mm_end = iv.end + search_radius - int(start_no_bulge)

			if bulged_offtarget_sequence and chosen_alignment_strand_b == '+':
				b_start = iv.start - search_radius + int(bulged_start)
				b_end = iv.start - search_radius + int(bulged_end)
			if bulged_offtarget_sequence and chosen_alignment_strand_b == '-':
				b_start = iv.end + search_radius - int(bulged_end)
				b_end = iv.end + search_radius - int(bulged_start)

			#  define overall start, end and strand. For bed annotation purposes
			if offtarget_sequence_no_bulge:
				target_start_absolute = mm_start
				target_end_absolute = mm_end
				target_strand_absolute = chosen_alignment_strand_m
			elif not offtarget_sequence_no_bulge and bulged_offtarget_sequence:
				target_start_absolute = b_start
				target_end_absolute = b_end
				target_strand_absolute = chosen_alignment_strand_b
			else:
				target_start_absolute = iv.start
				target_end_absolute = iv.end
				target_strand_absolute = '*'

			name = iv.chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
			read_count = int(max(set(narrow_ga[iv])))
			filename = os.path.basename(bam_filename)
			full_name = str(target_name) + '_' + str(target_cells) + '_' + str(name) + '_' + str(read_count)
			"""
			control_ga_overlap_bp_list = []
			nuclease_ga_overlap_bp_list = []
			if nuclease_ga_overlap != None:
				for item in nuclease_ga_overlap[ HTSeq.GenomicInterval( iv.chrom, target_start_absolute,target_end_absolute, "+" ) ].steps():
					if item[1] == None:
						continue
					if len(item[1]) > len(nuclease_ga_overlap_bp_list):
						nuclease_ga_overlap_bp_list = item[1]

				for item in control_ga_overlap[ HTSeq.GenomicInterval( iv.chrom, target_start_absolute,target_end_absolute, "+" ) ].steps():
					if item[1] == None:
						continue
					if len(item[1]) > len(control_ga_overlap_bp_list):
						control_ga_overlap_bp_list = item[1]
			nuclease_ga_overlap_bp_list = ",".join([str(x) for x in nuclease_ga_overlap_bp_list])
			control_ga_overlap_bp_list = ",".join([str(x) for x in control_ga_overlap_bp_list])
			"""
			nuclease_ga_overlap_bp_list=str(int(max(set(narrow_ga_converted[iv]))))+","+str(int(max(set(narrow_ga_noise[iv])))) # tmp
			control_ga_overlap_bp_list= str(int(max(set(narrow_ga_converted[iv]))))+","+str(int(max(set(narrow_ga_noise[iv])))) # tmp
			# control_ga_overlap_bp_list= str(int(max(set(control_ga_converted[iv]))))+","+str(int(max(set(control_ga_noise[iv])))) # tmp
			# print (iv.chrom,iv.start,read_count,nuclease_ga_overlap_bp_list)
			if offtarget_sequence_no_bulge or bulged_offtarget_sequence:
				tag = iv.chrom + ':' + str(target_start_absolute)
				if tag not in reads_dict.keys():
					reads_dict[tag] = read_count
					window_min[tag] = [iv.start]
					window_max[tag] = [iv.end]
					matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name, read_count, target_strand_absolute,
										 iv.start, iv.end, iv, window_sequence,
										 offtarget_sequence_no_bulge, mismatches,
										 chosen_alignment_strand_m, mm_start, mm_end,
										 bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
										 chosen_alignment_strand_b, b_start, b_end,
										 filename, target_cells, target_name, full_name, target_sequence, realigned_target,nuclease_ga_overlap_bp_list,control_ga_overlap_bp_list]
				else:
					current_read_count = reads_dict[tag]
					reads_dict[tag] = max(current_read_count, read_count) 
					window_min[tag].append(iv.start)
					window_max[tag].append(iv.end)
					matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name, reads_dict[tag], target_strand_absolute,
										 min(window_min[tag]), max(window_max[tag]), iv, window_sequence,
										 offtarget_sequence_no_bulge, mismatches,
										 chosen_alignment_strand_m, mm_start, mm_end,
										 bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
										 chosen_alignment_strand_b, b_start, b_end,
										 filename, target_cells, target_name, full_name, target_sequence, realigned_target,nuclease_ga_overlap_bp_list,control_ga_overlap_bp_list]
			else:
				untag = iv.chrom + ':' + str(iv.start)
				unmatched_dict[untag] = [iv.chrom, target_start_absolute, target_end_absolute, name, read_count, target_strand_absolute,
										 iv.start, iv.end, iv, window_sequence,
										 offtarget_sequence_no_bulge, mismatches,
										 chosen_alignment_strand_m, mm_start, mm_end,
										 bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
										 chosen_alignment_strand_b, b_start, b_end,
										 filename, target_cells, target_name, full_name, target_sequence, 'none',nuclease_ga_overlap_bp_list,control_ga_overlap_bp_list]

	# Write matched table
	
	# Yichao, add control reads
	# print(f"Writing matched table", file=sys.stderr)
	tags_sorted = matched_dict.keys()
	tags_sorted = sorted(tags_sorted)
	outfile_matched = '{0}_identified_matched.txt'.format(out)

	o1 = open(outfile_matched, 'w')
	# print('Chromosome', 'Start', 'End', 'Name', 'ReadCount', 'Strand',  # 0:5
		  # 'MappingPositionStart', 'MappingPositionEnd', 'WindowName', 'WindowSequence',  # 6:9
		  # 'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',  # 10:11
		  # 'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',  # 12:14
		  # 'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',  # 15:17
		  # 'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',  # 18:20
		  # 'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',  #21:23
		  # 'FileName', 'Cell', 'Targetsite', 'FullName', 'TargetSequence', 'RealignedTargetSequence',  # 24:29
		  # 'Position.Pvalue', 'Narrow.Pvalue', 'Position.Control.Pvalue', 'Narrow.Control.Pvalue','control_position_counts','control_window_counts',  # 30:33
		  # sep='\t', file=o1)
	# Yichao Redefine output
	print('#Chromosome', 'Start', 'End', 'Genomic Coordinate', 'Nuclease_Read_Count', 'Strand',  # 0:5 bed6 format
		  'Control_Read_Count','Site_Sequence','Site_Substitution_Number','Site_Sequence_Gaps_Allowed', # contron window count, # 10:11, 15
		  'File_Name', 'Cell', 'Target_site', 'Full_Name', 'Target_Sequence', 'Realigned_Target_Sequence',  # 24:29
		  'Nuclease_overlap_bp_list', 'Control_overlap_bp_list',  # which column, -2 -1
		  sep='\t', file=o1)
	o1.close()

	with open(outfile_matched, 'a') as o1:
		for key in tags_sorted:
			row = matched_dict[key]
			
			control_position_counts, control_window_counts = list(), list()
			
			iv_pval = HTSeq.GenomicInterval(row[0], int(row[1]), int(row[2]), '.')
			for interval, value in ga_pval[iv_pval].steps():
				 if value is not None:
					 control_position_counts.append(value[3])
					 control_window_counts.append(value[5])
			

			control_position_counts = np.mean(control_position_counts)
			control_window_counts = np.mean(control_window_counts)

			outline = [row[row_index] for row_index in [0,1,2,3,4,5]]
			outline += [control_window_counts]
			# outline += [row[row_index] for row_index in [10,11,15,24,25,26,27,28,29]]
			outline += [row[row_index] for row_index in [10,11,15,24,25,26,27,28,29,-2,-1]]
			print(*(outline), sep='\t', file=o1)


	# Write unmatched table
	# print("Writing unmatched table", file=sys.stderr)
	untags_sorted = unmatched_dict.keys()
	untags_sorted = sorted(untags_sorted)
	outfile_unmatched = '{0}_identified_unmatched.txt'.format(out)
	with open(outfile_unmatched, 'w') as o2:
		for unkey in untags_sorted:
			unrow = unmatched_dict[unkey]

			# direct code copy from matched dict
			control_position_counts, control_window_counts = list(), list()
			
			iv_pval = HTSeq.GenomicInterval(unrow[0], int(unrow[1]), int(unrow[2]), '.')
			for interval, value in ga_pval[iv_pval].steps():
				 if value is not None:
					 control_position_counts.append(value[3])
					 control_window_counts.append(value[5])
			

			control_position_counts = np.mean(control_position_counts)
			control_window_counts = np.mean(control_window_counts)
			unrow += [control_window_counts]
			# un_pos_pval_list, un_nar_pval_list = list(), list()
			# un_control_pos_pval_list, un_control_nar_pval_list = list(), list()
			#
			# iv_pval = HTSeq.GenomicInterval(unrow[0], int(unrow[1]), int(unrow[2]), '.')
			# for interval, value in ga_pval[iv_pval].steps():
			#	 if value is not None:
			#		 un_pos_pval_list.append(value[0])
			#		 un_nar_pval_list.append(value[1])
			#		 un_control_pos_pval_list.append(value[2])
			#		 un_control_nar_pval_list.append(value[3])
			#
			# un_pval_pos = min(un_pos_pval_list)
			# un_pval_nar = min(un_nar_pval_list)
			# un_control_pval_pos = min(un_control_pos_pval_list)
			# un_control_pval_nar = min(un_control_nar_pval_list)

			print(*(unrow), sep='\t', file=o2)


""" Reverse complement DNA sequence
"""
def reverseComplement(seq):
	compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'X': 'X',
				  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', 'x': 'x',
				  'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K',
				  'y': 'r', 'r': 'y', 'w': 'w', 's': 's', 'k': 'm', 'm': 'k',
				  'D': 'H', 'H': 'D', 'V': 'B', 'B': 'V',
				  'd': 'h', 'h': 'd', 'v': 'b', 'b': 'v',
				  '.': '.', '-': '-', '_': '_'})
	out_list = [compl[bp] for bp in seq]
	return ''.join(out_list[::-1])


def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
	seq = seq.upper()
	"""
	Given a sequence with ambiguous base characters, returns a regex that matches for
	the explicit (unambiguous) base characters
	"""
	IUPAC_notation_regex = {'N': '[ATCGN]',
							'Y': '[CTY]',
							'R': '[AGR]',
							'W': '[ATW]',
							'S': '[CGS]',
							'H': '[ACTH]',
							'M': '[ACM]',
							'B': '[CGTB]',
							'K': '[GTK]',
							'D': '[AGTD]',
							'V': '[ACGV]',
							'A': 'A',
							'T': 'T',
							'C': 'C',
							'G': 'G'}

	pattern = ''

	for c in seq:
		pattern += IUPAC_notation_regex[c]

	if lookahead:
		pattern = '(?b:' + pattern + ')'

	pattern_standard = pattern + '{{s<={0}}}'.format(errors)
	pattern_gap = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
	return pattern_standard, pattern_gap

"""
Allow for '-' in our search, but do not allow insertions or deletions. 
"""
def extendedPattern(seq, errors):
	IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
	realign_pattern = ''
	for c in seq:
		realign_pattern += IUPAC_notation_regex_extended[c]
	return '(?b:' + realign_pattern + ')' + '{{s<={0}}}'.format(errors)

"""
Recreate A!!! sequence in the window_sequence that matches the conditions given for the fuzzy regex. 
Currently only working for off-targets with at most one bulge !!! 
"""
## Yichao: Fix the code to get the 'realigned target sequence' to use the new feature in regex library that provides fuzzy_changes
def realignedSequences(targetsite_sequence, chosen_alignment, errors):

	m = chosen_alignment.group()
	q = targetsite_sequence
	substitutions, insertions, deletions = chosen_alignment.fuzzy_changes
	
	start = chosen_alignment.span()[0]	
	indels = {}
	# deletion index is for targetsite_sequence
	# insertion index is for match sequence
	insertions = [i-start for i in insertions]
	deletions = [i-start for i in deletions]
	q_list = []
	m_list = []
	for i in insertions:
		indels[i] = True
	for d in deletions:
		indels[d] = False
	
	count = 0
	q_index = 0
	m_index = 0
	count = 0
	while q_index < len(q) or m_index < len(m):
		if q_index in indels:	
			if not indels[q_index]:
				m_list.append("-")
				q_list.append(q[q_index])
				q_index+=1	
				count+=1
				continue
		if m_index in indels:
			if indels[m_index]:
				q_list.append("-")
				m_list.append(m[m_index])
				
				## update indel dict, I found the deletion positions are relative if insertion occurs first
				update_indels=[]
				for k in indels.keys():
					if k > m_index:
						if not indels[k]:
							del indels[k]
							update_indels.append(k-1)
				for u in update_indels:
					indels[u] = False
				m_index+=1	
				count+=1
				continue
		m_list.append(m[m_index])
		q_list.append(q[q_index])		
		q_index+=1
		m_index+=1
		count+=1
	
	realigned_target_sequence = "".join(q_list)
	realigned_offtarget_sequence =  "".join(m_list)
	
	
	return realigned_target_sequence, realigned_offtarget_sequence


"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match(es).
"""
def alignSequences(targetsite_sequence, window_sequence, max_score=7):

	window_sequence = window_sequence.upper()
	query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)

	# Try both strands
	alignments_mm, alignments_bulge = list(), list()
	alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
	alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
	alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
	alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

	lowest_distance_score, lowest_mismatch = 100, max_score + 1
	chosen_alignment_b, chosen_alignment_m, chosen_alignment_strand_b, chosen_alignment_strand_m = None, None, '', ''

	# Use regex to find the best match allowing only for mismatches
	for aln_m in alignments_mm:
		strand_m, alignment_type_m, match_m = aln_m
		if match_m != None:
			mismatches, insertions, deletions = match_m.fuzzy_counts
			if mismatches < lowest_mismatch:
				chosen_alignment_m = match_m
				chosen_alignment_strand_m = strand_m
				lowest_mismatch = mismatches

	# Use regex to find the best match allowing for gaps, so that its edit distance is strictly lower than the
	# total number of mismatches of the sequence founded (if any) allowing only for mismatches.
	for aln_b in alignments_bulge:
		strand_b, alignment_type_b, match_b = aln_b
		if match_b != None:
			substitutions, insertions, deletions = match_b.fuzzy_counts
			if insertions or deletions:
				distance_score = substitutions + (insertions + deletions) * 3
				edistance = substitutions + insertions + deletions
				if distance_score < lowest_distance_score and edistance < lowest_mismatch:
					chosen_alignment_b = match_b
					chosen_alignment_strand_b = strand_b
					lowest_distance_score = distance_score

	if chosen_alignment_m:
		offtarget_sequence_no_bulge = chosen_alignment_m.group()
		mismatches = chosen_alignment_m.fuzzy_counts[0]
		start_no_bulge = chosen_alignment_m.start()
		end_no_bulge = chosen_alignment_m.end()
	else:
		offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''

	bulged_offtarget_sequence, score, length, substitutions, insertions, deletions, bulged_start, bulged_end, realigned_target = \
		'', '', '', '', '', '', '', '', 'none'
	if chosen_alignment_b:
		realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_b, max_score)
		if bulged_offtarget_sequence:
			length = len(chosen_alignment_b.group())
			substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
			score = substitutions + (insertions + deletions) * 3
			bulged_start = chosen_alignment_b.start()
			bulged_end = chosen_alignment_b.end()
		else:
			chosen_alignment_strand_b = ''

	return [offtarget_sequence_no_bulge, mismatches, len(offtarget_sequence_no_bulge), chosen_alignment_strand_m, start_no_bulge, end_no_bulge,
			realigned_target,
			bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end]


""" Get sequences from some reference genome
"""
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
	if strand == "+":
		seq = reference_genome[chromosome][int(start):int(end)]
	elif strand == "-":
		seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
	return str(seq)

def bam_to_dict(bam,MAPQ=0):
	"""Store all reads, paired/single, into dict
	"""
	output = {}
	for read in bam.fetch(region=None):
		qname = read.query_name
		if read.is_paired:
			if read.is_unmapped:
				continue
			if read.mapping_quality<MAPQ:
				continue
			if not qname in output:
				output[qname]={0:[], 1:[]}
			if read.is_read1:
				output[qname][0].append(read)
			else:
				output[qname][1].append(read)
	return output



def compare(reference_genome=None, bam=None,label=None, control=None, targetsite=None, search_radius=30, window_size=30, mismatch_threshold=6,
			output_dir=None,read_count_cutoff=6,edited_read_cutoff=1,**kwargs):

	# housekeeping variables
	reference_genome_pyfaidx = pyfaidx.Fasta(reference_genome)
	output_list = list()
	bg_position = list()  # List to store nuclease_position_counts that were observed at least once
	bg_narrow = list()  # List to store the sum of nuclease_position_counts in the narrow window

	combined_ga = HTSeq.GenomicArray("auto", stranded=False)  # Store the union of control and nuclease positions
	offtarget_ga_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites
	ga_narrow_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites narrow windows read counts
	ga_narrow_windows_converted = HTSeq.GenomicArray("auto", stranded=False)
	ga_narrow_windows_noise = HTSeq.GenomicArray("auto", stranded=False)
	# output files
	output_filename = label + '_count.txt' # counts per position
	output_count_list_of_dict = []



	# print (label,output_dir)
	nuclease_ga, nuclease_ga_coverage_single, nuclease_ga_converted, nuclease_ga_overlap, nuclease_ga_noise,total_nuclease_count = \
		tabulate_start_positions_BE(bam=bam,label=label,output_dir=output_dir, **kwargs)

	control_ga, control_ga_coverage_single, control_ga_converted, control_ga_overlap, control_ga_noise,total_control_count = \
		tabulate_start_positions_BE(bam=control,label="Control_"+label,output_dir=output_dir, **kwargs)
	logger.info("Finished tabulating read start")
	# For all positions with detected read mapping positions, put into a combined genomicArray
	for iv, value in nuclease_ga.steps():
		if value:
			combined_ga[iv] = 1
	for iv, value in control_ga.steps():
		if value:
			combined_ga[iv] = 1
	print (combined_ga)
	logger.info("Finished combined_ga")
	for iv, value in combined_ga.steps():
		if value:
			for position in iv.range(step=1):
				# Define the windows
				window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - window_size),
											   position.pos + window_size + 1)

				# Start mapping positions, at the specific base position
				nuclease_position_counts = nuclease_ga[position]

				control_position_counts = control_ga[position]
				# Store control_position_counts for which it was observed at least one read
				if control_position_counts >= 0:
					bg_position.append(control_position_counts)

				# In the narrow (parameter-specified) window
				nuclease_window_counts = sum(nuclease_ga[window])
				control_window_counts = sum(control_ga[window])
				# new vars
				nuclease_window_converted_counts = sum(nuclease_ga_converted[window])
				nuclease_window_noise_counts = max(nuclease_ga_noise[window])
				# print (position.chrom,position.pos,nuclease_position_counts,nuclease_window_converted_counts,nuclease_window_noise_counts)
				# Store control_window_counts greater than zero
				if control_window_counts >= 0:
					bg_narrow.append(control_window_counts)

				# A list of the outputs
				row = [position.chrom, position.pos, nuclease_position_counts, control_position_counts,
					   nuclease_window_counts, control_window_counts,nuclease_window_converted_counts,nuclease_window_noise_counts]
				output_list.append(row)
	
	try:
		ecdf_pos = ECDF(bg_position)
		ecdf_nar = ECDF(bg_narrow)
	except Exception as e:
		print (e)
		ecdf_pos = ECDF([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
		ecdf_nar = ECDF([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
	
	# Genomic array to store the p-values for every chromosome:position object
	ga_pval = HTSeq.GenomicArray("auto", typecode='O', stranded=False)

	# Ratio to be used in scaling the nuclease count
	scale_factor = total_control_count/float(total_nuclease_count)
	logger.info("Finished forloop")
	logger.info(len(output_list))
	
	for idx, fields in enumerate(output_list):
		position_p_val = 1 - ecdf_pos(fields[2]*scale_factor)
		narrow_p_val = 1 - ecdf_nar(fields[4]*scale_factor)

		control_position_p_val = 1 - ecdf_pos(fields[3])
		control_narrow_p_val = 1 - ecdf_nar(fields[5])

		if fields[2] >= read_count_cutoff or fields[4] >= read_count_cutoff:  # fields[2] is nuclease_position_counts and fields[4] is nuclease_window_counts, 4 should be alwasy higher
			read_chr = fields[0]
			read_position = fields[1]
			offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
			ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[4] # this is the sum of number reads around a given position
			ga_narrow_windows_converted[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[6] # this is the sum of number reads around a given position
			ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[7] # this is the sum of number reads around a given position
		elif fields[4] >= edited_read_cutoff*2 and fields[6] >= edited_read_cutoff:
			# print ("using edited reads")
			read_chr = fields[0]
			read_position = fields[1]
			offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
			ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[4] # this is the sum of number reads around a given position
			ga_narrow_windows_converted[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[6] # this is the sum of number reads around a given position
			ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[7] # this is the sum of number reads around a given position

		chr_pos = HTSeq.GenomicPosition(fields[0], int(fields[1]), '.')
		# logger.info(chr_pos)
		ga_pval[chr_pos] = [position_p_val, narrow_p_val, control_position_p_val, control_narrow_p_val,fields[3],fields[5]]

	ga_consolidated_windows = find_windows(offtarget_ga_windows, window_size)	# consolidate windows within 3 bp
	logger.info(f"Start off-target sequence fuzzy alignment using mismatch cutoff of {mismatch_threshold}")
	output_alignments(ga_narrow_windows, ga_narrow_windows_converted,ga_narrow_windows_noise,ga_consolidated_windows, reference_genome_pyfaidx, targetsite, label, label, bam,
						mismatch_threshold, ga_pval, search_radius, f"{output_dir}/{label}",nuclease_ga_overlap,control_ga_overlap)

def main():
	parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
	parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
	parser.add_argument('--bam', help='Sorted BAM file', required=True)
	parser.add_argument('--control', help='Control BAM file', required=True)
	parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
	parser.add_argument('--search_radius', help='Search radius around the position window', default=20, type=int)
	parser.add_argument('--window_size', help='Windowsize', default=3, type=int)
	parser.add_argument('--mapq', help='mapq threshold', default=50, type=int)
	parser.add_argument('--gap', help='Gap threshold', default=3, type=int)
	parser.add_argument('--start', help='Start threshold', default=1 , type=int)
	parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=6, type=int)
	parser.add_argument('--read_count_cutoff', help='read_count threshold', default=6, type=int)
	parser.add_argument('--read_length', help='read_length', default=151, type=int)	
	parser.add_argument('--merged', dest='merged', action='store_true', default=True)
	parser.add_argument('--all_chromosomes', dest='all_chromosomes', action='store_true', default=False)
	parser.add_argument('--name', help='Targetsite Name', required=False)
	parser.add_argument('--cells', help='Cells', required=False)
	parser.add_argument('--out', help='Output file base', required=True)
	args = parser.parse_args()

	# Run the comparison if the control bam is specified, otherwise run the standard site identification routine.
	print("Nuclease: {0}\nControl: {1}".format(args.bam, args.control), file=sys.stderr)
	# compare(args.ref, args.bam, args.control, args.targetsite, args.search_radius, args.window_size, args.mapq, args.gap,
	# 		args.start, args.mismatch_threshold, args.name, args.cells, args.out, args.all_chromosomes, args.merged,args.read_count_cutoff)

if __name__ == "__main__":
	main()
