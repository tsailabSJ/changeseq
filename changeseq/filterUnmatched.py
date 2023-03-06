from __future__ import print_function


import os
import logging
import argparse
import pandas as pd
import subprocess
from skbio.alignment import global_pairwise_align_nucleotide 
from skbio.sequence import DNA
logger = logging.getLogger('root')
logger.propagate = False

"""Parse identified unmatched table and relax mismatch threshold

The goal is to perform a sensitive local alignment and output possible off-target-gRNA alginment



"""



def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	mainParser.add_argument('-o',"--output", help="output file", default="off-target_high_mismatch.tsv")

	mainParser.add_argument('-f',"--input",  help="identified unmatched txt",required=True)
	mainParser.add_argument('--PAM',  help="identified unmatched txt",required=True)
	mainParser.add_argument('--gRNA',  help="gRNA without PAM",required=True)
	mainParser.add_argument('--min_count',  help="identified unmatched txt",required=True)


	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def revcomp(seq):
	tab = str.maketrans("ACTG", "TGAC")
	return seq.translate(tab)[::-1]
		
def get_edited_reads(df):
	col=30
	return [int(x.split(",")[0]) for x in df[col]]

def force_realign(input,gRNA,PAM,min_count,output):
	df = pd.read_csv(input,sep="\t")
	df['edit'] = get_edited_reads(df)
	df = df[df.edit>=min_count]
	off_Target_seq_col = 9
	df[["Site_Sequence","Site_Substitution_Number","Site_Sequence_Gaps_Allowed","Target_Sequence","Realigned_Target_Sequence","Strand"]] = df.apply(lambda x:pairwise_alginment_both_strand(gRNA,x[9],PAM), axis=1).apply(pd.Series)
	output_to_vis(df,output)

def pairwise_alginment_both_strand(gRNA,off_target,PAM):
	off_target_revcomp = revcomp(off_target)
	Site_Sequence,Site_Substitution_Number,Site_Sequence_Gaps_Allowed,Target_Sequence,Realigned_Target_Sequence,score = pairwise_alginment(gRNA,off_target,PAM)
	Site_Sequence_revcomp,Site_Substitution_Number_revcomp,Site_Sequence_Gaps_Allowed_revcomp,Target_Sequence_revcomp,Realigned_Target_Sequence_revcomp,score_revcomp = pairwise_alginment(gRNA,off_target,PAM)
	if score> score_revcomp:
		return [Site_Sequence,Site_Substitution_Number,Site_Sequence_Gaps_Allowed,Target_Sequence,Realigned_Target_Sequence,"+"]
	else:
		return [Site_Sequence_revcomp,Site_Substitution_Number_revcomp,Site_Sequence_Gaps_Allowed_revcomp,Target_Sequence_revcomp,Realigned_Target_Sequence_revcomp,"-"]


def pairwise_alginment(gRNA,off_target,PAM):
	Site_Sequence = ""
	Site_Substitution_Number = -1
	Site_Sequence_Gaps_Allowed = ""
	Target_Sequence = gRNA+ PAM
	Realigned_Target_Sequence = ""
	alignment, score, start_end_positions = global_pairwise_align_nucleotide(DNA(gRNA),DNA(off_target))
	a,b = alignment[0]._string.decode("utf-8"),alignment[1]._string.decode("utf-8")
	aligned_gRNA = []
	aligned_off_target = []
	count = 0
	for i in range(len(a)):
		if a[i]==b[i]=="-":
			continue
		if a[i]!=b[i]:
			count+=1
		aligned_gRNA.append(a[i])
		aligned_off_target.append(b[i])
	aligned_off_target = "".join(aligned_off_target)
	current_PAM_pos = b.index(aligned_off_target)+len(aligned_off_target)
	current_PAM = b[current_PAM_pos:current_PAM_pos+len(PAM)]
	Site_Sequence = aligned_off_target.replace("-","")+current_PAM
	Site_Substitution_Number = count
	if "-" in aligned_gRNA or "-" in aligned_off_target:
		Site_Sequence_Gaps_Allowed = aligned_off_target+current_PAM
		Realigned_Target_Sequence = "".join(aligned_gRNA) + PAM
	return Site_Sequence,Site_Substitution_Number,Site_Sequence_Gaps_Allowed,Target_Sequence,Realigned_Target_Sequence,score
		
def output_to_vis(df,output):
	df['Control_Read_Count'] = 0
	columns_to_keep = [0,1,2,3,4,"Strand",'Control_Read_Count','Site_Sequence','Site_Substitution_Number','Site_Sequence_Gaps_Allowed',25,26,27,'Target_Sequence','Realigned_Target_Sequence',30,31]
	df = df[columns_to_keep]
	rename_columns = ['#Chromosome','Start','End','Genomic Coordinate','Nuclease_Read_Count','Strand','Control_Read_Count','Site_Sequence','Site_Substitution_Number','Site_Sequence_Gaps_Allowed','File_Name','Cell','Target_site','Full_Name','Target_Sequence','Realigned_Target_Sequence','Nuclease_overlap_bp_list','Control_overlap_bp_list']
	df.columns = rename_columns
	df.to_csv(output,sep="\t",index=False)

def main():

	args = my_args()
	force_realign(args.input,args.gRNA,args.PAM,args.min_count,args.output)


if __name__ == "__main__":
	main()



