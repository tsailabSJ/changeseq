"""
alignReads
"""

import subprocess
import os
import logging
import sys

logger = logging.getLogger('root')
logger.propagate = False


def trim_reads(cutadapt=None, R1=None, R2=None, label=None,output_dir=None, njobs=None, **kwargs):
	"""Trimming PE reads and return trimmed file name"""
	logger.info(f'Start trimming for {label}')
	Tn5="CTGTCTCTTATACACATCTACGTAGATGTGTATAAGAGACAG"
	trimming_input = f"-a {Tn5} -A {Tn5}"
	trimmed_reads_PE_output = f"-o {output_dir}/{label}.trim.R1.fastq.gz -p {output_dir}/{label}.trim.R2.fastq.gz"
	# too_short_reads_PE_output = f"--too-short-output {output_dir}/{label}.too_short.R1.fastq.gz --too-short-paired-output {output_dir}/{label}.too_short.R2.fastq.gz"
	# untrimmed_reads_PE_output = f"--untrimmed-output {output_dir}/{label}.untrimmed.R1.fastq.gz --untrimmed-paired-output {output_dir}/{label}.untrimmed.R2.fastq.gz"
	# other_options = f"-j {njobs} -Z -m 10 --overlap 20 -e 0.1 --info-file {output_dir}/{label}.cutadapt_info_file.tsv --pair-filter both"
	other_options = f"-j {njobs} -Z --overlap 40 -e 0.15 --info-file {output_dir}/{label}.cutadapt_info_file.tsv"

	# cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {too_short_reads_PE_output} {untrimmed_reads_PE_output} {other_options} {R1} {R2}"
	cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {other_options} {R1} {R2}"
	logger.info(cutadapt_command)
	subprocess.call(cutadapt_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
	return f"{output_dir}/{label}.trim.R1.fastq.gz",f"{output_dir}/{label}.trim.R2.fastq.gz"

def alignReads(bwa=None, samtools=None, # Tools
			   reference_genome=None, R1=None, R2=None, label=None, output_dir=None, # Inputs
			   trim=None, trimmed_fastq_output = None, njobs=None, **kwargs): # other parameters
	"""Only for paired-end reads, return bam file name for next analysis

	Tools and inputs, e.g., reference_genome, R1, R2, can be absolute or relative path
	"""

	# output files
	bam_file = f"{output_dir}/{label}.bam" # tmp
	sorted_bam_file = f"{output_dir}/{label}.st.bam"
	name_sorted_bam_file = f"{output_dir}/{label}.st.name_sorted.bam"


	# Check if genome is already indexed by bwa
	index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']
	genome_indexed = True
	for extension in index_files_extensions:
		if not os.path.isfile(reference_genome + extension):
			genome_indexed = False
			break
	# If the genome is not already indexed, index it
	if not genome_indexed:
		logger.info('Genome index files not detected. Running BWA to generate indices.')
		bwa_index_command = f'{bwa} index {reference_genome}'
		logger.info('Running bwa command:', bwa_index_command)
		subprocess.call(bwa_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		logger.info('BWA genome index generated')

	# Trimming
	if trim:
		R1,R2 = trim_reads(R1=R1, R2=R2, label=label, output_dir=trimmed_fastq_output, njobs=njobs, **kwargs)

	# BWA alignment
	opts = f"-t {njobs}"
	bwa_alignment_command = f'{bwa} mem {opts} {reference_genome} {R1} {R2} | {samtools} view -bS - > {bam_file}'
	logger.info("BWA alignment")
	logger.info(bwa_alignment_command)
	subprocess.call(bwa_alignment_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)

	# samtools sort
	samtools_sort_command = f'{samtools} sort -@ {njobs} -o {sorted_bam_file} {bam_file}'
	samtools_index_command = f"{samtools} index {sorted_bam_file};rm {bam_file}"
	samtools_sort_name_command = f"{samtools} sort -@ {njobs} -o {name_sorted_bam_file} -n {sorted_bam_file};{samtools} index {name_sorted_bam_file}"

	# sort bam and index bam
	logger.info('samtools sort bam file')
	logger.info(samtools_sort_command)
	logger.info(samtools_index_command)
	logger.info(samtools_sort_name_command)
	subprocess.call(samtools_sort_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	subprocess.call(samtools_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	# subprocess.call(samtools_sort_name_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)

	return sorted_bam_file


def alignReads_bk(BWA_path, HG19_path, read1, read2, outfile):
	sample_name = os.path.basename(outfile).split('.')[0]
	output_folder = os.path.dirname(outfile)
	base_name = os.path.join(output_folder, sample_name)
	sam_filename = outfile
	bam_filename = base_name + '.bam'

	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	# Check if genome is already indexed by bwa
	index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']

	genome_indexed = True
	for extension in index_files_extensions:
		if not os.path.isfile(HG19_path + extension):
			genome_indexed = False
			break

	# If the genome is not already indexed, index it
	if not genome_indexed:
		logger.info('Genome index files not detected. Running BWA to generate indices.')
		bwa_index_command = '{0} index {1}'.format(BWA_path, HG19_path)
		logger.info('Running bwa command: %s', bwa_index_command)
		subprocess.call(bwa_index_command.split())
		logger.info('BWA genome index generated')
	else:
		logger.info('BWA genome index found.')

	# Run paired end alignment against the genome
	logger.info('Running paired end mapping for {0}'.format(sample_name))
	# bwa_alignment_command = '{0} mem -t 6 -B 2 -L 9999,5 {1} {2} {3} > {4}'.format(BWA_path, HG19_path, read1, read2, sam_filename)
	bwa_alignment_command = '{0} mem -t 6 {1} {2} {3} > {4}'.format(BWA_path, HG19_path, read1, read2, sam_filename)
	samtools_sam_to_bam_command = 'samtools sort -o {0} {1}'.format(bam_filename, sam_filename)
	samtools_index_command = 'samtools index {0}'.format(bam_filename)
	samtools_sort_by_name_command = 'samtools sort -o {0} -n {1}'.format("".join([base_name, '_sorted.bam']),
																		 bam_filename)

	# Open the outfile and redirect the output of the alignment to it.
	logger.info(bwa_alignment_command)
	subprocess.call(bwa_alignment_command, shell=True)
	logger.info('Paired end mapping for {0} completed.'.format(sample_name))

	# Convert SAM to BAM file
	logger.info(samtools_sam_to_bam_command)
	subprocess.call(samtools_sam_to_bam_command, shell=True)
	logger.info('Sorting by coordinate position for {0} complete.'.format(sample_name))

	# Index BAM file
	logger.info(samtools_index_command)
	subprocess.call(samtools_index_command, shell=True)
	logger.info('Indexing for {0} complete.'.format(sample_name))

	# Sort BAM file by name
	logger.info(samtools_sort_by_name_command)
	subprocess.call(samtools_sort_by_name_command, shell=True)
	logger.info('Sorting for {0} by name complete.'.format(sample_name))


def alignReads_bowtie2(BWA_path, index_path, read1, read2, outfile):
	sample_name = os.path.basename(outfile).split('.')[0]
	output_folder = os.path.dirname(outfile)
	base_name = os.path.join(output_folder, sample_name)
	sam_filename = outfile
	bam_filename = base_name + '.bam'

	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	# Run paired end alignment against the genome
	logger.info('Running paired end mapping for {0}'.format(sample_name))
	bowtie_alignment_command = 'bowtie2 -p 6 --very-sensitive -x {0} -1 {1} -2 {2} -S {3}'.format(index_path, read1,
																								  read2, sam_filename)
	samtools_sam_to_bam_command = 'samtools sort -o {0} {1}'.format(bam_filename, sam_filename)
	samtools_index_command = 'samtools index {0}'.format(bam_filename)
	samtools_sort_by_name_command = 'samtools sort -o {0} -n {1}'.format("".join([base_name, '_sorted.bam']),
																		 bam_filename)

	# Open the outfile and redirect the output of the alignment to it.
	logger.info(bowtie_alignment_command)
	subprocess.call(bowtie_alignment_command, shell=True)
	logger.info('Paired end mapping for {0} completed.'.format(sample_name))

	# Convert SAM to BAM file
	logger.info(samtools_sam_to_bam_command)
	subprocess.call(samtools_sam_to_bam_command, shell=True)
	logger.info('Sorting by coordinate position for {0} complete.'.format(sample_name))

	# Index BAM file
	logger.info(samtools_index_command)
	subprocess.call(samtools_index_command, shell=True)
	logger.info('Indexing for {0} complete.'.format(sample_name))

	# Sort BAM file by name
	logger.info(samtools_sort_by_name_command)
	subprocess.call(samtools_sort_by_name_command, shell=True)
	logger.info('Sorting for {0} by name complete.'.format(sample_name))
