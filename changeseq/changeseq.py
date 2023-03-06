#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
circleseq.py as the wrapper for CIRCLE-seq analysis
"""

from alignReads import alignReads
from visualization import visualizeOfftargets
import argparse
import os
import sys
import subprocess
import traceback
import log
import findCleavageSites
from utility import get_parameters
from copy import deepcopy as dp
logger = log.createCustomLogger('root')

class CircleSeq:

	def __init__(self):
		self.parameters = {}
		self.output_dir = {}
		self.findCleavageSites_input_bam = {}
		self.vis_input_tsv = {}
	def parseManifest(self, manifest_path, sample='all'):
		logger.info('Loading manifest...')
		try:
			parameters = get_parameters(manifest_path)
			self.parameters = dp(parameters)
			# print (parameters)
			if sample != 'all':
				self.parameters['samples'] = {}
				self.parameters['samples'][sample] = parameters['samples'][sample]
			# print (self.parameters)
			# Make folders for output
			for folder in ['aligned', 'identified', 'fastq', 'visualization']:
				self.output_dir[folder] = os.path.join(self.parameters["analysis_folder"], folder)
				if not os.path.exists(self.output_dir[folder]):
					os.makedirs(self.output_dir[folder])
			# Just to initialize some default input file names for the identify and visualization steps
			for sample in self.parameters['samples']:
				self.findCleavageSites_input_bam[sample] = [f"{self.output_dir['aligned']}/{sample}.st.bam",f"{self.output_dir['aligned']}/Control_{sample}.st.bam"]
			# print (self.findCleavageSites_input_bam)
		except Exception as e:
			logger.error(
				'Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
			logger.error(traceback.format_exc())
			sys.exit()

	def alignReads(self):
		"""BWA mapping
		"""
		logger.info('Aligning reads...')
		for sample in self.parameters['samples']:
			try:
				bam=alignReads(trim=True, output_dir=self.output_dir['aligned'], trimmed_fastq_output=self.output_dir['fastq'],
						   R1=self.parameters['samples'][sample]["read1"], R2=self.parameters['samples'][sample]["read2"], label=sample,
						   **self.parameters)
				control_bam = alignReads(trim=True, output_dir=self.output_dir['aligned'], trimmed_fastq_output=self.output_dir['fastq'],
						   R1=self.parameters['samples'][sample]["controlread1"], R2=self.parameters['samples'][sample]["controlread2"], label="Control_"+sample,
						   **self.parameters)
				logger.info('Finished aligning reads to genome.')
				self.findCleavageSites_input_bam[sample] = [bam,control_bam]
			except Exception as e:
				logger.error('Error aligning for sample %s.'%(sample))
				logger.error(traceback.format_exc())

	def findCleavageSites(self):
		logger.info('Identifying off-target cleavage sites...')
		for sample in self.parameters['samples']:
			try:
				# logger.info(sample)
				bam,control_bam = self.findCleavageSites_input_bam[sample]
				findCleavageSites.compare(bam=bam, control=control_bam, label=sample, output_dir=self.output_dir['identified'],
							targetsite=self.parameters['samples'][sample]["target"],**self.parameters)
			except Exception as e:
				logger.error('Error findCleavageSites for sample %s.'%(sample))
				logger.error(traceback.format_exc())

	def visualize(self):
		logger.info('Visualizing off-target sites')

		for sample in self.parameters['samples']:
			try:
				infile = os.path.join(self.parameters["analysis_folder"], 'identified', sample + '_identified_matched.txt')
				outfile = os.path.join(self.parameters["analysis_folder"], 'visualization', sample + '_offtargets')
				visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters["PAM"])
			except Exception as e:
				logger.error('Error visualizing off-target sites: %s'%(sample))
				logger.error(traceback.format_exc())
		logger.info('Finished visualizing off-target sites')

	def parallel(self, manifest_path, lsf, run='all'):
		logger.info('Submitting parallel jobs')
		current_script = __file__

		try:
			for sample in self.parameters['samples']:
				cmd = ' python {0} {1} --manifest {2} --sample {3}'.format(current_script, run, manifest_path, sample)
				logger.info(cmd)
				# subprocess.call(lsf.split() + [cmd])
				# print (lsf+cmd)
				subprocess.call(lsf + cmd,shell=True)
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())

	def skip_align_parallel(self, manifest_path, lsf, run='skip_align'):
		logger.info('Submitting parallel jobs')
		current_script = __file__

		try:
			for sample in self.parameters['samples']:
				cmd = 'python {0} {1} --manifest {2} --sample {3}'.format(current_script, run, manifest_path, sample)
				logger.info(cmd)
				subprocess.call(lsf.split() + [cmd])
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())

def parse_args():
	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(description='Individual Step Commands',
									   help='Use this to run individual steps of the pipeline',
									   dest='command')

	all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
	all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	skip_align_parser = subparsers.add_parser('skip_align', help='Run all steps of the pipeline')
	skip_align_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	skip_align_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)',required=True)

	parallel_parser = subparsers.add_parser('parallel', help='Run all steps of the pipeline in parallel')
	parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -n 6 -R "rusage[mem=20000] span[hosts=1]" -P ABE -q priority ')
	parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='all')

	skip_align_parallel_parser = subparsers.add_parser('skip_align_parallel', help='Run all steps of the pipeline in parallel')
	skip_align_parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	skip_align_parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=100000] -P ABE -q priority')
	skip_align_parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='skip_align')

	align_parser = subparsers.add_parser('align', help='Run alignment only')
	align_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	align_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	merge_parser = subparsers.add_parser('merge', help='Merge paired end reads')
	merge_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	merge_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	identify_parser = subparsers.add_parser('identify', help='Run identification only')
	identify_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	identify_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	visualize_parser = subparsers.add_parser('visualize', help='Run visualization only')
	visualize_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	visualize_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	variants_parser = subparsers.add_parser('variants', help='Run variants analysis only')
	variants_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	variants_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	reference_free_parser = subparsers.add_parser('reference-free', help='Run reference-free discovery only')
	reference_free_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	reference_free_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	return parser.parse_args()

def main():
	args = parse_args()

	if args.command == 'all':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.alignReads()
		c.findCleavageSites()
		c.visualize()
		# c.callVariants()
		# c.extract_deamination_position()
		# c.extract_outward_reads()
	elif args.command == 'skip_align':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.findCleavageSites()
		c.visualize()
		# c.callVariants()
	elif args.command == 'parallel':
		c = CircleSeq()
		c.parseManifest(args.manifest)
		c.parallel(args.manifest, args.lsf, args.run)
	elif args.command == 'skip_align_parallel':
		c = CircleSeq()
		c.parseManifest(args.manifest)
		c.skip_align_parallel(args.manifest, args.lsf, args.run)
	elif args.command == 'align':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.alignReads()
	elif args.command == 'identify':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.findCleavageSites()
		c.visualize()
	elif args.command == 'merge':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.mergeAlignReads()
	elif args.command == 'visualize':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.visualize()
	elif args.command == 'variants':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.callVariants()

if __name__ == '__main__':
	main()
