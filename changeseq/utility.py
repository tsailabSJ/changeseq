
import re
import gzip
import yaml
import validation
import logging
import string
import os
logger = logging.getLogger('root')
logger.propagate = False


def get_parameters(manifest_data):
	default_yaml = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"
	with open(default_yaml, 'r') as f:
		default = yaml.load(f)
	with open(manifest_data, 'r') as f:
		return_dict = yaml.load(f)
	default['analysis_folder'] = os.getcwd()
	validation.validateManifest(return_dict)
	for p in default:
		if not p in return_dict:
			if p == "samples":
				logger.error("No samples are found in the yaml file, please provide samples (fastq) to start with!")
				exit()
			return_dict[p] = default[p]
	return return_dict


def fq(file):
	if re.search('.gz$', file):
		fastq = gzip.open(file, 'rt')
	else:
		fastq = open(file, 'r')
	with fastq as f:
		while True:
			l1 = f.readline().strip()
			if not l1:
				break
			l2 = f.readline().strip()
			l3 = f.readline().strip()
			l4 = f.readline().strip()
			yield [l1, l2, l3, l4]


def reverseComplement(seq):
	tab = str.maketrans("ACTG", "TGAC")
	return seq.translate(tab)[::-1]