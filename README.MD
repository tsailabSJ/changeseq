
[![Version][version-shield]][version-url]
[![Python versions][python-shield]][python-url]
[![Platforms][platform-shield]][python-url]


# CHANGE-seq-BE enables simultaneously sensitive and unbiased in vitro profiling of base editor genome-wide activity

This is a repository for CHANGE-seq-BE analytical software, which takes sample-specific paired-end FASTQ files as input and produces a list of CHANGE-seq-BE detected off-target cleavage sites as output.

# Note

Currently, code for CHANGE-seq and CHANGE-seq-BE are different. `CHANGE-seq-BE` is in the `dev` branch only.

# Summary

This package implements a pipeline that takes in reads from the CHANGE-seq-BE assay and returns detected cleavage sites as output. The individual pipeline steps are:

1. **Tn5 adapter trimming**: 
2. **Read Alignment**: Merged paired end reads from the assay are aligned to the reference genome using the BWA-MEM algorithm with default parameters (Li. H, 2009).
3. **Cleavage Site Identification**: Mapped sites are analyzed to determine which represent high-quality cleavage sites.
4. **Visualization of Results**: Identified on-target and off-target cleavage sites are rendered as a color-coded alignment map for easy analysis of results.

# Installation

The most easiest way to install change-seq pipeline is via conda.

```

mamba create -n changeseqbe tsailabsj::changeseqbe

mamba activate changeseqbe

git clone -b dev https://github.com/tsailabSJ/changeseq

python changeseq/changeseq/changeseq.py -h

## BWA 0.7.17 and samtools 1.9 are automatically installed

```


# Usage

The change-seq-be pipeline requires a manifest yaml file specifying input files, output directory, and pipeline parameters. Once the yaml file is created, users can simply run ``change_seq.py all --manifest /path/to/manifest.yaml``


Below is an example ``manifest.yaml`` file::

	reference_genome: /research/groups/tsaigrp/projects/Genomics/common/genomes/hg38/hg38_chroms_only.fa
	analysis_folder: /research/groups/tsaigrp/projects/Genomics/common/projects/CHANGE_seq/CHANGE_seq_BE/test

	bwa: bwa
	samtools: samtools

	read_threshold: 3
	window_size: 30
	mapq_threshold: 50
	start_threshold: 1
	gap_threshold: 30
	mismatch_threshold: 6
	search_radius: 30
	merged_analysis: False
	BEmodel_min_overlap: -30
	BEmodel_max_overlap: 30
	BE_analysis: True

	PAM: NGG
	genome: hg38

	samples:
	    test:
	        target: sgRNA-seq-NGG
	        read1: GM255_S2_R1_001.fastq.gz
	        read2: GM255_S2_R2_001.fastq.gz
	        controlread1: GM254_S1_R1_001.fastq.gz
	        controlread2: GM254_S1_R2_001.fastq.gz
	        description: GM255_CD40L_E1

## Quickstart

(TODO)

```

git clone https://github.com/tsailabSJ/changeseq

cd changeseq/test

changeseq.py all --manifest CIRCLEseq_MergedTest.yaml

```

# Writing A Manifest File
When running the end-to-end analysis functionality of the CHANGEseq package a number of inputs are required. To simplify the formatting of these inputs and to encourage reproducibility, these parameters are inputted into the pipeline via a manifest formatted as a YAML file. YAML files allow easy-to-read specification of key-value pairs. This allows us to easily specify our parameters. The following fields are required in the manifest:

- `reference_genome`: The absolute path to the reference genome FASTA file.
- `output_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- `bwa`: The absolute path to the `bwa` executable
- `samtools`: The absolute path to the `samtools` executable
- `read_threshold`: The minimum number of reads at a location for that location to be called as a site. We recommend leaving it to the default value of 4.
- `window_size`: Size of the sliding window, we recommend leaving it to the default value of 3.
- `mapq_threshold`: Minimum read mapping quality score. We recommend leaving it to the default value of 50.
- `start_threshold`: Tolerance for breakpoint location. We recommend leaving it to the default value of 1.
- `gap_threshold`: Distance between breakpoints. We recommend leaving it to the default value of 3 for Cas9.
- `mismatch_threshold`: Number of tolerated gaps in the fuzzy target search setp. We recommend leaving it to the default value of 6.
- `read_length`: Fastq file read length, default is 151.
- `PAM`: PAM sequence, default is NGG.
- `merged_analysis`: Whether or not the paired read merging step should takingTrue
- `samples`: Lists the samples you wish to analyze and the details for each. Each sample name should be nested under the top level samples key, and each sample detail should be nested under the sample name. See the sample manifest for an example.
    - For each sample, you must provide the following parameters:
        - `target`: Target sequence for that sample. Accepts degenerate bases.
        - `read1`: The absolute path to the .FASTQ(.gz) file containing the read1 reads.
        - `read2`: The absolute path to the .FASTQ(.gz) file containing the read2 reads.
        - `controlread1`: The absolute path to the .FASTQ(.gz) file containing the control read1 reads.
        - `controlread2`: The absolute path to the .FASTQ(.gz) file containing the control read2 reads.
        - `description`: A brief description of the sample


# Pipeline Output
When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:

- `output_folder/aligned`: Contains an alignment `.sam`, alignment `.bam`, sorted `bam`, and `.bai` index file for each sample.
- `output_folder/fastq`: Merged `.fastq.gz` files for each sample.
- `output_folder/identified`: Contains tab-delimited `.txt` files for each sample containing the identified DSBs, control DSBs, filtered DSBs, and read quantification.
- `output_folder/visualization`: Contains a `.svg` vector image representing an alignment of all detected off-targets to the targetsite for each sample.

# FAQ

None yet, we will keep this updated as needed.

[version-shield]: https://img.shields.io/conda/v/tsailabsj/changeseq.svg
[version-url]: https://anaconda.org/tsailabSJ/changeseq
[python-shield]: https://img.shields.io/pypi/pyversions/changeseq.svg
[python-url]: https://pypi.python.org/pypi/changeseq
[platform-shield]: https://anaconda.org/tsailabsj/changeseq/badges/platforms.svg
