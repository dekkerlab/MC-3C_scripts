# Cwalk analysis pipeline README

Author: Filipe Tavares-Cadete

## Introduction

The pipeline to analyse the Dekker lab Cwalk data consists of several steps:
	1) Processing of raw PacBio data into fastq files;
	2) Processing of fastq files to separate into interaction fragments;
	3) Mapping of interaction fragments;
	4) Assembly of alignments into walks;
	5) Preparing data frames with detailed walk information;
	6) Preparing walk permutations;
	7) Scripts for plotting.

## Step requirements

All steps can be achieved on a Unix environment on a normal workstation, unless specifically noted.

### 1) Processing of raw PacBio data into fastq files;

This step requires the SMRT Analysis software by Pacific Biosystems running on a Unix environment.

## 2) Processing of fastq files to separate into interaction fragments;

This step uses the 'digest_roi.py' script and requires Python 2.7 with the Bio package installed.

## 3) Mapping of interaction fragments

This step requires bwa-mem version 0.7.12 and samtools version 1.3 installed. Exact parameters are found on 'launch_bwa_mem.sh'. For faster run-time, a machine with a large number of cores (32 or above) and large memory (32Gb or above) is recommended.

## Assembly of alignments into walks

This step is done with the 'reduce_frag_mappings.R' script, running R 3.5.0 or later, with the BioConductor GenomicRanges package installed.

## 5) Preparing data frames with detailed walk information

This step is done with the 'interactions_to_usable_frame_stricter.R' and 'interactions_to_usable_frame_keep_NAs.R' scripts. They require R 3.5.0 or later, with the GenomicRanges, rtracklayer, and tidyverse packages installed.

## 6) Preparing walk permutations

This step is done through the 'launch_permutations.sh' script. For faster results the use of a machine with 32 cores and 64Gb of RAM is recommended.

## 7) Scripts for plotting

Plotting was done in R, version 3.5.0 or later, with the tidyverse, cowplot and gridExtra packaged installed. 