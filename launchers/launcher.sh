#!/usr/bin/env bash
cd ../


targets_ref="real-data/*.ref"
targets_mut="real-data/*.mut"
output_directory="real-data/"

nextflow run compare-miRNA-pairs.nf \
	--targets_ref $targets_ref \
	--targets_mut	$targets_mut \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html
