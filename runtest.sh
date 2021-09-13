targets_ref="test/data/sample.targets.ref"
targets_mut="test/data/sample.targets.mut"
output_directory="$(dirname $targets_ref)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run plot-microRNA-targets.nf \
	--targets_ref $targets_ref \
  --targets_mut $targets_mut \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
