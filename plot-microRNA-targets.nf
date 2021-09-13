#!/usr/bin/env nextflow

/*================================================================
The Aguilar Lab presents...

The miRNome comparisson plot Pipeline

- A miRNA pairs comparisson plot tool

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)


- Bioinformatics Development
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

=============================
Pipeline Processes In Brief:

Pre-processing:
_pre1_compare-mirnatargets
_pre2_convert-target-file

Core-processing:
_001_butterfly-plot-target-changes
_002_plot-target-changes-count

Pos-processing:

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
	The miRNome comparisson plot Pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --targets_ref <path to input 1> --targets_mut <path to input 2> [--output_dir path to results ]

	  --targets_ref	<- miRNA reference targets file;

	  --targets_mut	<- miRNA mutate targets file;

	  --output_dir     <- directory where results, intermediate and log files will bestored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "plot-microRNA-targets"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.targets_ref = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.targets_mut = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.targets_ref | !params.targets_mut ) {
  log.error " Please provide the --targets_ref AND --targets_mut \n\n" +
  " For more information, execute: nextflow run compare-miRNA-pairs.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.targets_ref).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The miRNome comparisson plot Pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input reference targets file']			= params.targets_ref
pipelinesummary['Input mutate targets file']			= params.targets_mut
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/*
	READ INPUTS_A
*/

/* Load  miRNA reference targets file into channel*/
Channel
	.fromPath( "${params.targets_ref}" )
	// .view()
	.set{ targets_ref_input }

/* Load miRNA mutate targets file into channel */
Channel
	.fromPath( "${params.targets_mut}" )
	// .view()
	.set{ targets_mut_input }

/* _pre1_compare_mirnatargets */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mk-modules/mk-compare-mirnatargets/*")
	.toList()
	.set{ mkfiles_pre1 }

process _pre1_compare_mirnatargets {

	publishDir "${results_dir}/_pre1_compare_mirnatargets/",mode:"copy"

	input:
	file ref_targets from targets_ref_input
	file mut_targets from targets_mut_input
	file mk_files from mkfiles_pre1

	output:
	file "*.changes" into results_pre1_compare_mirnatargets
	file "*.png" into results_pre1_compare_mirnatargets_png
	"""
	bash runmk.sh
	"""

}


/* _pre2_convert-target-file */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mk-modules/mk-convert-target-file/*")
	.toList()
	.set{ mkfiles_pre2 }

process _pre2_convert_target_file {

	publishDir "${results_dir}/_pre2_convert_target_file/",mode:"copy"

	input:
	file changes from results_pre1_compare_mirnatargets
	file mk_files from mkfiles_pre2

	output:
	file "*.changes.tsv" into results_A_pre2_convert_target_file, results_B_pre2_convert_target_file

	"""
	bash runmk.sh
	"""

}

/* 001_butterfly-plot-target-changes */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mk-modules/mk-butterfly-plot-target-changes/*")
	.toList()
	.set{mkfiles_core1}


process _001_butterfly_plot_target_changes {

	publishDir "${results_dir}/001_butterfly_plot_target_changes/",mode:"copy"

	input:
	file changes_tsv from results_A_pre2_convert_target_file
	file mk_files from mkfiles_core1

	output:
	file "*.png" into results_001_butterfly_plot_target_changes

	"""
	bash runmk.sh
	"""

}

/* 001_butterfly-plot-target-changes */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mk-modules/mk-plot-target-changes-count/*")
	.toList()
	.set{mkfiles_core2}


process _002_plot_target_changes_count {

	publishDir "${results_dir}/_002_plot-target-changes-count/",mode:"copy"

	input:
	file changes_tsv from results_B_pre2_convert_target_file
	file mk_files from mkfiles_core2

	output:
	file "*.png" into results_002_plot_target_changes_count

	"""
	bash runmk.sh
	"""

}
