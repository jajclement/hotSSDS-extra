#!/bin/bash

#################################################################################
#		BASH WRAPPER FOR SSDSPOSTPROCESS PIPELINE			#
#################################################################################
# 2020, Pauline Auffret
# This bash wrapper is designed to launch ssdsnextflowpipeline.
# See https://gitlab.igh.cnrs.fr/pauline.auffret/ssdsnextflowpipeline/-/tree/master

#Set slurm environment variables
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

#Singularity settings
#Set environmental variable to mount the non canonical directories on IGH cluster
#Not used currently (as of 2021-12-17)
export SINGULARITY_BINDPATH="/work,/poolzfs"

#Pipeline default parameters
ANALYSIS_NAME="SSDS_postprocess"
PIPELINE_DIRECTORY="/home/${USER}/work/ssdspostprocess"
BASE_DIRECTORY="/home/${USER}/work/results"
CONF="${PIPELINE_DIRECTORY}/conf/igh.config"
GENOME_PROFILE="${PIPELINE_DIRECTORY}/conf/mm10.json"
CENV="/home/${USER}/work/bin/miniconda3/envs/nextflow_dev"
now=`date +"%FT%H%M%S"`
INPUT=""
OPTIONS=""
TEST="0"
FORCE="0"
TOWER_TOKEN="None"

#Get command line arguments
while getopts hp:b:n:c:a:o:t:f:g: flag
do
	case "${flag}" in
		h) echo ""; echo "Usage: bash `basename $0` -i input_file [options] "; \
		   echo "Options : "; echo "-h display help message"; \
		   echo "-g Absolute path to the genome config file (default : ${GENOME_PROFILE})" ; \
		   echo "-p Absolute path to ssds postprocess pipeline base directory (default : ${PIPELINE_DIRECTORY})"; \
		   echo "-b Absolute path to base directory where the output directory will be created (default : ${BASE_DIRECTORY})"; \
		   echo "-n Analysis name (default : ${ANALYSIS_NAME}) INFO : by default, this parameter will match the --name option in nextflow command line"; \
		   echo "-c Absolute path to IGH cluster configuration file (default : ${CONF})"; \
		   echo "-a Absolute path to conda environment for nextflow (default : ${CENV})"; \
		   echo "-o Optional arguments for the pipeline (for example \"--corMethod spearman --with_heatmap false\" ;  default : \"${OPTIONS}\")"; \
		   echo "-w Valid Nextflow Tower token (default : ${TOWER_TOKEN} ; if not None, then the option -with-tower has to be added in -o parameter)"; \
		   echo "-t set to 1 if running pipeline on test data located in ${PIPELINE_DIRECTORY}/tests/fastq (default : ${TEST})"; \
		   echo "-f set to 1 to force pipeline to run without checking resume/output directory (default : ${FORCE})" ; \
		   echo "INFO : the output directory will be located in the base directory and will be named after the analysis name parameter with the .outdir suffix (default ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir)"; \
		   echo ""; exit 0;;
		p) PIPELINE_DIRECTORY=${OPTARG};if [ ! -d ${PIPELINE_DIRECTORY} ]; then echo "Directory ${PIPELINE_DIRECTORY} not found!" ; exit 0; fi;;
		b) BASE_DIRECTORY=${OPTARG};;
		n) ANALYSIS_NAME=${OPTARG};;
		g) GENOME_PROFILE=${OPTARG};if [ ! -f ${GENOME_PROFILE} ]; then echo "File ${GENOME_PROFILE} not found!" ; exit 0; fi;;
		c) CONF=${OPTARG};if [ ! -f ${CONF} ]; then echo "File ${CONF} not found!" ; exit 0; fi;;
		a) CENV=${OPTARG};if [ ! -d ${CENV} ]; then echo "Environment ${CENV} not found!" ; exit 0; fi;;
		o) OPTIONS=${OPTARG};;
		w) TOWER_TOKEN=${OPTARG};;
		t) TEST=${OPTARG};;
		f) FORCE=${OPTARG};;
		\? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        	:  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        	*  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
	esac
done

#Create base directory if not existing
echo "Checking base directory..."
mkdir -p ${BASE_DIRECTORY}

#Test if a project with the same name already exists and if the option -resume is properly set
if [[ $FORCE == 0 ]] 
then
	echo "Checking output directory."
	if [[ -d ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir && $OPTIONS == *"-resume"* ]]
	then
		echo "WARNING : Output directory already exists for project ${ANALYSIS_NAME} in ${BASE_DIRECTORY}, are you sure you want to resume this run ? (y/n)"
		read answer1
		case $answer1 in
			[yYoO]*) echo "Ok; run will be resumed.";;
			[nN]*) echo "Ok; quitting. Bye, see you soon !"; exit 0;;
			*) echo "ABORT. Please enter y (yes) or n (no) next time. Bye, see you soon !"; exit 1;;
		esac
	elif [[ -d ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir && ! $OPTIONS == *"-resume"* ]]
	then
		echo "WARNING : Output directory already exists for project ${ANALYSIS_NAME} in ${BASE_DIRECTORY}, are you sure you want to run the pipeline from the beginning (this will erase previous results from this project) ? (y/n)"
		read answer2
		case $answer2 in
			[yYoO]*) echo "Ok; previous results located in ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir will be erased.";;
			[nN]*) echo "Ok; quitting. Consider using another analysis name (-n argument) or using the option -resume in the -o argument next time. Bye, see you soon !"; exit 0;;
			*) echo "ABORT. Please enter y (yes) or n (no) next time. Bye, see you soon !"; exit 1;;
		esac
	else
		echo "Ok, creating ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir output project directory."
		mkdir ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir
	fi
else
	echo "Forcing the creation of ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir output project directory."
	mkdir -p ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir
fi

#Go to the output directory
cd ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir

#Create output directory for slurm log files
mkdir -p ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir/slurm

#Settings if running in tst mode
if [ ${TEST} == "0" ]; then
	echo "Running SSDS postprocess pipeline from ${PIPELINE_DIRECTORY} data within ${CENV##*/} conda environment. Check output directory ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir/"
else 
	GENOME_PROFILE=${PIPELINE_DIRECTORY}/tests/mm10_test.json
	echo "Running SSDS postprocess pipeline from ${PIPELINE_DIRECTORY} on test data within ${CENV##*/} conda environment. Check output directory ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir/"
fi 

#Check conda environment.
echo "Checking conda environment..."
if [ ! -d ${CENV} ]; 
then echo "ABORT : environment ${CENV} not found! Check -a argument. Bye, see you soon !" ; exit 0; 
else echo "Ok.";
fi

#Check pipeline directory
echo "Checking pipeline directory..."
if [ ! -d ${PIPELINE_DIRECTORY} ];
then echo "ABORT : pipeline directory ${PIPELINE_DIRECTORY} not found! Check -p argument. Bye, see you soon !" ; exit 0;
else echo "Ok.";
fi

#Checking configuration files
echo "Checking configuration files..."
if [ ! -f ${CONF} ] ;
then echo "ABORT : configuration file ${CONF} not found! Check -c argument. Bye, see you soon !" ; exit 0;
else echo "Ok.";
fi
if [ ! -f ${GENOME_PROFILE} ];
then echo "ABORT : configuration file ${GENOME_PROFILE} not found! Check -g argument. Bye, see you soon !" ; exit 0; 
else echo "Ok.";
fi 

#Activate conda environment
echo "Activate conda environment ${CENV}."
eval "$(conda shell.bash hook)"
conda activate ${CENV}

#If tower token is set, then TOWER_ACCESS_TOKEN variable must be added to the environment
if [ ${TOWER_TOKEN} != "None" ];
then
    export TOWER_ACCESS_TOKEN=${TOWER_TOKEN} ;
    export NXF_VER=21.10.0 ;
fi


JOBNAME="SSDS_main_${ANALYSIS_NAME}_${now}"

sbatch -p computepart -J ${JOBNAME} -o ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir/slurm/%x.%j.out --export=ALL -n 1 --mem 7G -t 5-0:0  \
--wrap "export MKL_NUM_THREADS=1 ; export NUMEXPR_NUM_THREADS=1 ; export OMP_NUM_THREADS=1 ; \
nextflow run ${PIPELINE_DIRECTORY}/main.nf -c ${CONF} -params-file ${GENOME_PROFILE} --name ${ANALYSIS_NAME} --outdir ${BASE_DIRECTORY}/${ANALYSIS_NAME}.outdir ${OPTIONS}"

#Deactivate conda environment
conda deactivate

