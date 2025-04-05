#!/bin/bash
#SBATCH -p medium
#SBATCH -t 1-00:00
#SBATCH --mem=200
#SBATCH -c 1
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL_ADDRESS

##usage: 
# bash run_pipeline.sh [single|pooled] SAMPLE_ID
# sbatch -o LOG_DIR/SAMPLE_ID.%j run_pipeline.sh [single|pooled] SAMPLE_ID

MODE=$1
SAMPLE_ID=$2

if [ ${MODE} == "single" ]; then
	snakefile=scripts/Snakefile.single_cell
	profile=scripts/slurm.single_cell
else if [ ${MODE} == "pooled" ]; then
	snakefile=scripts/Snakefile.pooled_cell
	profile=scripts/slurm.pooled_cell
else
	echo 'ERROR: please indicate "single" or "pooled" run mode'
    exit
fi

configfile=/path/to/config/${SAMPLE_ID}.yaml
directory=/output_dir/${SAMPLE_ID}

# snakemake -n --reason -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} #dry run

snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --unlock #unlock a directory if there was a kill signal
snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} #start new runs or continue unfinished runs



