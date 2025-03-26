#!/bin/bash
#SBATCH -p medium
#SBATCH -t 5-00:00
#SBATCH --mem=200
#SBATCH -c 1
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=EMAIL_ADDRESS

##usage: 
# sbatch -o log_dir/SAMPLE_ID.%j run_pipeline.sh SAMPLE_ID

module load gcc/6.2.0 conda2/4.2.13
source activate duplex-indel
export PATH=/path/to/duplex-indel:/path/to/seqtk:$PATH

snakefile=Snakefile.single_cell #CHANGE for Snakefile.pooled_cell
profile=/path/to/duplex-indel/slurm.single_cell #CHANGE for slurm.pooled_cell
configfile=/path/to/config.$1.yaml #$1 is SAMPLE_ID provided when running the sbatch command
directory=/output_dir/$1

# snakemake -n --reason -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} #dry run

snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --unlock #unlock a directory if there was a kill signal
snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} #start new runs or continue unfinished runs



