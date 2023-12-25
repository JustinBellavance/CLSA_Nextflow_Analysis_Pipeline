#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=35:58:58
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

module load nextflow
nextflow run main.nf -resume
#nextflow run test.nf

#bash run_pheweb.sh