# An analysis pipeline for CLSA-PheWeb project (December 2023).
## Includes the analysis of quantitative traits into GWAS catalog pheweb via nextflow, including a sex-stratified analysis. All done via slurm scheduler (more specificially, the Alliance's Graham cluster).

## Files present
- main.nf -> main nextflow pipeline script, includes all processes and workflow.
- nextflow.config -> config file used for main analysis
- run_pipeline.sh -> bash script used to launch pipeline

## After pipeline is done, user has to run:
 $ pheweb process
 $ pheweb serve --open

## to properly see their pheweb results! 