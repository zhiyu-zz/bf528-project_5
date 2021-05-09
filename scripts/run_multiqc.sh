#!/bin/bash -l
#$ -P bf528
#$ -cwd
#$ -o /projectnb/bf528/users/dreadlocks/project_5/sarah/run_mqc_o
#$ -e /projectnb/bf528/users/dreadlocks/project_5/sarah/run_mqc_e

module load python3/3.7.9
module load multiqc/1.10.1

PATH1=/projectnb2/bf528/users/dreadlocks/project_5/sarah/multiqc/fastqc/
PATH2=/projectnb2/bf528/users/dreadlocks/project_5/sarah/multiqc/star/
PATH=/projectnb2/bf528/users/dreadlocks/project_5/sarah
multiqc -o $PATH1 $PATH/fastqc/
multiqc -o $PATH2 $PATH/STAR/
