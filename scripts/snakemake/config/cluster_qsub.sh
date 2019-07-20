#!/bin/bash -l

#$ -q broad
#$ -cwd
#$ -V
#$ -l h_vmem=4G,h_rt=6:00:00
#$ -pe smp 1
#$ -binding linear:1
#$ -o logs/cluster/
# join stdout and stderr output
#$ -j y
#$ -sync y
#$ -b y
#$ -r y

source /broad/software/scripts/useuse
reuse Bcftools
reuse Samtools
reuse Python-3.4
reuse R-3.4

{exec_job}

