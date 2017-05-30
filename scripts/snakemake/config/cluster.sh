#!/bin/bash -l

#$ -q short
#$ -cwd
#$ -V
#$ -l h_vmem=12G
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
reuse R-3.3

{exec_job}

