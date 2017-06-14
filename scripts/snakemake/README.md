# *Snakemake workflow for ichorCNA*

## Description
This workflow will run the ichorCNA pipeline for a , starting from the BAM files and generating ichorCNA outputs. 

## Requirements
### Software packages or libraries
 - R-3.3
   - HMMcopy
   - optparse
 - Python 3.4 
   - snakemake-3.12.0
   - PyYAML-3.12
 - HMMcopy Suite (<http://compbio.bccrc.ca/software/hmmcopy/>).  
 		-In particular, `readCounter` is used.

### Scripts/executables
1. readCounter (C++ executable; HMMcopy Suite)
2. runIchorCNA.R

## cfDNA sample list
The list of cfDNA samples should be defined in a YAML file.  See `config/samples.yaml` for an example.  The field `samples` must to be provided.  
```
samples:
  tumor_sample_1:  /path/to/bam/tumor.bam
```

## snakefiles
1. `ichorCNA.snakefile`


Invoking the full snakemake workflow for ichorCNA
```
# show commands and workflow
snakemake -s ichorCNA.snakefile -np
# run the workflow locally using 5 cores
snakemake -s ichorCNA.snakefile --cores 5
# run the workflow on qsub using a maximum of 50 jobs. 
# Broad UGER cluster parameters can be set directly in config/cluster.sh. 
snakemake -s ichorCNA.snakefile --cluster-sync "qsub" -j 50 --jobscript config/cluster.sh
```
