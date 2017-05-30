# *Snakemake workflow for TITAN*

## Description
This workflow will run the TITAN a set of tumour-normal pairs, starting from the BAM files and generating TitanCNA outputs.  It will also perform model selection at the end of the workflow to choose the optimal ploidy and clonal cluster solutions.

## Requirements
### Software packages or libraries
 - R-3.3
   - TitanCNA (v1.15.0)
   		- TitanCNA imports: GenomicRanges, dplyr, data.table, doMC
   - HMMcopy
   - optparse
   - stringr
   - SNPchip
 - Python 3.4 
   - snakemake-3.12.0
   - PySAM-0.11.2.1
   - PyYAML-3.12
 - samtools v1.3.1
 - bcftools v1.1
 - HMMcopy Suite (<http://compbio.bccrc.ca/software/hmmcopy/>).  
 		-In particular, `readCounter` is used.

### Scripts/executables
 - readCounter (C++ executable; HMMcopy Suite)
 - ulp_wgs_hmm.R (ichorCNA tool for normalizing and correcting read coverage)
 - countPysam.py (generates input allele counts)
 - titanCNA.R (main R script to run TitanCNA)
 - selectSolution.R (R script to select optimal solution for each sample)

## Tumour-Normal sample list
The list of tumour-normal paired samples should be defined in a YAML file.  See `config/samples.yaml` for an example.  Both fields `samples` and `pairings` must to be provided.  `pairings` key must match the tumour sample while the value must match the normal sample.
```
samples:
  01115_466.TM.EX:  /xchip/bloodbiopsy/gavinha/cfDNA_WES/PRAD/samples/01115_466.TM.EX.bam
  01115_466.BN.EX:  /xchip/bloodbiopsy/gavinha/cfDNA_WES/PRAD/samples/01115_466.BN.EX.bam


pairings:
  01115_466.TM.EX:  01115_466.BN.EX
```


## snakefiles
1. `ichorCNA.snakefile`
2. `getAlleleCounts.snakefile`
3. `TitanCNA.snakefile`

Invoking the full snakemake workflow for TITAN
```
# show commands and workflow
snakemake -s TitanCNA.snakefile -np
# run the workflow locally using 5 cores
snakemake -s TitanCNA.snakefile --cores 5
# run the workflow on qsub using a maximum of 50 jobs. Broad UGER cluster parameters can be set directly in config/cluster.sh. 
snakemake -s TitanCNA.snakefile --cluster-sync "qsub" -j 50 --jobscript config/cluster.sh
```
This will also run both `ichorCNA.snakefile` and `getAlleleCounts.snakefile` which generate the necessary inputs for `TitanCNA.snakfile`.

`ichorCNA.snakefile` and `getAlleleCounts.snakefile` can also be invoked separately. If only one but not both results are needed, then you can invoke the snakefiles independently.
```
snakemake -s ichorCNA.snakefile --cores 5
# OR
snakemake -s getAlleleCounts.snakefile --cores 5
``` 





