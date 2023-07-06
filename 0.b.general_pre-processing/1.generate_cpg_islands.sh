#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

#make sure you have access to the R script with makecgi installed
#have a folder rawdata/${genome}.fa 
#confirmed with version makeCGI_1.3.4.tar.gz

SCRATCH=/tmp/$SLURM_JOB_ID

Rscript makeCGI.R

#afterwards filter at least 200bp, GC > 50%, obs/expected CG ratio > 60% 
awk '($4 > 200) && ($7 > 0.5) && ($8 > 0.6){print $1, $2, $3}' CGI-Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.CHR.txt | \
	sed '1d' | \
	tr ' ' \\t  > ../Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.CGI.bed

