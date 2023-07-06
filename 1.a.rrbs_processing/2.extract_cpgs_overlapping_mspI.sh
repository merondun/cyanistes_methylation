#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

cuts=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/Genomes/Cyanistes.caeruleus.v1/CutsitesCyanistesMspI100_400.bed
gfile=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/Genomes/Cyanistes.caeruleus.v1/Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.genome

for RUN in $(cat PE_RRBS.list); do 
zcat ${RUN}*cov.gz | bedtools sort -i - | bedtools closest -g $gfile -d -a - -b $cuts | awk -v s=${RUN} '{OFS="\t"}{print s, $1, $10}' > ${RUN}.cuts.tmp;
done
