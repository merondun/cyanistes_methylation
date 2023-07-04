#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=12
#SBATCH --time=200:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/Genomes/Cyanistes.caeruleus.v1
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0024/rawdata/Bisulfite-Seq/WGBS

basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing
trimdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing/scratch
workdir=${basedir}/scratch

#positional arguments
RUN=$1

mkdir ${workdir}
cd ${workdir}

F1=$( realpath ${rawdataWGBS}/${RUN}*__R1__*gz)
F2=$( realpath ${rawdataWGBS}/${RUN}*__R2__*gz)

#Trim adapters
trim_galore -j 8 --rrbs --quality 20 --output_dir ${trimdir} --basename ${RUN} --paired ${F1} ${F2}

#Map reads
bismark --parallel 9 --output_dir ${workdir} --genome ${genomedir} -1 ${trimdir}/${RUN}_val_1.fq.gz -2 ${trimdir}/${RUN}_val_2.fq.gz

#Deduplicate
deduplicate_bismark --paired ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam --output_dir ${workdir}

#Extract methylation
bismark_methylation_extractor --parallel 9 --gzip --bedgraph --buffer_size 60G --cytosine_report --genome_folder ${genomedir} --output ${workdir} ${workdir}/${RUN}_val_1_bismark_bt2_pe.deduplicated.bam
