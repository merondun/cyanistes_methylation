#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=6
#SBATCH --time=200:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/Genomes/Cyanistes.caeruleus.v1
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0024/rawdata/Bisulfite-Seq/RRBS

basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing
trimdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing/scratch
workdir=${basedir}/scratch

#positional arguments
RUN=$1

mkdir ${workdir}

cd ${workdir}

#Trim adapters
trim_galore -j 6 --rrbs --quality 20 --output_dir ${trimdir} --basename ${RUN} $rawdataRRBS/${RUN}*__SE__*fastq.gz

#Map reads
bismark --parallel 6 --output_dir ${workdir} --genome ${genomedir} ${trimdir}/${RUN}_trimmed.fq.gz
mv ${workdir}/${RUN}_trimmed_bismark_bt2.bam ${workdir}/${RUN}_val_1_bismark_bt2_pe.bam
#Extract methylation
bismark_methylation_extractor --parallel 6 --gzip --bedgraph --buffer_size 60G --cytosine_report --genome_folder ${genomedir} --output ${workdir} ${workdir}/${RUN}_val_1_bismark_bt2_pe.bam
