#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

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
trim_galore -j 4 --rrbs --quality 20 --output_dir ${trimdir} --basename ${RUN} --paired $rawdataRRBS/${RUN}*__R1__*fastq.gz $rawdataRRBS/${RUN}*__R2__*fastq.gz

#Map reads
bismark --parallel 4 --output_dir ${workdir} --genome ${genomedir} -1 ${trimdir}/${RUN}_val_1.fq.gz -2 ${trimdir}/${RUN}_val_2.fq.gz

#Extract methylation
bismark_methylation_extractor --parallel 4 --gzip --bedgraph --buffer_size 14G --cytosine_report --genome_folder ${genomedir} --output ${workdir} ${workdir}/${RUN}_val_1_bismark_bt2_pe.bam

