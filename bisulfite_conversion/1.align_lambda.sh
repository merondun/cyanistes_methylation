#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/lambda
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/lambda
trimdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing/scratch
outdir=${basedir}/scratch

#positional arguments, add fastq prefix from trimmed reads
RUN=$1

mkdir ${outdir}

cd ${outdir}

#Map reads
bismark --parallel 2 --output_dir ${outdir} --genome ${genomedir} -1 ${trimdir}/${RUN}_val_1.fq.gz -2 ${trimdir}/${RUN}_val_2.fq.gz
