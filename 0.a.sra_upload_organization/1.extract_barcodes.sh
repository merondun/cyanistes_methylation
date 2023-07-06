#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

#have seqkit available, submit with fastq as positional argument
i=$1
zcat ${i}.fastq.gz | grep '@' | sed 's/.*://g' | sort | uniq -c > ${i}.BARCODES.txt
sort -nk1 ${i}.BARCODES.txt | tail -n1 | awk '{print $2}' | tr '+' '\n' | awk '$0=">ID"NR"\n"$0' > ${i}.BARCODE.txt
seqkit seq --seq-type DNA -v --complement --reverse ${i}.BARCODE.txt | grep -v '>' | tr '\n' '\t' | awk -v s=${i} '{OFS="\t"}{print s,$1, $2}'> ${i}.Lookup.txt
