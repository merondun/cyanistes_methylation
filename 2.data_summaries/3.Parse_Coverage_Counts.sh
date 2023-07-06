#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

rm Coverage_Statistics.txt

for i in $(ls *cov.gz | sed 's/.cov.gz//g'); do

abv0=$(zcat ${i}.cov.gz | wc -l)
abv10=$(zcat ${i}.cov.gz | awk '$5+$6>=10' | wc -l)
echo -e "${i}\t${abv0}\t${abv10}" >> Coverage_Statistics.txt

done
