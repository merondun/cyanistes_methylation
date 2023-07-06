#extract flowcell ids from fastq
for i in $(ls *fastq.gz | sed 's/.fastq.gz//g' ); do
zcat ${i}.fastq.gz | head | grep '@' | head -n 1 | tr ':' '\t' | awk -v s=${i} '{print s,$3}' > ${i}.FLOWCELL.txt
done
