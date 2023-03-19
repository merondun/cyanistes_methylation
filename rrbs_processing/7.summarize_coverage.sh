for RUN in $(cat PE_RRBS.list); do
zcat ${RUN}*.cov.gz | awk '{print $5+$6}' | datamash -H mean 1 median 1 sstdev 1 | awk -v s=${RUN} '{OFS="\t"}{print $0, s}' > ${RUN}.coverage
done
