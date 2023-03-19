#within directory which has all the alignment and mbias files
for RUN in $(cat PE_RRBS.list ); do 

sed -n '/CHG context (R1)/q;p' ${RUN}*M-bias.txt | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R1","CpG"}' > ${RUN}.MBIAS1.tmp

sed -n '/CpG context (R2)/,$p' ${RUN}*M-bias.txt | sed -n '/CHG context (R2)/q;p'  | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R2","CpG"}' > ${RUN}.MBIAS2.tmp

cat ${RUN}.MBIAS1.tmp ${RUN}.MBIAS2.tmp > ${RUN}.M-Bias.txt
rm ${RUN}*MBI*tmp
done

#merge together
*M-Bias.txt > M_Bias.txt

#messy script to extract counts data from alignment reports
egrep 'efficiency|CpG' *E_report.txt | sed 's/_val_1_bismark_bt2_PE_report.txt:/\t/g' | sed 's/_trimmed_bismark_bt2_SE_report.txt:/\t/g' | sed "s/Total methylated C's in CpG context:/Methylated_Cs/g" | sed "s/Total unmethylated C's in CpG context:/Unmethylated_Cs/g" | sed "s/C methylated in CpG context:/Global_Methylation/g" | sed "s/Mapping efficiency:/Alignment_Rate/g" > Counts_Summary.txt

#and extract raw sequencing effort
egrep 'Sequence pairs analysed|unique best hit|s analysed' *E_report.txt | sed 's/_val_.*txt:/\t/g' | sed "s/'//g" | sed 's/-/_/g' | sed 's/ /_/g' | sed 's/://g' > Effort_Summary.txt

