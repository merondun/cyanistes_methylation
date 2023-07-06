# Takes bismark alignment_report.txt as input. The $ID will be the string preceding '.alignment_report.txt', so rename your files accordingly
echo -e "ID\ttotal_sequences\tunique_hits\tmapping_rate\ttotal_cytosines\tcpg_rate\tchg_rate\tchh_rate" > processed_reports.txt

for file in *.alignment_report.txt
do
    ID=$(basename "$file" .alignment_report.txt)
    total_sequences=$(grep "analysed in total:" "$file" | awk '{print $NF}')
    unique_hits=$(grep "unique best hit" "$file" | awk '{print $NF}') # 
    mapping_rate=$(grep "Mapping efficiency:" "$file" | awk '{print $3}')
    total_cytosines=$(grep "Total number of C's analysed:" "$file" | awk '{print $6}')
    cpg_rate=$(grep "C methylated in CpG context:" "$file" | awk '{print $6}')
    chg_rate=$(grep "C methylated in CHG context:" "$file" | awk '{print $6}')
    chh_rate=$(grep "C methylated in CHH context:" "$file" | awk '{print $6}')

    # Append to the output file
    echo -e "$ID\t$total_sequences\t$unique_hits\t$mapping_rate\t$total_cytosines\t$cpg_rate\t$chg_rate\t$chh_rate" >> processed_reports.txt
done
