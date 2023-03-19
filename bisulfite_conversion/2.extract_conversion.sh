#simply extract the methylated in CpG context from the alignment reports
egrep 'methylated in CpG' *txt |sed 's/_val.*context://g' > Lambda.txt
