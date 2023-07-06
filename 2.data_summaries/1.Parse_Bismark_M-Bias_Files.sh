for file in *M-Bias.txt
do
    outfile="${file%.*}_processed.txt"
    context=""
    readtype=""
    prefix=$(echo ${file} | sed 's/\..*//g' )
    rm -f "$outfile"
    while IFS= read -r line
    do
        if [[ $line == *"context"* ]]; then
            context=$(echo $line | awk '{print $1}')
            if [[ $line == *"R1"* ]]; then
                readtype="R1"
            elif [[ $line == *"R2"* ]]; then
                readtype="R2"
            else
                readtype="SE"
            fi
        elif [[ $line == [1-9]* ]]; then
            echo -e "$line\t$context\t$readtype\t$prefix" >> "$outfile"
        fi
    done < "$file"
done


