# Basic Summary Statistics

This repository contains a series of scripts for summarizing M-bias and counts, global methylation rates, and alignment rates. Mainly parsing bismark outputs. 

1. `1.Parse_Bismark_M-Bias_Files.sh` - Converts the bismark output file to a tidy format, works with both PE and SE outputs. 
2. `2.Parse_Bismark_Alignment_Reports.sh` - Parses the alignment reports and outputs a tidy format with lots of details 
3. `3.Parse_Coverage_Counts.sh` - Parses the cov.gz outputs with $ID $RawCpGs $CpGs>10X
4. `4.plot_summaries.R` - Old script for plotting summaries, could use some parts likely.  
5. `M_Bias_Input-n467.txt.gz` - Output of (1) with M-bias information with n=467 individuals. 
6. `Alignment_Counts_Summary-n467.txt` - Output of (2) with alignment information with n=467 individuals. 
7. `Coverage_Summary-n467.txt` - Output of (3) with coverage information with n=467 individuals. 

