# RRBS and EM-seq Data Processing Scripts

This repository contains a series of scripts for the processing of RRBS and EM-seq data, specifically aligned to the blue tit reference genome.

1. `1.align_to_cyanistes_reference.sh` - Aligns the PE RRBS reads to the Cyanistes reference genome provided a prefix with sample ID and SRR accession.
2. `1.b_align_to_cyanistes_reference_single_end.sh` - Aligns the single-end RRBS reads to the Cyanistes reference genome (e.g. Saarland SE libs).
3. `1.c_align_to_cyanistes_reference_EMseq-PE.sh` - Aligns paired-end EM-seq reads to the Cyanistes reference genome (e.g. Kiel PE libs). I specify R1 and R2 with a sample prefix explicitly as no accessions exist atm.
4. `2.extract_cpgs_overlapping_mspI.sh` - Extracts 5mC CpG sites that overlap with MspI restriction sites.
5. `3.summarize_cutsites.R` - Summarizes the 5mC CpG overlap with the cut sites.
6. `4.extract_cpgs_overlapping_CpGislands.sh` - Extracts CpG sites that overlap with annotated CpG islands.
7. `5.summarize_cgis.R` - Summarizes the 5mC CpGs that overlap CpG islands.



