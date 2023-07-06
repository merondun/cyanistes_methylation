# Blue Tit 5mC Pilot

This project aims to investigate 5mC methylation varition among a pedigree of Blue Tits and assess reproducibility across batches and protocols. The directory structure of this project is as follows:

## Directory Structure

- `1.a.rrbs_processing`: Contains scripts and data related to Reduced Representation Bisulfite Sequencing (RRBS) and additionally EM-Seq with an MspI application.

- `1.b.wgbs_processing`: Contains scripts and data related to WGBS and additionally EM-Seq with WGS.

- `2.data_summaries`: Scripts for parsing bismark output files for plotting QC. 

- `3.cpg_extraction`: Contains scripts for extracting 5mC CpG methylation for each sample and overlapping them into a large matrix.

- `4.technical_replicates`: Scripts to assess reproducibility among samples.

- `bisulfite_conversion`: Barebones scripts for examining bisulfite conversion success.

- `conda_environments`: Contains the Conda environments.

- `0.a.sra_upload_organization`: Contains barebones scripts for parsing for upload to the Sequence Read Archive (SRA).

- `0.b.general_pre-processing`: Contains general utility scripts and files that don't belong to any specific part of the analysis but are used throughout the project.

- `metadata_n467_pilots.txt`: This text file contains metadata for the 467 pilot samples used in the project. This also includes some of the 2.datasummaries summary stats.  

- `pedigree`: This directory contains information about the pedigree of the Blue Tits.


## Contact Information

For further information or if you have any questions, please feel free to reach out: Justin merondun heritabilities@gmail.com

