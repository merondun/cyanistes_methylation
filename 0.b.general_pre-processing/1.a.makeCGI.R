#generate CpG island track with makecgi, you can modify '5000' to your desired window size (largely insensitive). ensure you have a 
library(makeCGI);
.CGIoptions=CGIoptions();
.CGIoptions$rawdat.type="txt";
.CGIoptions$species="Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1";
.CGIoptions;
makeCGI(.CGIoptions);
5000
