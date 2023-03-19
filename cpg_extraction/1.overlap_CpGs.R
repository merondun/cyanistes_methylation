library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(forcats)
library(matrixStats)
library(ggpubr)
library(viridis)

#takes bismark .cov.gz files and identifies overlaps, creates a massive text file with %5mC, modify to get count data 
setwd('/dss/dssfs03/pn69za/pn69za-dss-0001/di39dux/files/tits/alldata')
#grab all files
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)
allnames <- NULL
cov <- NULL

#loop through them..
for (samp in files) {

  cat('Reading in file: ',samp,'\n')
  #read data
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);

  #extract name which will be our new variable
  name <- gsub('./','',samp)
  name <- gsub('_val_1_bismark_bt2_pe.bismark.cov.gz','',name)
  name <- gsub('_trimmed_bismark_bt2.bismark.cov.gz','',name)

  #save coverage and %5mc statistics per sample
  c <- tab %>% mutate(coverage=V5+V6,ID=name) %>% select(coverage,V4,ID)
  names(c) <- c('Coverage','Percent_5mC','ID')
  c1 <- c %>% mutate(Discrete = cut(Coverage, breaks = c(0,10,25,100,Inf), right = F, labels = c('Below10x','Below25x','Below100x','Above100x'))) %>% group_by(Discrete,ID) %>% dplyr::count()
  c2 <- c %>% mutate(Discrete = cut(Percent_5mC, breaks = c(-1,10,90,101), right = F, labels = c('Below10p','Below90p','Above90p'))) %>% group_by(Discrete,ID) %>% dplyr::count()
  c1$Variable <- 'Coverage'
  c2$Variable <- 'Percent_5mC'
  cc <- rbind(c1,c2)
  cov <- rbind(cov,cc)

  #only keep positions with 10x coverage or higher
  tab <- subset(tab, (V5 + V6) > 9)
  tab$site <- paste0(tab$V1,'_',tab$V2)
  tab$V4 <- tab$V4/100

  #grab only site / countM / countU / %5mC
  #tab2 <- tab[,c(7,5,6,4)]
  #names(tab2) <- c('site',paste0('Cs_',name),paste0('Ts_',name),name);
  #or for quick, just grab site and %5mC
  tab2 <- tab[,c(7,4)]
  names(tab2) <- c('site',name);

  cat('Saved to variable: ',name,'\n')
  #assign it to a new variable
  assign(paste0(name),tab2)

  #add this sample to our total list
  allnames <- unique(c(allnames,name))
}

#get names
vars <- mget(allnames)
# merge the points, keeping all
masterCG <- Reduce(function(x,y) merge(x = x, y = y, by=c('site'),all=TRUE),vars)

#filter according to missingness, if there's 10% missing data (n=20) more more missing individuals, drop the site
covdat <- masterCG[rowSums(is.na(masterCG[grepl('^B', names(masterCG))])) <= 10, ]

write.table(masterCG,file='5mC_Data_Pilot_10x.txt',quote=F,sep='\t',row.names=FALSE)
write.table(cov,file='Coverage-Sample_Statistics.txt',quote=F,sep='\t',row.names=F)

#write file with missingness filters
write.table(covdat,file='5mC_Data_Pilot_10x-MM1.txt',quote=F,sep='\t',row.names=FALSE)
