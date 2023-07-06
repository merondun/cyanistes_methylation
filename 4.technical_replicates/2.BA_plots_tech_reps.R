setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing/scratch')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(forcats)
library(matrixStats)
library(ggpubr)
library(viridis)
library(data.table)
library(gmodels)

masterCG = fread('5mC_Data_Pilot_10x.txt',header=T)
#plot missingness
d1 = masterCG %>% dplyr::select(!(starts_with(c('site'))))

### #examine technical replicate correlations ####
kept =  masterCG %>% dplyr::select(-c(site)) %>% names %>% data.frame
names(kept) <- 'ID'
md = read.table('../metadata.txt',header=T)
md1 = merge(md,kept,by='ID')
samp_count = md1 %>% dplyr::count(Biosample) %>% subset(n>1)
techs = samp_count$Biosample
cordat=NULL
trun = techs[1:12]
#Here, try subsetting each rep one at a time 
counter = 0
for (tech in trun) {
  cat('Working on technical replicate: ',tech,'\n')
  dd = masterCG %>% dplyr::select(site,contains(c(tech))) %>% na.omit #grab only those individuals, each column is an individual
  dd = dd %>% separate(site,into=c('chr'))
  #remove reads on scaffolds...only on chromosomes then 
  dd = dd %>% filter(chr!= 'NW') %>% select(-chr)
  dd = dd %>% slice_sample(n=2500)
  names(dd) = c('T1','T2')
  dd$avg = rowMeans(dd)
  dd = dd %>% mutate(diff = T1 - T2,
                     summed = T1 + T2)
  
  mean_diff <- mean(dd$diff)
  mean_diff
  lower <- mean_diff - 1.96*sd(dd$diff)
  upper <- mean_diff + 1.96*sd(dd$diff)

  #create Bland-Altman plot
  gp = ggplot(dd, aes(x = avg, y = diff)) +
    geom_hex(bins=25)+
    scale_fill_continuous(low='yellow',high='red')+
    geom_hline(yintercept = mean_diff) +
    geom_hline(yintercept = lower, color = "red", linetype="dashed") +
    geom_hline(yintercept = upper, color = "red", linetype="dashed") +
    ggtitle(tech) +
    ylab("") +
    xlab("")

  counter = counter + 1
  assign(paste0('p',counter),gp)
}
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,common.legend = TRUE)

