#!/usr/bin/env Rscript
#ensure you have your file formatted as families, all samples must be present
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/pedigree')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(kinship2)
library(dplyr)
ped = read.table('families.txt',header=TRUE)
ped = ped %>% mutate(can_use = ifelse(good_data == 1 & count > 0 & Greater500 == 1,1,0),
                     RE_ID = ifelse(Reextracted == 1,ID,''),
                     UNRE_ID = ifelse(Reextracted == 0,ID,''))
pedAll <- pedigree(id=ped$id, 
                   dadid=ped$father, momid=ped$mother, 
                   sex=ped$sex, famid=ped$ped,
                   affected=ped$can_use)
ped1basic <- pedAll['1']
pdf('PedigreeAll.pdf',height=8,width=30)
plot(ped1basic, col=ped$count,id=ped$UNRE_ID,cex=0.5)
dev.off()

