setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/lambda/scratch')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(tidyverse)
library(forcats)
library(ggpubr)
library(viridis)

d = read.table('Lambda.txt',header=FALSE)
names(d) = c('ID','Methylation_Rate')
bs = d %>% mutate(Methylation_Rate = gsub('%','',Methylation_Rate)) %>% mutate_at(vars(Methylation_Rate),as.numeric) %>%
  mutate(Bisulfite_Conversion_Rate = 100-Methylation_Rate)

#import metadata
md = read.table('../../processing/metadata.txt',header=T)
bs = left_join(bs,md)

bs %>% group_by(Batch) %>% summarize(median = median(Bisulfite_Conversion_Rate))

#plot
bsp = bs %>% 
  ggplot(aes(x=Batch,y=Bisulfite_Conversion_Rate,fill=Batch))+
  geom_violin(alpha=0.5,draw_quantiles = c(0.025,0.5,0.975))+
  scale_fill_manual(values=c('cadetblue','cadetblue2'))+
  geom_jitter()+
  theme_bw()

png('Bisulfite_Conversion_Rate.png',height=5,width=6,units='in',res=600)
bsp
dev.off()

