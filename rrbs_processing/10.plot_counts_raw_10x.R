#will need to install merondun/meRo R package from github
setwd('/dss/dsshome1/lxc07/di39dux/merondun/cyanistes_methylation')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(meRo)

d = read.table('metadata_n467_pilots.txt',header=TRUE,sep='\t',comment.char='') %>% as_tibble
d = d %>% mutate_at(c('Number10XCpGs','NumberRawCpGs','MappingRate'),as.numeric) %>% drop_na(Number10XCpGs) 
dp = d %>% 
  pivot_longer(c(MappingRate,NumberRawCpGs,Number10XCpGs)) %>%
  group_by(Batch,name) %>% 
  sum_stats(value) %>% 
  ggplot(aes(x=Batch,y=mean,ymax=mean+sd,ymin=mean-sd,col=Batch))+
  geom_point()+
  ylab('Mean +/- SD')+
  scale_color_viridis(discrete=TRUE)+
  facet_grid(name~.,scales='free')+
  geom_errorbar()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

png('Summary_Statistics.png',height=5,width=9,units='in',res=300)
dp
dev.off()

