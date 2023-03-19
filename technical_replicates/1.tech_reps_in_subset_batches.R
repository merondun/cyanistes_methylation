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

masterCG = fread('5mC_Data_Pilot_10x.txt',header=T)
cov = read.table('Coverage-Sample_Statistics.txt',header=T)
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
#Here, try subsetting 6 reps at a time 
for (grp in seq(1,length(techs),6)) {
  cat('Working on technical replicate: ',tech,'\n')
  tech = techs[grp:(grp+5)]
  dd = masterCG %>% dplyr::select(contains(c(tech))) %>% na.omit #grab only those individuals, each column is an individual
  cd = cor(dd,method='spearman') %>% data.frame 
  cd$ID = row.names(cd)  
  co = cd %>% pivot_longer(!c(ID))
  names(co) = c('ID_A','ID_B','Correlation')
  #remove redundant comparisons (eg ID_A v ID_B and ID_B v ID_A as well as ID_A = ID_B) 
  cof = co %>% 
    rowwise() %>% 
    mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
    group_by(pair) %>%
    distinct(pair, .keep_all = T) %>% 
    separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair) %>% filter(ID_A != ID_B)
  #add biosample
  cof = cof %>% mutate(Biosample_A = gsub('__.*','',ID_A),
                       Biosample_B = gsub('__.*','',ID_B))
  cof = cof %>% mutate(Technical = ifelse(Biosample_A == Biosample_B,'Technical','Inter-individual'))
  cordat = rbind(cordat,cof)
}
cs = cordat %>% group_by(Technical) %>% summarize(mean=mean(Correlation),sd=sd(Correlation))
cp = cordat %>% 
  ggplot(aes(x=Technical,y=Correlation,fill=Technical))+
  geom_violin(alpha=0.5,draw_quantiles = c(0.025,0.5,0.975))+
  scale_fill_manual(values=c('darkblue','darkorchid2'))+
  geom_text(data=cs,aes(x=Technical,y=Inf,label=paste0('rho: ',round(mean,2), '+/- ',round(sd,2))),vjust=2)+
  geom_jitter(alpha=0.2)+
  theme_bw()+
  theme(axis.text.x=element_blank())

png('../Correlations_TechnicalReplicates_N6.png',units='in',res=600,height=5,width=6)
cp
dev.off()

