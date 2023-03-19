setwd('F:/Research/scratch/epi_tit/pilot_kiel/')

library(ggplot2)
library(gridExtra)
library(ggsci)
library(viridis)
library(dplyr)
library(ggpubr)

#mbias 
md <- read.table('metadata.txt',header=T)
bias <- read.table('M_Bias.txt',header=F)
names(bias) <- c('ID','Positions','Methylated','Nonmethylated','Percent_Methylated','Coverage','Read','Context')
bb <- merge(bias,md,by=Reduce(intersect, list(names(bias),names(md))))

bp <- ggplot(bb,aes(x=Positions,col=ID))+
  geom_line(aes(y=Percent_Methylated),stat="identity",show.legend=F,size=1)+ylab('Percent Methylation')+
  theme_minimal(base_size=16)+facet_grid(Read~Experiment)+
  scale_color_viridis(discrete=T,option='turbo')+
  coord_cartesian(ylim=c(25,75))+
  ggtitle('Methylation Along Read Length')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
bp

png('MBias_Statistics.png',height=6,width=10,res=300,units='in')
bp
dev.off()

#counts / global
md <- read.table('metadata.txt',header=T)
cc <- read.table('Counts_Summary.txt',header=F)
names(cc) <- c('ID','Variable','Value')
c <- merge(cc,md,by=Reduce(intersect, list(names(bias),names(md))))
c <- c %>% mutate(Value = gsub('%','',Value)) %>% mutate_at(grep("Value",colnames(.)),funs(as.numeric))
c <- c %>% mutate(Value = ifelse(Variable == 'Methylated_Cs',Value/1000000,
                                 ifelse(Variable == 'Unmethylated_Cs',Value/1000000,Value)))
options(scipen=999)

c1 <- ggplot(c,aes(x=ID,col=Experiment,y=Value,shape=Age))+
  facet_wrap(Variable~.,scales='free')+
  geom_point(size=3)+
  theme_minimal(base_size=16)+
  scale_color_manual(values=viridis(5,option='turbo'))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

c1

png('Summary_Statistics.png',height=6,width=10,res=300,units='in')
c1
dev.off()

#Cutsites
library(data.table)
ct = fread('Distance-Cutsites.txt.gz')
ct1 = merge(ct,md %>% select(ID,Experiment),by='ID')
ct1 = subset(ct1,Experiment == 'Saarland')
cr = fread('Distance-Cutsites_CROW.txt.gz')
kr = fread('Distance-Cutsites_Kiel.txt.gz')
kr$Experiment <- 'IKMB'

ct3 = rbind(ct1,cr,kr)
ct3 = subset(ct3,distance != -1)
#alternatively, visualize bases within 100bp, or bases not within 100bp
ct2 = ct3
ct2 = ct2 %>% mutate_at(vars(distance,count),as.integer)
ct2 = ct2 %>% mutate(Discrete = ifelse(distance == 0,'WithinFragment',
                                       ifelse(distance < 100, 'Within100bp','Outside100bp'))) %>% na.omit
ct2 = ct2 %>% group_by(ID,Discrete,Experiment) %>% summarize(Sites = sum(count))
ct2 = ct2 %>% mutate(XP = ifelse(Experiment == 'CG','Crow_RRBS1',
                                 ifelse(Experiment=='HZ','Crow_RRBS2',
                                        ifelse(Experiment =='WGBS','Crow_WGBS',Experiment))))
ct2$XP = factor(ct2$XP,levels=c('Crow_WGBS','Crow_RRBS1','Crow_RRBS2','IKMB','Kiel','Saarland'))
ctp = ct2 %>% 
  group_by(ID,XP) %>%
  mutate(value = Sites/sum(Sites) * 100) %>%
  ggplot(aes(x=ID,y=value,fill=Discrete))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=c('black','yellow','seagreen'))+
  facet_grid(.~XP,scales='free',space='free')+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
ctp

#as total
ctp2 = ct2 %>% 
  ggplot(aes(x=ID,y=Sites,fill=Discrete))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=c('black','yellow','seagreen'))+
  facet_grid(.~XP,scales='free',space='free')+
  theme_classic()+
  geom_hline(yintercept=1000000,lty=2)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
ctp2

png('Cut-Sites.png',height=9,width=18,units='in',res=300)
ggarrange(ctp,ctp2,nrow=2)
dev.off()

#CpG islands
library(data.table)
cg = fread('CGI-DistanceHists.txt.gz')
cg1 = merge(cg,md %>% select(ID,Experiment),by='ID')
cr = fread('CGI-DistanceHists-Crow.txt.gz')
names(cr) <- c('ID','distance','count','Experiment')
cg3 = rbind(cg1,cr)
cg3 = subset(cg3,distance != -1)
#alternatively, visualize bases within 100bp, or bases not within 100bp
cg2 = cg3
cg2 = cg2 %>% mutate_at(vars(distance,count),as.integer)
cg2 = cg2 %>% mutate(Discrete = ifelse(distance == 0,'WithinIsland',
                                       ifelse(distance < 100, 'Within100bp','Outside100bp'))) %>% na.omit
cg2 = cg2 %>% group_by(ID,Discrete,Experiment) %>% summarize(Sites = sum(count))
cg2 = cg2 %>% mutate(XP = ifelse(Experiment == 'CG','Crow_RRBS1',
                                 ifelse(Experiment=='HZ','Crow_RRBS2',
                                        ifelse(Experiment =='WGBS','Crow_WGBS',Experiment))))
cg2$XP = factor(cg2$XP,levels=c('Crow_WGBS','Crow_RRBS1','Crow_RRBS2','IKMB','Kiel','Saarland'))

cgp = cg2 %>% 
  group_by(ID,XP) %>%
  mutate(value = Sites/sum(Sites) * 100) %>%
  ggplot(aes(x=ID,y=value,fill=Discrete))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=viridis(3))+
  facet_grid(.~XP,scales='free',space='free')+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
cgp

#as total
cgp2 = cg2 %>% 
  ggplot(aes(x=ID,y=Sites,fill=Discrete))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=viridis(3))+
  facet_grid(.~XP,scales='free',space='free')+
  theme_classic()+
  geom_hline(yintercept=1000000,lty=2)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
cgp2


png('CGI-Sites.png',height=9,width=18,units='in',res=300)
ggarrange(cgp,cgp2,nrow=2)
dev.off()

### coverage
library(tidyverse)
cc <- read.table('Coverage.txt',header=F)
names(cc) <- c('mean','median','sd','ID')
cc$ID <- gsub('.coverage','',cc$ID)
ccm <- merge(cc,md,by='ID') %>% select(ID,mean,median,Experiment)
ccmp =ccm %>% 
  pivot_longer(cols=c(mean,median)) %>%
  ggplot(aes(x=ID,y=value,col=Experiment))+
  facet_grid(name~Experiment,scales='free',space='free')+
  geom_point()+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
png('Coverage.png',height=6,width=8,res=300,units='in')                    
ccmp
dev.off()
