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

#import metadata
md = read.table('../metadata.txt',header=T)
s1 = subset(md,Batch=='Saarland_First') %>% select(ID)
s1 = intersect(names(d1),s1$ID)
s2 = subset(md,Batch=='Saarland_Second') %>% select(ID)
s2 = intersect(names(d1),s2$ID)

#divide samples by experiment
s1d = d1 %>% select(c(s1))
s2d = d1 %>% select(c(s2))
s1d$NAs = rowSums(is.na(s1d))
s1d = subset(s1d,NAs < 191)
s2d$NAs = rowSums(is.na(s2d))
s2d = subset(s2d,NAs < 192)
options(scipen=999)

missdat = rbind(s1d %>% select(NAs) %>% mutate(Experiment = 'First_Batch'),s2d %>% select(NAs) %>% mutate(Experiment = 'Second_Batch'))
misssum = missdat %>% group_by(Experiment) %>% dplyr::count()
missp = missdat %>% ggplot(aes(x=NAs,fill=Experiment))+
  geom_histogram(bins=50)+
  geom_text(data=misssum,aes(x=10,y=Inf,label=paste0(round(n/1000000,2),'M CpGs')),vjust=1)+
  scale_fill_manual(values=c('cadetblue','cadetblue2'))+
  xlab('Number of Samples with Missing Data')+
  facet_grid(Experiment~.)+
  ylab('Number of CpGs')+
  ggtitle(paste0('Missingness Within Each Batch, Calculated at each CpG Site'))+
  theme_classic()

png('../Missingness_Sample.png',height=7,width=10,res=300,units='in')
missp
dev.off()

#to compare, grab the 48 individuals with the least NAs in Saarland
missd = s1d %>%
  select(-(NAs)) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.)))) 
missd1 = pivot_longer(missd,cols = names(missd),names_to = 'ID')
missd2 = missd1 %>% arrange(value) %>% head(n=48)
sah = d1 %>% select(c(missd2$ID))
sah$NAs = rowSums(is.na(sah))
sah = subset(sah,NAs < 48)

missplot3 <- sah %>% ggplot(aes(x=as.factor(NAs)))+
  geom_histogram(fill='cadetblue',stat='count')+
  xlab('Number of Samples with Missing Data')+
  ylab('Number of CpGs')+
  ggtitle(paste0('Missingness, total sites: ',round(nrow(sah)/1000000,2),'M CpGs, \nFirst Batch HQ Individuals n=48'))+
  theme_classic()
missplot3

missd = s2d %>%
  select(-(NAs)) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.)))) 
missd1 = pivot_longer(missd,cols = names(missd),names_to = 'ID')
missd2 = missd1 %>% arrange(value) %>% head(n=48)
sah = d1 %>% select(c(missd2$ID))
sah$NAs = rowSums(is.na(sah))
sah = subset(sah,NAs < 48)

missplot4 <- sah %>% ggplot(aes(x=as.factor(NAs)))+
  geom_histogram(fill='cadetblue2',stat='count')+
  xlab('Number of Samples with Missing Data')+
  ylab('Number of CpGs')+
  ggtitle(paste0('Missingness, total sites: ',round(nrow(sah)/1000000,2),'M CpGs, \nSecond Batch HQ Individuals n=48'))+
  theme_classic()
missplot4

png('../Missingness_HQ-Saarland.png',height=6,width=8,res=300,units='in')
ggarrange(missplot3,missplot4,nrow=2)
dev.off()

#### Plot sample stats, discrete ####
md <- read.table('../metadata.txt',header=T)
covm <- merge(md,cov,by='ID')

covm$Discrete = factor(covm$Discrete,levels=c('Below10p','Below90p','Above90p','Below10x','Below25x','Below100x','Above100x'))
#plot coverage x sample 
cm1 <- covm %>% subset(Variable == 'Coverage' ) %>% 
  group_by(ID) %>% 
  mutate(count = sum(n)) %>% ungroup %>% 
  mutate(ord = fct_reorder(ID,count)) %>%
  ggplot(aes(x=ord,fill=Discrete,y=n))+
  geom_bar(stat='identity')+
  facet_grid(.~Batch,scales='free',space='free')+
  scale_fill_manual(values=viridis(5,option='inferno'))+
  theme_classic() + xlab('Sample')+
  ylab('Number of CpGs')+ggtitle('Total CpGs by Sample Classified by Coverage')+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

#plot % 5mC x sample 
cm2 <- covm %>% subset(Variable == 'Percent_5mC' ) %>% 
  group_by(ID) %>% 
  mutate(count = sum(n)) %>% ungroup %>% 
  mutate(ord = fct_reorder(ID,count)) %>%
  ggplot(aes(x=ord,fill=Discrete,y=n))+
  geom_bar(stat='identity')+
  facet_grid(.~Batch,scales='free',space='free')+
  scale_fill_manual(values=viridis(5,option='cividis'))+
  theme_classic() + xlab('Sample')+
  ylab('Number of CpGs')+ggtitle('Total CpGs by Sample Classified by % Methylation')+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

png('../Sample_Stats.png',height=8,width=10,units='in',res=300)
ggarrange(cm1,cm2,nrow = 2)
dev.off()

### #examine technical replicate correlations ####
library(GGally)
kept =  masterCG %>% dplyr::select(-c(site)) %>% names %>% data.frame
names(kept) <- 'ID'
md = read.table('../metadata.txt',header=T)
md1 = merge(md,kept,by='ID')
samp_count = md1 %>% dplyr::count(Biosample) %>% subset(n>1)
techs = samp_count$Biosample
rounding = c(1,2,3,4)
cordat=NULL
for (rd in rounding) { 
  cat('Estmating correlations if rounding to this level: ',rd,'\n')
  for (tech in techs) {
    cat('Working on technical replicate: ',tech,'\n')
    dd = masterCG %>% dplyr::select(contains(tech)) %>% na.omit
    dd = dd %>% mutate(across(everything(), round, rd))
    cd = cor(dd,method='spearman') %>% data.frame 
    cd$ID = row.names(cd)  
    co = cd %>% pivot_longer(!c(ID))
    names(co) = c('ID_A','ID_B','Correlation')
    cof = co %>% 
      rowwise() %>% 
      mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
      group_by(pair) %>%
      distinct(pair, .keep_all = T) %>% 
      separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair) %>% filter(ID_A != ID_B)
    cof$Sensitivity = rd
    cof$Biosample = tech
    cordat = rbind(cordat,cof)
  }
}
cs = cordat %>% group_by(Sensitivity) %>% summarize(mean=mean(Correlation))
cp = cordat %>% 
  ggplot(aes(x=Sensitivity,y=Correlation,fill=as.factor(Sensitivity)))+
  geom_violin(alpha=0.5,draw_quantiles = c(0.025,0.5,0.975))+
  scale_fill_manual(values=viridis(4))+
  geom_text(data=cs,aes(x=-Inf,y=Inf,label=paste0('rho: ',round(mean,2))),vjust=2,hjust=-.2)+
  geom_jitter()+
  facet_wrap(Sensitivity~.,scales='free')+
  theme_bw()+
  theme(axis.text.x=element_blank())

png('../Correlations_TechnicalReplicates.png',units='in',res=600,height=5,width=6)
cp
dev.off()

#### PCA on samples with not much missing data ####
library(ade4)
dk <- d1[rowSums(is.na(d1 %>% dplyr::select(contains('SRR')))) < 76, ]
ad = d1 %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.)))) 
ad1 = pivot_longer(ad,cols = names(ad),names_to = 'ID')
ad2 = ad1 %>% mutate(Missing = value/nrow(d1)) %>%  arrange(value)
ad2 = left_join(ad2,md)
#find which technical replicates both have low missing data
ad3 = ad2 %>% group_by(Biosample) %>% summarize(miss = mean(Missing)) %>% arrange(miss) %>% head(n=12)
adf = d1 %>% dplyr::select(contains(c(ad3$Biosample)))
adf$NAs = rowSums(is.na(adf))
adf = subset(adf,NAs < 1) %>% select(-c(NAs))

#perform ordination which allows for missing data 
#and traditional pca
adf$var <- apply(adf, 1, var, na.rm=TRUE)
pcv <- subset(adf,var > 0)
pcz <- pcv %>% select(-c(var)) %>% t(.)
pca <- prcomp(pcz,scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
vec <- data.frame(pca$x)
vec$ID <- row.names(vec)
vec <- merge(vec,md,by='ID')
samps = data.frame(Biosample = unique(vec$Biosample))
samps = samps %>% mutate(RepID = row_number())
vec = left_join(vec,samps)

np = vec %>% 
  ggplot(aes(x=PC1,y=PC2, shape = Batch,colour = Age,fill = Age,label=RepID)) +
  geom_point(size=6) +
  geom_text(col='white')+
  scale_color_manual(values=rev(viridis(3,option='turbo')))+
  scale_fill_manual(values=rev(viridis(3,option='turbo')))+
  scale_shape_manual(values=c(21,24))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

png('../Technical_PCA-12Reps.png',height=5,width=6,res=300,units='in')
np
dev.off()

#perform lPCA with missing data
library(missMDA)
lp = imputePCA(dk,method="EM",ncp=2)
pca <- prcomp(t(lp$completeObs),scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
vec <- data.frame(pca$x)
vec$ID <- row.names(vec)
vec <- merge(vec,md,by='ID')
samps = data.frame(Biosample = unique(vec$Biosample))
samps = samps %>% mutate(RepID = row_number())
vec = left_join(vec,samps)

np = vec %>% 
  ggplot(aes(x=PC1,y=PC2, shape = Batch,colour = Age,fill = Age,label=RepID)) +
  geom_point(size=6) +
  geom_text(col='white')+
  scale_color_manual(values=rev(viridis(3,option='turbo')))+
  scale_fill_manual(values=rev(viridis(3,option='turbo')))+
  scale_shape_manual(values=c(21,24))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

png('../Imputed_PCA-.png',height=5,width=6,res=300,units='in')
np
dev.off()


