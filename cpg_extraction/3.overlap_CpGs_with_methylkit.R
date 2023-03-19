#this takes bismark cov.gz files and imports into methylkit. Here the 'treatment' vector is simply biosample since this was used to analyze technical replicates
#if you want to analyze the entire all files together, unhash the initial import samples lines
#the below script was used to analyze technical replicates while also including 3 external random samples (1 with 1 cov, 1 with medium cov, 1 with low cov). 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023FEB/processing/scratch')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(tidyverse)
.libPaths('~/miniconda3/envs/methylation/lib/R/library')
library(methylKit)

#find your files
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)

#extract sample lists, we have 4 treatments : age, sex, tissue, taxon ( NO TISSUE IN RRBS ) 
name <- gsub('./','',files)
name <- gsub('_val_1_bismark_bt2_pe.bismark.cov.gz','',name) %>% as.data.frame()
names(name) <- c('ID')
name <- name %>% separate(ID,remove=F,into=c('Ring','Tissue','Age','Sex','Accession'))
name = name %>% mutate(Biosample = paste0(Ring,'_',Tissue,'_',Age,'_',Sex))
#give each sample an original file order number to re-match it with afterwrads
name <- name %>% mutate(Order = row_number())
#extract treatments based on unique taxon/tissue/age/sex
vars <- name %>% dplyr::select(Biosample) %>% unique
nrow(vars)
#rebind them with the order number, so we now have ID / treatment vectors
vars <- vars %>% mutate(Treatment = row_number())
name <- merge(name,vars,by=Reduce(intersect, list(names(name),names(vars))))
name$Treatment <- as.numeric(name$Treatment)
name <- name %>% arrange(Order)

# #Import samples, start with a low min coverage so they aren't filtered out 
# CpG.bismark.covM=methRead(location=as.list(files), 
#                           sample.id=as.list(name$ID),
#                           treatment=name$Treatment,
#                           assembly="ccornix5.7",context="CpG",pipeline="bismarkCoverage",
#                           mincov = 10)
# 
# covm=filterByCoverage(CpG.bismark.covM,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99)
# covm
# 
# #Require at least 1 representative of each treatment (TAXON / SEX / AGE, 12 treatments)
# meth=methylKit::unite(covm,min.per.group = 1L)
# meth
# 
# getCorrelation(meth,plot=TRUE)
# clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
# PCASamples(meth)

#OR, extract based on technical replicates 
samps = unique(name$Biosample)
outs = c('B2X7380_BL_CHK_F__SRR22581142','B3Y5802_BL_CHK_F__SRR20216566','B4H0631_BL_CHK_F__SRR20216481')
covs = c(10,25,50)
cordat = NULLG
for (out in outs) {
  for (cov in covs) {
    for (tech in samps) {
      cat('Working on sample: ',tech,'\n')
      sf = c(tech,out)
      fileuse = files[grepl(paste0(sf,collapse='|'),files)]
      if (length(fileuse) > 3 || length(fileuse) < 3) { next } 
      nameuse = name %>% filter(Biosample %in% sf | ID == out) %>% 
        mutate(Treatment = row_number()) 
      CpG.bismark.covM=methRead(location=as.list(fileuse), 
                                sample.id=as.list(nameuse$ID),
                                treatment=nameuse$Treatment,assembly="bluetit",context="CpG",pipeline="bismarkCoverage",mincov = 10)
      covm=filterByCoverage(CpG.bismark.covM,lo.count=cov,lo.perc=NULL,hi.count=250,hi.perc=NULL)
      meth=methylKit::unite(covm)
      sink('tmp.txt')
      getCorrelation(meth,method='spearman',plot=FALSE,nrow=2e1000)
      sink()
      cd = read.table('tmp.txt',comment.char='')
      cd$ID = row.names(cd)  
      co = cd %>% pivot_longer(!c(ID))
      names(co) = c('ID_A','ID_B','Correlation')
      cc = gsub('.*__','',out)
      if (out == 'B2X7380_BL_CHK_F__SRR22581142') { outcov = 'High' } else if (out == 'B3Y5802_BL_CHK_F__SRR20216566') { outcov = 'Intermediate' } else { outcov = 'Low'} 
      cof = co %>% 
        rowwise() %>% 
        mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
        group_by(pair) %>%
        distinct(pair, .keep_all = T) %>% 
        separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% dplyr::select(-pair) %>% filter(ID_A != ID_B) %>% 
        mutate(Comparison = ifelse(grepl(cc,ID_A),'Interindividual',ifelse(grepl(cc,ID_B),'Interindividual','Technical')))
      cof$CpGs = nrow(meth)
      cof$Min_Coverage = cov
      cof$Out_Coverage = outcov
      cordat = rbind(cordat,cof)
    }
  }
}
write.table(cordat,file='../Correlations_Methylkit.txt',quote=F,sep='\t',row.names=F)
cordat = read.table('../Correlations_Methylkit.txt',header=TRUE)
cs = 
rep = cordat %>% ggplot(aes(x=Comparison,y=Correlation,fill=Comparison))+
  geom_violin(alpha=0.5,draw_quantiles = c(0.025,0.5,0.975))+
  geom_jitter(alpha=0.2)+
  scale_fill_manual(values=c('darkblue','darkorchid2'))+
  facet_grid(Min_Coverage ~ Out_Coverage,scales='free')+
  theme_bw()+
  theme(legend.position='none')
rep

png('../Repeatability_Interindividual_Variations.png',units='in',res=600,height=5,width=6)
rep
dev.off()


