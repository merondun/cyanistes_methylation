library(tidyverse)
library(matrixStats)

# Path with cov.gz files
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2022_EpigeneticsTits/2023JULY/bismark_cov')

# Get a list of all .cov.gz files in the working directory
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)

# Initialize empty lists to store names and coverage data
master = NULL
cov = NULL

# Loop through all .cov.gz files
for (samp in files) {
  cat('Reading in file: ',samp,'\n') # Display current file being read

  # Read in .cov.gz file
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);
  names(tab) = c('chr','start','end','Percent_5mC','countM','countU')

  # Extract name for new variable, trim all the garbage, modify as necessary depending on your file extensions....
  name <- gsub('./','',samp); name <- gsub('_val_1.*','',name); name <- gsub('_trimmed.*','',name); name <- gsub('.cov.gz','',name)

  # Calculate quick coverage and %5mc summary statistics per sample and store them
  tab = tab %>% mutate(Coverage=countM+countU,ID=name)
  c1 <- tab %>% mutate(Discrete = cut(Coverage, breaks = c(0,10,25,100,Inf), right = F, labels = c('Below10x','Below25x','Below100x','Above100x'))) %>% group_by(Discrete,ID) %>% dplyr::count()
  c2 <- tab %>% mutate(Discrete = cut(Percent_5mC, breaks = c(-1,10,90,101), right = F, labels = c('Below10p','Below90p','Above90p'))) %>% group_by(Discrete,ID) %>% dplyr::count()
  c1$Variable <- 'Coverage'
  c2$Variable <- 'Percent_5mC'
  cc <- rbind(c1,c2)
  cov <- rbind(cov,cc)

  # Keep only positions with 10x coverage or higher ### IMPORTANT!
  tab <- subset(tab, Coverage > 9)
  tab$site <- paste0(tab$chr,'_',tab$start)
  tab$Percent_5mC <- tab$Percent_5mC/100

  # Grab only site / countM / countU / %5mC
  tab2 <- tab %>% select(site,Percent_5mC)
  names(tab2) <- c('site',name);

  cat('Saved to variable: ',name,'\n') # Indicate current file has been processed, then bind with all the data
  if(is.null(master)){
    master <- tab2
  }else{
    master <- full_join(master, tab2, by = "site")
  }
}

# Filter according to missingness. If there are more than 'missing_samples' NAs at a site, drop the site
missingness_threshold = 0.2
colcount = ncol(master) - 1
filtered = master %>% mutate(
  NAs = rowSums(is.na(.)),
  NAp = NAs / colcount) %>%
  as_tibble %>%
  filter(NAp <= 0.2) %>%
  select(-c(NAs,NAp))

write.table(master,file='../overlapped_cpgs/5mC_Data_Pilot_10x.txt',quote=F,sep='\t',row.names=FALSE)
write.table(cov,file='../overlapped_cpgs/Coverage-Sample_Statistics_10x.txt',quote=F,sep='\t',row.names=F)

#write file with missingness filters
write.table(filtered,file='../overlapped_cpgs/5mC_Data_Pilot_10x-MM2.txt',quote=F,sep='\t',row.names=FALSE)

