library(plyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(gdata)
library(grid)
library(readxl)
library(xlsx)

setwd("/Users/")

#Read in imaging phenotype data
pheno <- read.csv("20170321/processed_FreeSurfer_phenotypes.csv", header=T)
pheno <- subset(pheno, TIMEPT==1)

#Subset by ancestry (specify ancestry)
#ancestry <- read.table("20170321/FID_by_Ancestry/caucasians_ids.txt")
#colnames(ancestry) <- c("name", "IID")
#ancestry <- c(as.character(ancestry$IID))
#pheno <- subset(pheno, IID %in% ancestry)

#Read in combined ped files
geno_raw <- read.table("Documents/Projects/MGH/genetic_association_analyses/plot_data/combined_peds_imputed_bgn.txt")
colnames(geno_raw) <- c("FID", "IID", "MID", "PID", "sex", "pheno", "g1", "g2")

#Merge ped files and phenotype files by FID and IID
merged = merge(pheno,geno_raw,by=c("FID","IID"))
#Remove duplicates
merged_2 = subset(merged, !duplicated(IID))
#Select only (patient) groups 2,3 & 4 --> combined under group code 2 in processed data
merged_3 = subset(merged_2, merged_2$GROUP=="Schizophrenia" | merged_2$GROUP=="Bipolar Disorder" | merged_2$GROUP=="Major Depressive Disorder")
#Certain alleles have 0; remove these; might have been removed in filtering if subject had too many SNP missing
#merged_4 = subset(merged_3, merged_3$g1=="C" | merged_3$g1=="T")
merged_4 <- merged_3
#Do not include A,G,& 0 in levels of g1 and g2
merged_4$g1 <- droplevels(merged_4$g1)
merged_4$g2 <- droplevels(merged_4$g2)
#merged_4 <- merged_4a

#Combine TC & TT (rename g1 col values)
merged_4$genotype <- ifelse(merged_4$g1=="C" & merged_4$g2=="C", "CC",
                            ifelse(merged_4$g1=="C" & merged_4$g2=="T", "TC",
                                   ifelse(merged_4$g1=="T" & merged_4$g2=="C", "TC",
                                          ifelse(merged_4$g1=="T" & merged_4$g2=="T", "TT", NA))))

merged_4$genotype2 <- factor(as.character(merged_4$genotype), levels = c("CC", "TC", "TT"), labels = c("CC", "TC/TT", "TC/TT"))
#merged_4$g1 <- factor(as.character(merged_4$g1), levels = c("C", "T"), labels = c("CC", "TC/TT"))

#Import imaging dictionary
dict = read_excel("20170321/GENUS_processed_FreeSurfer_phenotypes_data_dictionary.xlsx")
#List of phenotypes to test
var_list = read.table("20170321/imaging_variable_list.txt")[,2]
#remove variables with rh, lh, and tot prefixes
var_list = var_list[!grepl('rh_|lh_|tot_|rh|lh',var_list)]

#Select phenotypes containing string
phenos = dict$Variable[grepl(paste(var_list,collapse="|"),dict$Variable)]
phenos = phenos[-c(length(phenos), length(phenos)-1)]

### for all phenotypes ###
#phenos = dict$Variable
#phenos <- subset(phenos, phenos %in% names(merged_4))
##########################

#Include only phenotype and allele columns
merged_5 <- subset(merged_4[,c("genotype", "genotype2", phenos)], !is.na(merged_4$genotype), stringsAsFactors=FALSE)

#Count number of TT & TC in combined TT/TC group
num_tt = c()
num_tc = c()
for (i in 1:length(phenos)){
  tt = 0
  tc = 0
  pheno = phenos[i]
  subtt = subset(merged_5, merged_5$genotype=="TT" & !is.na(merged_5[pheno]))
  subtc = subset(merged_5, merged_5$genotype=="TC" & !is.na(merged_5[pheno]))
  num_tt[i] = length(subtt$genotype)
  num_tc[i] = length(subtc$genotype)
}

#phenos <- c("\n avg_G_temp_sup.G_T_transv_area_D \n","avg_G_temp_sup.G_T_transv_curvind_D","avg_G_temp_sup.G_T_transv_foldind_D",
#            "avg_G_temp_sup.G_T_transv_gauscurv_D","avg_G_temp_sup.G_T_transv_meancurv_D","avg_G_temp_sup.G_T_transv_thickness_D",
#            "avg_G_temp_sup.G_T_transv_thicknessstd_D","avg_G_temp_sup.G_T_transv_volume_D","avg_G_temp_sup.Lateral_area_D",
#            "avg_G_temp_sup.Lateral_curvind_D","avg_G_temp_sup.Lateral_foldind_D","avg_G_temp_sup.Lateral_gauscurv_D",
#            "avg_G_temp_sup.Lateral_meancurv_D","avg_G_temp_sup.Lateral_thickness_D","avg_G_temp_sup.Lateral_thicknessstd_D","avg_G_temp_sup.Lateral_volume_D",
#            "avg_G_temp_sup.Plan_polar_area_D","avg_G_temp_sup.Plan_polar_curvind_D","avg_G_temp_sup.Plan_polar_foldind_D")

#Reshape dataframe
merged_5 <- melt(merged_5)
#Remove if phenotype value is NA
merged_5 <- merged_5[!is.na(merged_5$value), ]

#Compute p-values for t-test 
test_results <- numeric(length(phenos))
count=1
for (i in phenos){
  x1 = subset(merged_5, merged_5$variable==i)
  if (!all(is.na(x1$value))){
    out <- t.test(value ~ genotype2, data = x1)
    test_results[count] = out$p.value
    count = count+1
  } else {
    count = count+1
    next
  }
}
test_results <- data.frame(test_results)

#Append significant phenotype p-values to phenos_b and highly significant values to phenos_a
phenos_a = c()
phenos_b = c()
counta=1
countb=1
for (i in 1:length(phenos)) {
  if (test_results[i,]<0.005) {
    phenos_a[counta] = phenos[i]
    counta=counta+1
  }
  else if (test_results[i,]<0.05 & test_results[i,]>0.005) {
    phenos_b[countb] = phenos[i]
    countb=countb+1
  }
}

#Write phenotypes with corresponding p-values to excel file
t_data <- data.frame(phenos, test_results)
write.xlsx(t_data, file='data.xlsx')

#load phenos with only statistically significant phenotypes
phenos <- c(phenos_a,phenos_b)

# For each phenotype, count number of patients in CC group & TC/TT group
a1 = c()
for (i in 1:length(phenos)){
  x = 0
  for (j in 1:length(merged_5$variable)){
    if (merged_5$variable[j] == phenos[i] & merged_5$genotype2[j] == "CC"){
      x = x+1
    }
  }
  a1[i]=x
}

a2 = c()
for (i in 1:length(phenos)){
  x = 0
  for (j in 1:length(merged_5$variable)){
    if (merged_5$variable[j] == phenos[i] & merged_5$genotype2[j] == "TC/TT"){
      x = x+1
    }
  }
  a2[i]=x
}

### use to remove phenotypes values lying beyond range of the rest ###

#phenos_a <- c()
#phenos_b <- phenos_b[11]

######################################################################

#Drop statistically insignificant phenotypes from merged_5
merged_5 <- subset(merged_5, variable %in% phenos)
merged_5$variable <- droplevels(merged_5$variable)

#Change variable names
final_phenos_a <- phenos_a
if (length(phenos_a) != 0){
  phenos_a <- as.character(1:length(final_phenos_a))
}
final_phenos_b <- phenos_b
phenos_b <- as.character((1+length(final_phenos_a)):(length(final_phenos_a) + length(final_phenos_b)))

final_phenos = c(final_phenos_a, final_phenos_b)

#Labels of statistically significant phenotypes
final = dict$Label[grepl(paste(final_phenos,collapse="|"),dict$Variable)]
final_a = dict$Label[grepl(paste(final_phenos_a,collapse="|"),dict$Variable)]
final_b = dict$Label[grepl(paste(final_phenos_b,collapse="|"),dict$Variable)]

#manually adjusted labels
final = c('Avg Frontal Pole CT','Avg Inferior Frontal CT','Avg Lateral Orbitofrontal CT','Avg Lingual Volume', 'AVG Medial Orbitofrontal CT', 'Avg Pars Opercularis CT', 'Avg Rostral Middle Frontal Volume')
final_a = c('Avg Frontal Pole CT','Avg Inferior Frontal CT','Avg Lingual Volume')
final_b = c('Avg Lateral Orbitofrontal CT','Avg Lingual Volume', 'AVG Medial Orbitofrontal CT', 'Avg Rostral Middle Frontal Volume')

#levels(merged_5$variable) <- 1:length(levels(merged_5$variable))

levels(merged_5$variable) <- final

#Plot data
if (length(phenos_a) != 0){
  coordinates_a <- data.frame(variable=sort(rep(final_a,4)), x=rep(c(1,1,2,2),length(final_a)), y=rep(c(4,4.22,4.22,4),length(final_a)))
}
if (length(phenos_b) != 0){
  coordinates_b <- data.frame(variable=sort(rep(final_b,4)), x=rep(c(1,1,2,2),length(final_b)), y=rep(c(4,4.22,4.22,4),length(final_b)))
}

size <- data.frame(variable=final, labela=a1, labelb=a2)

#Remove datapoints not in range -3 to 3
#merged_5 <- subset(merged_5, merged_5$value>-3 & merged_5$value<3)

g <- ggplot(merged_5, aes_string(x="genotype2", y="value"))+
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE)+
  coord_cartesian(ylim = c(-5.5, 5))+ 
  facet_wrap(~variable, ncol = 6, scales="free_x") +
  theme_bw()+ 
  #theme(axis.title.x = element_blank(), strip.background = element_rect(fill = "white"), 
  #      strip.text.x=element_text(margin = margin(0,0,0,0, "cm"),size=10),legend.position = c(4, 4), legend.justification = c(4, 4))+
  labs(caption = "*domain Z-Score")+
  geom_path(data = coordinates_a, aes(x = x, y = y), linetype = 1, size = 0.3)+
  geom_text(data = coordinates_a, aes(x=1.5, y=4.5, label="p<0.005"), colour="black", 
            inherit.aes=FALSE, parse=FALSE, size = 3)+
  geom_path(data = coordinates_b, aes(x = x, y = y), linetype = 1, size = 0.3)+
  geom_text(data = coordinates_b, aes(x=1.5, y=4.5, label="p<0.05"), colour="black", 
            inherit.aes=FALSE, parse=FALSE, size = 3)+
  geom_text(data = size, aes(x=1, y=-5.8, label=labela), colour="black", 
            inherit.aes=FALSE, parse=FALSE, size = 3)+
  geom_text(data = size, aes(x=2, y=-5.8, label=labelb), colour="black", 
            inherit.aes=FALSE, parse=FALSE, size = 3)

ggsave(filename = "plot.pdf", g, height = 24, width = 40, units="cm", dpi=300)
