library(plyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(gdata)
library(grid)

#Read in neuropsych phenotype data
pheno_raw <- read.csv("20170321/GENUS_neuropsych_data_raw.csv", header=T)
pheno_processed <- read.csv("20170321/GENUS_neuropsych_data_processed.csv", header=T)
pheno_combined <- merge(pheno_raw, pheno_processed, by="GENUS_ID")
pheno_combined <- subset(pheno_combined, TIMEPT==1)

#Replace repeated FID & IID column names
pheno_combined <- rename.vars(pheno_combined, from=c("FID.y", "IID.y"), to=c("FID", "IID"))

#Subset by ancestry
ancestry <- read.table("20170321/FID_by_Ancestry/eastasians_ids.txt")
colnames(ancestry) <- c("name", "IID")
ancestry <- c(as.character(ancestry$IID))
pheno_combined <- subset(pheno_combined, IID %in% ancestry)

#Read in combined ped files
geno_raw <- read.table("documents/projects/MGH/plot_data/combined_peds_imputed_bgn.txt")
colnames(geno_raw) <- c("FID", "IID", "MID", "PID", "sex", "pheno", "g1", "g2")

#Merge ped files and phenotype files by FID and IID
merged = merge(pheno_combined,geno_raw,by=c("FID","IID"))
#Remove duplicates
merged_2 = subset(merged, !duplicated(IID))
#Select only (patient) groups 2,3 & 4 --> combined under group code 2 in processed data
merged_3 = subset(merged_2, merged_2$GROUP.y==2 | merged_2$GROUP.y==3 | merged_2$GROUP.y==4)
#Certain alleles have 0; remove these; might have been removed in filtering if subject had too many SNP missings
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

phenos = c("BD_COMPOSITE_AGESEX","BVMT_COMPOSITE_AGESEX","CPTIP_COMPOSITE_AGESEX",
           "LNS_COMPOSITE_AGESEX","SC_COMPOSITE_AGESEX","SPSP_COMPOSITE_AGESEX",
           "TMTA_COMPOSITE_AGESEX","TMTB_COMPOSITE_AGESEX","VF_COMPOSITE_AGESEX",
           "WLLT_COMPOSITE_AGESEX","ATVI_DOMAIN_AGESEX","NVLM_DOMAIN_AGESEX",
           "NVWM_DOMAIN_AGESEX","RPS_DOMAIN_AGESEX","SOP_DOMAIN_AGESEX","VISPA_DOMAIN_AGESEX",
           "VLM_DOMAIN_AGESEX","VWM_DOMAIN_AGESEX","g_unimp_AGESEX")

#Include only phenotype and allele columns
merged_5 <- subset(merged_4[,c("genotype", "genotype2", phenos)], !is.na(merged_4$genotype))

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

#Rename phenotype names 
colnames(merged_5) <- c("genotype", "genotype2", "\n WAIS Block \n Design Scaled \n","Brief \n Visiospatial \n Memory Test","Continuous \n Performance \n Test",
                        "WAIS Letter- \n Number \n Sequencing","BACS \n Symbol \n Coding","WMS \n Spatial \n Span",
                        "Trial \n Making \n Test A","Trial \n Making \n Test B","Category \n Fluency \n Composite",
                        "Word List \n Learning \n Tests","Attention/ \n Vigilance*","Non-Verbal \n Learning/ \n Memory*",
                        "Non-Verbal \n Working \n Memory*","Reasoning/ \n Problem \n Solving*","Speed of \n Processing*","Visuo- \n Spatial \n Ability*",
                        "Verbal \n Learning/ \n Memory*","Verbal \n Working \n Memory*",  "Spearman's \n g")

phenos <- c("\n WAIS Block \n Design Scaled \n","Brief \n Visiospatial \n Memory Test","Continuous \n Performance \n Test",
            "WAIS Letter- \n Number \n Sequencing","BACS \n Symbol \n Coding","WMS \n Spatial \n Span",
            "Trial \n Making \n Test A","Trial \n Making \n Test B","Category \n Fluency \n Composite",
            "Word List \n Learning \n Tests","Attention/ \n Vigilance*","Non-Verbal \n Learning/ \n Memory*",
            "Non-Verbal \n Working \n Memory*","Reasoning/ \n Problem \n Solving*","Speed of \n Processing*","Visuo- \n Spatial \n Ability*",
            "Verbal \n Learning/ \n Memory*","Verbal \n Working \n Memory*",  "Spearman's \n g")

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

#Append significant phenotype p-values to phenos_a and insignificant values to phenos_b
phenos_a = c()
phenos_b = c()
counta=1
countb=1
for (i in 1:length(phenos)) {
  if (test_results[i,]<0.05) {
    phenos_a[counta] = phenos[i]
    counta=counta+1
  }
  else {
    phenos_b[countb] = phenos[i]
    countb = countb+1
  }
}

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

coordinates_a <- data.frame(variable=sort(rep(phenos_a,4)), x=rep(c(1,1,2,2),length(phenos_a)), y=rep(c(3.5,4.5,4.5,3.5),length(phenos_a)))
coordinates_b <- data.frame(variable=sort(rep(phenos_b,4)), x=rep(c(1,1,2,2),length(phenos_b)), y=rep(c(3.5,4.5,4.5,3.5),length(phenos_b)))
size <- data.frame(variable=phenos, labela=a1, labelb=a2)

#Remove datapoints not in range -3 to 3
merged_5 <- subset(merged_5, merged_5$value>-3 & merged_5$value<3)

if (nrow(coordinates_a)==0){
  g <- ggplot(merged_5, aes_string(x="genotype2", y="value"))+
    geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE)+
    coord_cartesian(ylim = c(-4, 5.2))+ 
    facet_wrap(~variable, ncol = 6, scales="free_x") +
    theme_bw()+ 
    theme(axis.title.x = element_blank(), strip.background = element_rect(fill = "white"), 
          strip.text.x=element_text(margin = margin(0,0,0,0, "cm"),size=10),legend.position = c(4, 4), legend.justification = c(4, 4))+
    labs(caption = "*domain Z-Score")+
    #geom_path(data = coordinates_a, aes(x = x, y = y), linetype = 1, size = 0.3)+
    #geom_text(data = coordinates_a, aes(x=1.5, y=5, label="p<0.05"), colour="black", 
    #         inherit.aes=FALSE, parse=FALSE, size = 3)+
    geom_path(data = coordinates_b, aes(x = x, y = y), linetype = 1, size = 0.3)+
    geom_text(data = coordinates_b, aes(x=1.5, y=5, label="ns"), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)+
    geom_text(data = size, aes(x=1, y=-3.8, label=labela), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)+
    geom_text(data = size, aes(x=2, y=-3.8, label=labelb), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)
} else {
  g <- ggplot(merged_5, aes_string(x="genotype2", y="value"))+
    geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE)+
    coord_cartesian(ylim = c(-4, 5.2))+ 
    facet_wrap(~variable, ncol = 6, scales="free_x") +
    theme_bw()+ 
    theme(axis.title.x = element_blank(), strip.background = element_rect(fill = "white"), 
          strip.text.x=element_text(margin = margin(0,0,0,0, "cm"),size=10),legend.position = c(4, 4), legend.justification = c(4, 4))+
    labs(caption = "*domain Z-Score")+
    geom_path(data = coordinates_a, aes(x = x, y = y), linetype = 1, size = 0.3)+
    geom_text(data = coordinates_a, aes(x=1.5, y=5, label="p<0.05"), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)+
    geom_path(data = coordinates_b, aes(x = x, y = y), linetype = 1, size = 0.3)+
    geom_text(data = coordinates_b, aes(x=1.5, y=5, label="ns"), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)+
    geom_text(data = size, aes(x=1, y=-3.8, label=labela), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)+
    geom_text(data = size, aes(x=2, y=-3.8, label=labelb), colour="black", 
              inherit.aes=FALSE, parse=FALSE, size = 3)
}

ggsave(filename = "plot.png", g, height = 24, width = 18, units="cm", dpi=300)
