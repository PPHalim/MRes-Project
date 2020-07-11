library(Rlabkey)
library(tidyverse, lib="~/re_gecip/immune/phalim/R_packages")
library(curl)
library(rjson)
library(RCurl)
library(dplyr)
library(purrr)
library(stringr)

#setting working directory
setwd("/home/phalim/re_gecip/immune/phalim/")

#match the CNV in control group to the PID genes
#PID gene ranges
PIDgenerange <-read.csv("SVCNV/CNV_in_PID_immune_gene/PID_genes_range.txt", sep = "\t")
CNV.in.control <-read.csv(file = "SVCNV/CNV_in_Control_unrelated/Control_Group_18592peoplewithCNV.txt",sep = "\t")

#take only the chromosome information needed for genomic range analysis
##seperate the rare disease and cancer germline samples
CNV.in.control_list <- split(CNV.in.control, CNV.in.control$Type)

CNV.in.control_cancer <- as.data.frame(CNV.in.control_list[[1]])#377319 CNVs
length(unique(CNV.in.control_cancer$Participant.Id))#5783 control from cancer germline

CNV.in.control_rare <- as.data.frame(CNV.in.control_list[[2]])#830321 CNVs
length(unique(CNV.in.control_rare$Participant.Id))#12809 control from rare disease

###For cancer germline samples
CNVcontrolcancer <- CNV.in.control_cancer[,c(2,4:13)]
CNVcontrolcancer$Chromosome <- gsub("chr","",as.character(CNVcontrolcancer$Chromosome))

#overlap within the genes
Immunegene_CNV_cancer1 <- inner_join(CNVcontrolcancer, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x)) #33
Immunegene_CNV_cancer1 <- unique(Immunegene_CNV_cancer1)#33
Immunegene_CNV_cancer1$match <- rep("within")
#overlap at the start of the genes
Immunegene_CNV_cancer2 <- inner_join(CNVcontrolcancer, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x)) #4465
Immunegene_CNV_cancer2 <- unique(Immunegene_CNV_cancer2)#4465

#overlap at the end of the genes
Immunegene_CNV_cancer3 <- inner_join(CNVcontrolcancer, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
Immunegene_CNV_cancer3 <- unique(Immunegene_CNV_cancer3)#4802

#overlap counted twice in Immunegene_CNV_cancer2 and Immunegene_CNV_cancer3
ImmuneGene_CNV_cancerwithinCNV <- semi_join(Immunegene_CNV_cancer2,Immunegene_CNV_cancer3)#403
#4254 CNVs larger than gene
Immunegene_CNV_canceratstart <- anti_join(Immunegene_CNV_cancer2,Immunegene_CNV_cancer3) #8
#221 CNVs at the start of the gene
Immunegene_CNV_canceratend <- anti_join(Immunegene_CNV_cancer3,Immunegene_CNV_cancer2) #33
#548 CNVs at the end of the gene
Immunegene_CNV_canceratstart$match <- rep("start")
Immunegene_CNV_canceratend$match <- rep("end")

splitdeladdCNVcancer <- split(ImmuneGene_CNV_cancerwithinCNV, ImmuneGene_CNV_cancerwithinCNV$copy.number)
sc0 <- as.data.frame(splitdeladdCNVcancer[[1]])
sc1 <- as.data.frame(splitdeladdCNVcancer[[2]])
sc2 <- as.data.frame(splitdeladdCNVcancer[[3]])
sc3 <- as.data.frame(splitdeladdCNVcancer[[4]])
sc4 <- as.data.frame(splitdeladdCNVcancer[[5]])
sc5 <- as.data.frame(splitdeladdCNVcancer[[6]])
ImmunegeneCNVdelc <- full_join(sc0,sc1)
ImmunegeneCNVdelc$match <- rep("del")
ImmunegeneCNVaddc <- full_join(sc3,full_join(sc4,sc5))
ImmunegeneCNVaddc$match <- rep("add")
sc2$match <- rep("CNV>gene, no change of copy number")
####check total match with immunegeneCNV####
ImmuneGene_CNV_cancertotal <- full_join(Immunegene_CNV_cancer1,full_join(ImmunegeneCNVaddc, full_join(sc2, full_join(ImmunegeneCNVdelc, full_join(Immunegene_CNV_canceratstart, Immunegene_CNV_canceratend)))))

#rearrange ImmuneGene_CNV_cancertotal row
ImmuneGene_CNV_cancertotal <- ImmuneGene_CNV_cancertotal[order(ImmuneGene_CNV_cancertotal$Plate.key, ImmuneGene_CNV_cancertotal$Chromosome, ImmuneGene_CNV_cancertotal$condition),c(1:15)]
length(unique(ImmuneGene_CNV_cancertotal$Plate.key))#2537
#rearrange ImmuneGene_CNV_cancertotal collumn
ImmuneGene_CNV_cancertotal <- ImmuneGene_CNV_cancertotal[c(1,4,7:11,5,6,12,13,2,3,14,15)]
names(ImmuneGene_CNV_cancertotal)[7:13] <- c("read.depth","start","end","gene.start","gene.end","sex.phenotype","sample.type")
#save
write.table(ImmuneGene_CNV_cancertotal, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_2537cancersample.txt", sep= "\t", row.names = F)


###For rare disease samples
CNVcontrolrare <- CNV.in.control_rare[,c(2,4:13)]
CNVcontrolrare$Chromosome <- gsub("chr","",as.character(CNVcontrolrare$Chromosome))

#overlap within the genes
Immunegene_CNV_rare1 <- inner_join(CNVcontrolrare, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x)) #110
Immunegene_CNV_rare1 <- unique(Immunegene_CNV_rare1)#110
Immunegene_CNV_rare1$match <- rep("within")
#overlap at the start of the genes
Immunegene_CNV_rare2 <- inner_join(CNVcontrolrare, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x)) #11118
Immunegene_CNV_rare2 <- unique(Immunegene_CNV_rare2)#11118

#overlap at the end of the genes
Immunegene_CNV_rare3 <- inner_join(CNVcontrolrare, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x)) #11700
Immunegene_CNV_rare3 <- unique(Immunegene_CNV_rare3)#11700

#overlap counted twice in Immunegene_CNV_rare2 and Immunegene_CNV_rare3
ImmuneGene_CNV_rarewithinCNV <- semi_join(Immunegene_CNV_rare2,Immunegene_CNV_rare3)#10515
#10515 CNVs larger than gene
Immunegene_CNV_rareatstart <- anti_join(Immunegene_CNV_rare2,Immunegene_CNV_rare3) #8
#603 CNVs at the start of the gene
Immunegene_CNV_rareatend <- anti_join(Immunegene_CNV_rare3,Immunegene_CNV_rare2) #33
#1185 CNVs at the end of the gene
Immunegene_CNV_rareatstart$match <- rep("start")
Immunegene_CNV_rareatend$match <- rep("end")

splitdeladdCNVrare <- split(ImmuneGene_CNV_rarewithinCNV, ImmuneGene_CNV_rarewithinCNV$copy.number)
sr0 <- as.data.frame(splitdeladdCNVrare[[1]])
sr1 <- as.data.frame(splitdeladdCNVrare[[2]])
sr2 <- as.data.frame(splitdeladdCNVrare[[3]])
sr3 <- as.data.frame(splitdeladdCNVrare[[4]])
sr4 <- as.data.frame(splitdeladdCNVrare[[5]])
sr5 <- as.data.frame(splitdeladdCNVrare[[6]])
ImmunegeneCNVdelr <- full_join(sr0,sr1)
ImmunegeneCNVdelr$match <- rep("del")
ImmunegeneCNVaddr <- full_join(sr3,full_join(sr4,sr5))
ImmunegeneCNVaddr$match <- rep("add")
sr2$match <- rep("CNV>gene, no change of copy number")
####check total match with immunegeneCNV####
ImmuneGene_CNV_raretotal <- full_join(Immunegene_CNV_rare1,full_join(ImmunegeneCNVaddr, full_join(sr2, full_join(ImmunegeneCNVdelr, full_join(Immunegene_CNV_rareatstart, Immunegene_CNV_rareatend)))))

#rearrange ImmuneGene_CNV_raretotal row
ImmuneGene_CNV_raretotal <- ImmuneGene_CNV_raretotal[order(ImmuneGene_CNV_raretotal$Plate.key, ImmuneGene_CNV_raretotal$Chromosome, ImmuneGene_CNV_raretotal$condition),c(1:15)]
length(unique(ImmuneGene_CNV_raretotal$Plate.key))#5786 patient
#rearrange ImmuneGene_CNV_raretotal collumn
ImmuneGene_CNV_raretotal <- ImmuneGene_CNV_raretotal[c(1,4,7:11,5,6,12,13,2,3,14,15)]
names(ImmuneGene_CNV_raretotal)[7:13] <- c("read.depth","start","end","gene.start","gene.end","sex.phenotype","sample.type")
#save
write.table(ImmuneGene_CNV_raretotal, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_5786raresample.txt", sep= "\t", row.names = F)

-----#R session aborted on and on. Thus, data must be imported from saved file
  
ImmuneGene_CNV_cancertotal <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_2537cancersample.txt", sep= "\t")
ImmuneGene_CNV_raretotal <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_5786raresample.txt", sep= "\t")
#rejoin data from rare and cancer
Immunegene_CNV_control <- full_join(ImmuneGene_CNV_raretotal, ImmuneGene_CNV_cancertotal)
write.table(Immunegene_CNV_control, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_8323_allsample.txt", sep= "\t", row.names = F)


