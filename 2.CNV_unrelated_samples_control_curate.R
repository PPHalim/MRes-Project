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

#data of unrelated individual have 2 versions: 
#1 from the home directory, 1 from labkey v8
#both from aggregate_illumina_gvcf database


#copy data from labkey #change API code everytime
labkey.setDefaults(apiKey="session|d064191cafe66db048cad1275cb2f8b7")
labkey.setDefaults(baseUrl="http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/")

UnRInd1 <- labkey.selectRows(
  baseUrl="http://emb-prod-mre-labkey-01.gel.zone:8080/labkey",
  folderPath="/main-programme/main-programme_v8_2019-11-28",
  schemaName="lists",
  queryName="aggregate_gvcf_sample_stats",
  viewName="",
  colFilter=makeFilter(c("set_of_unrelated_individuals", "EQUAL", "1")),
  containerFilter=NULL
)
#only 37959 data remained in labkey v8
UnRInd2 <- read.csv("SVCNV/aggregate_illumina_gvcf_relatedness/relatedness/60k_HWE_30k_random_unrelated_participants.txt", header = FALSE, col.names = "Plate.key")
                    
#crossref for unrelatedness
names(UnRInd1)[2] <- "Plate.key"
UnRIndconf <- inner_join(UnRInd1,UnRInd2, by= "Plate.key")
#all 37959 data in labkey v8 are still in accordance to previous version

#get platekey of patient/parents from immune domain
Immune_patient <- labkey.selectRows(
  baseUrl="http://emb-prod-mre-labkey-01.gel.zone:8080/labkey",
  folderPath="/main-programme/main-programme_v8_2019-11-28",
  schemaName="lists",
  queryName="domain_assignment",
  viewName="",
  colFilter=makeFilter(c("domain", "EQUAL", "Immune disorders")),
  containerFilter=NULL
)

Disease_analysis <- labkey.selectRows(
  baseUrl="http://emb-prod-mre-labkey-01.gel.zone:8080/labkey",
  folderPath="/main-programme/main-programme_v8_2019-11-28",
  schemaName="lists",
  queryName="rare_disease_analysis",
  viewName="",
  colFilter=makeFilter(c("genome_build", "EQUAL", "GRCh38")),
  containerFilter=NULL
)

Patient.data <- inner_join(Immune_patient, Disease_analysis)
names(Patient.data)[8] <- "Plate.key"
Patient.data <-Patient.data[,c(1,2,3,4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]

#crossref to eliminate any overlap with PID patients- match Plate.key
Control_Group <- anti_join(UnRIndconf, Patient.data)
#37173 people included in control group, 1110 people have paternal and maternal plate.key --> confirm unrelatedness of paternal and maternal DNA
names(Control_Group)[2]<- "patient.plate.key"
names(Control_Group)[7]<- "Plate.key"
Control_pa <- anti_join(Control_Group, Patient.data)#37173
names(Control_Group)[7:8]<- c("Paternal.plate.key","Plate.key")
Control_ma <- anti_join(Control_Group, Patient.data)#37173
names(Control_ma)[8]<- "Maternal.plate.key"
names(Control_ma)[2]<- "Plate.key"
names(Control_pa)[7:8]<- c("Paternal.plate.key","Maternal.plate.key")
names(Control_pa)[2]<- "Plate.key"
Control_Group <- inner_join(Control_ma, Control_pa)
write.table(Control_Group, file = "SVCNV/CNV_in_Control_unrelated/Control_Group_database.txt", sep = "\t")

Control_Group_clean <- Control_Group[,c(1,2,4,5,6)]

#importing CNV data list and adding header #count CNV loss/gain length
CNVdataloss <- read.csv("CNV-internal/RD_Canvas_LOSS_WEIGHTED.2018-07-26.bed", header = FALSE, col.names = c("Chromosome", "start", "end", "Plate.key", "copy.number", "dbl.1", "result", "algorithm:specification", "gender"), sep = "")
CNVdataloss$length.loss <- CNVdataloss$end-CNVdataloss$start
CNVdataloss <- CNVdataloss[order(CNVdataloss$Plate.key),c(4,1,2,3,10,5,6,7,8,9)]

CNVdatagain <- read.csv("CNV-internal/RD_Canvas_GAIN.2018-07-26.bed", header = FALSE, col.names = c("Chromosome", "start", "end", "Plate.key", "copy.number", "dbl.1", "result", "algorithm:specification", "gender"), sep = "")
CNVdatagain$length.gain <- CNVdatagain$end-CNVdatagain$start
CNVdatagain <- CNVdatagain[order(CNVdatagain$Plate.key),c(4,1,2,3,10,5,6,7,8,9)]

CNVdataloss <- unique(CNVdataloss)
CNVdatagain <- unique(CNVdatagain)

CNVdataloss$condition <- rep("loss")
CNVdatagain$condition <- rep("gain")
CNVdatafull <- full_join(CNVdatagain, CNVdataloss)
CNVdatafull <- CNVdatafull[order(CNVdatafull$Plate.key),c(1,2,3,4,11,5,12,6,7,9)]
#the starting position of nucleotide is still using 0 base coordinate -- need to be changed into 1 base coordinate just in case it affect the result
CNVdatafull$start <- (CNVdatafull$start+1)

#get CNV in control group
CNV.in.control <- inner_join(Control_Group_clean, CNVdatafull, by= "Plate.key")
length(unique(CNV.in.control$Plate.key))#18592
#1207640 CNVs found from 18592 patients matched
write.table (CNV.in.control, file = "SVCNV/CNV_in_Control_unrelated/Control_Group_18592peoplewithCNV.txt",sep = "\t")

