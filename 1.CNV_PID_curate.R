library(Rlabkey)
library(tidyverse, lib="~/re_gecip/immune/phalim/R_packages")
library(curl)
library(rjson)
library(RCurl)
library(dplyr)
library(purrr)
library(plyr)

#setting working directory
setwd("/home/phalim/re_gecip/immune/phalim/")

#copy data from labkey #change API code everytime
labkey.setDefaults(apiKey="session|711a0b4199f303041fa9a51789ef317b")
labkey.setDefaults(baseUrl="http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/")


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
#get Platekey data to be matched with patient.code
Patient.data <- inner_join(Immune_patient, Disease_analysis)
Patient.data <-Patient.data[,c(1,2,3,4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
names(Patient.data)[6]<-c("Plate.key")
#not all immune patient registered in rare disease analysis

####2.	If they are using two tools, there should be two CNV/SV output files per patient: can you locate and read these files into R?

#importing CNV data list and adding header #count CNV loss/gain length
CNVdataloss <- read.csv("CNV-internal/RD_Canvas_LOSS_WEIGHTED.2018-07-26.bed", header = FALSE, col.names = c("Chromosome", "start", "end", "Plate.key", "dbl", "dbl.1", "result", "algorithm:specification", "gender"), sep = "")
CNVdataloss$length.loss <- CNVdataloss$end-CNVdataloss$start
CNVdataloss <- CNVdataloss[order(CNVdataloss$Plate.key),c(4,1,2,3,10,5,6,7,8,9)]

CNVdatagain <- read.csv("CNV-internal/RD_Canvas_GAIN.2018-07-26.bed", header = FALSE, col.names = c("Chromosome", "start", "end", "Plate.key", "dbl", "dbl.1", "result", "algorithm:specification", "gender"), sep = "")
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

#-> no need because starting position also counted as length
CNVdatafull$length.gain <- (CNVdatafull$length.gain-1) 
CNVdatafull$length.loss <- (CNVdatafull$length.loss-1)


CNV.in.immune <- inner_join(Patient.data, CNVdatafull, by= "Plate.key")
CNV.in.immune <- CNV.in.immune[order(CNV.in.immune$Plate.key), c(1:28)]
Immune.CNV.patient <- data.frame(unique(CNV.in.immune$Plate.key),unique(CNV.in.immune$Path))

CNV.in.immune$Chrloc.start <- paste(CNV.in.immune$Chromosome,(CNV.in.immune$start), sep = ":")
CNV.in.immune$Chrloc.end <- paste(CNV.in.immune$Chromosome,(CNV.in.immune$end), sep = ":")

write.table(CNV.in.immune, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_immune_domain/List_immune_patient_with_CNV.txt", sep="\t", row.names = F)
write.table(Immune.CNV.patient, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_immune_domain/List_immune_patient_pathonly.txt", sep="\t", row.names = F)

#3.	Each patient should have a list of CNVs with their location. Calculate the size of each CNV and:
#(a) plot a histogram of the CNV sizes (with log10 transformation); (b) print a summary of the sizes using summary()

par(mfrow=c(1,1))
Plot1 = hist(log10(CNVdatafull$length.gain), freq = NULL, xlim = c(3,9), xlab="Data", breaks = 30, main = "Gained Nucleotides")
Plot2 = hist(log10(CNVdatafull$length.loss), freq = NULL, xlim = c(3,9), xlab="Data", breaks = 30, main = "Loss Nucleotides")

summary (CNVdatafull$length.gain)
summary (CNVdatafull$length.loss)

boxplot(20*log10(CNVdatafull[,6:7]), main = "Overall Size Distribution")

length(unique(CNVdatafull$Plate.key))
#[1] 28676
list(unique(CNVdatafull$Plate.key))
#first 10 -
CNVdatafull <- read.csv("SVCNV/CNVimmuneresult/List_immune_patient_with_CNV.txt", sep="\t")


#4.	Group the CNVs of 10 different patients in a list() and represent the size distribution of all of them as boxplots.
splitCNV <- split(CNVdatafull, CNVdatafull$Plate.key)
t1 <- as.data.frame(splitCNV[[1]])
t2 <- as.data.frame(splitCNV[[2]])
t3 <- as.data.frame(splitCNV[[3]])
t4 <- as.data.frame(splitCNV[[4]])
t5 <- as.data.frame(splitCNV[[5]])
t6 <- as.data.frame(splitCNV[[6]])
t7 <- as.data.frame(splitCNV[[7]])
t8 <- as.data.frame(splitCNV[[8]])
t9 <- as.data.frame(splitCNV[[9]])
t10 <- as.data.frame(splitCNV[[10]])
#combine the list of first 10 patient
# make list --> 
CNV10 <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)
CNV10 <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10) %>% reduce(full_join, by= c("Plate.key", "Chromosome", "start", "end", "condition", "algorithm.specification"))
CNV10 <- CNV10[c(1,2,3,4,5,8,6,7,9:26)]
names(CNV10)[7:26] <- c("gain.1","loss.1","gain.2","loss.2","gain.3","loss.3",
                         "gain.4","loss.4","gain.5","loss.5","gain.6","loss.6","gain.7","loss.7",
                         "gain.8","loss.8","gain.9","loss.9","gain.10","loss.10")

boxplot(CNV10[,7:26], main="10 Patients' boxplot")
