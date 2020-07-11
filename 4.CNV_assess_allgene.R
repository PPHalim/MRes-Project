library(Rlabkey)
library(tidyverse, lib="~/re_gecip/immune/phalim/R_packages")
library(curl)
library(rjson)
library(RCurl)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(egg)
library(grid)
library(gridExtra)

#setting working directory
setwd("/home/phalim/re_gecip/immune/phalim/")

##CASE##

#gene information downloaded from ensembl.com
Gene_information <- read.csv("SVCNV/Geneinfo_ensembl.txt")
GeneRangeDF <- Gene_information[,c(1:4)]
names(GeneRangeDF)[1:4] <- c("Chromosome","start","end","gene.name")
GeneRangeDF <- unique(GeneRangeDF)#67140

GeneRangeDF$Chromosome <- paste0("chr",GeneRangeDF$Chromosome)

#CNV in immune patient found from processing in SVCNVscript.R
PIDpatientCNV <- read.csv("SVCNV/CNV_in_immune_domain/List_immune_patient_with_CNV.txt", sep="\t")
PIDpsplit <- split.data.frame(PIDpatientCNV, PIDpatientCNV$condition)
PIDpatientloss <- as.data.frame(PIDpsplit[[2]])

#CHECK OVERLAP
#overlap within the genes
CNV_PIDp_LOSS1 <- inner_join(PIDpatientloss, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
CNV_PIDp_LOSS1 <- unique(CNV_PIDp_LOSS1)#13780
CNV_PIDp_LOSS1$match <- rep("within")
#overlap at the start of the genes
CNV_PIDp_LOSS2 <- inner_join(PIDpatientloss, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
CNV_PIDp_LOSS2 <- unique(CNV_PIDp_LOSS2)#196731

#overlap at the end of the genes
CNV_PIDp_LOSS3 <- inner_join(PIDpatientloss, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
CNV_PIDp_LOSS3 <- unique(CNV_PIDp_LOSS3)#196153

#overlap counted twice CNVs larger than gene
CNV_PIDp_LOSSwithinCNV <- semi_join(CNV_PIDp_LOSS2,CNV_PIDp_LOSS3)#179955
#CNVs at the start of the gene 
CNV_PIDp_LOSSatstart <- anti_join(CNV_PIDp_LOSS2,CNV_PIDp_LOSS3) #16776
#CNVs at the end of the gene
CNV_PIDp_LOSSatend <- anti_join(CNV_PIDp_LOSS3,CNV_PIDp_LOSS2) #16198

CNV_PIDp_LOSSatstart$match <- rep("start")
CNV_PIDp_LOSSatend$match <- rep("end")

CNV_PIDp_LOSStotal <- list(CNV_PIDp_LOSS1,CNV_PIDp_LOSSwithinCNV,CNV_PIDp_LOSSatstart,CNV_PIDp_LOSSatend) %>% reduce(full_join)
#226709
CNV_patients <- CNV_PIDp_LOSStotal[,c(1,6)]
CNV_patients <- unique(CNV_patients)#486

write.table(CNV_PIDp_LOSStotal, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_226709CNV_487patients.txt", sep="\t")
write.table(CNV_patients, file="/home/phalim/re_gecip/immune/phalim/SVCNV/List_of_486PIDpatients.txt", sep="\t")

###########################################
##CONTROL##
CNV.in.control <- read.table(file = "SVCNV/CNV_in_Control_unrelated/Control_Group_18592peoplewithCNV.txt",sep = "\t")
splitControl <- split.data.frame(CNV.in.control,CNV.in.control$condition)
CNV.control.loss <- as.data.frame(splitControl[[2]])#630646
CNV_control <- CNV.control.loss[,c(1:3)]
CNV_control <- unique(CNV_control)#18592
write.table(CNV_control, file="/home/phalim/re_gecip/immune/phalim/SVCNV/List_of_18592control_unrelated.txt", sep="\t")

Splitcontrolloss <- split.data.frame(CNV.control.loss, CNV.control.loss$Type)
CNV.loss.cancer <- as.data.frame(Splitcontrolloss[[1]])
CNV.loss.cancer1 <- CNV.loss.cancer[c(1:98000),]
CNV.loss.cancer2 <- CNV.loss.cancer[c(98001:197958),]
CNV.loss.rare <- as.data.frame(Splitcontrolloss[[2]])
CNV.loss.rare1 <- CNV.loss.rare[c(1:150000),]
CNV.loss.rare2 <- CNV.loss.rare[c(150001:300000),]
CNV.loss.rare3 <- CNV.loss.rare[c(300001:432688),]

#CHECK OVERLAP cancer
#overlap within the genes
CNV_cancer1_LOSS1 <- inner_join(CNV.loss.cancer1, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
CNV_cancer1_LOSS1 <- unique(CNV_cancer1_LOSS1)#42569
CNV_cancer1_LOSS1$match <- rep("within")
#overlap at the start of the genes
CNV_cancer1_LOSS2 <- inner_join(CNV.loss.cancer1, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
CNV_cancer1_LOSS2 <- unique(CNV_cancer1_LOSS2)#647278

#overlap at the end of the genes
CNV_cancer1_LOSS3 <- inner_join(CNV.loss.cancer1, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
CNV_cancer1_LOSS3 <- unique(CNV_cancer1_LOSS3)#643588

#overlap counted twice CNVs larger than gene
CNV_cancer1_LOSSwithinCNV <- semi_join(CNV_cancer1_LOSS2,CNV_cancer1_LOSS3)#590892
#CNVs at the start of the gene 
CNV_cancer1_LOSSatstart <- anti_join(CNV_cancer1_LOSS2,CNV_cancer1_LOSS3) #56386
#CNVs at the end of the gene
CNV_cancer1_LOSSatend <- anti_join(CNV_cancer1_LOSS3,CNV_cancer1_LOSS2) #52696

CNV_cancer1_LOSSatstart$match <- rep("start")
CNV_cancer1_LOSSatend$match <- rep("end")

CNV_cancer1_LOSStotal <- list(CNV_cancer1_LOSS1,CNV_cancer1_LOSSwithinCNV,CNV_cancer1_LOSSatstart,CNV_cancer1_LOSSatend) %>% reduce(full_join)

#overlap within the genes
CNV_cancer2_LOSS1 <- inner_join(CNV.loss.cancer2, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
CNV_cancer2_LOSS1 <- unique(CNV_cancer2_LOSS1)#42569
CNV_cancer2_LOSS1$match <- rep("within")
#overlap at the start of the genes
CNV_cancer2_LOSS2 <- inner_join(CNV.loss.cancer2, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
CNV_cancer2_LOSS2 <- unique(CNV_cancer2_LOSS2)#647278

#overlap at the end of the genes
CNV_cancer2_LOSS3 <- inner_join(CNV.loss.cancer2, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
CNV_cancer2_LOSS3 <- unique(CNV_cancer2_LOSS3)#643588

#overlap counted twice CNVs larger than gene
CNV_cancer2_LOSSwithinCNV <- semi_join(CNV_cancer2_LOSS2,CNV_cancer2_LOSS3)#590892
#CNVs at the start of the gene 
CNV_cancer2_LOSSatstart <- anti_join(CNV_cancer2_LOSS2,CNV_cancer2_LOSS3) #56386
#CNVs at the end of the gene
CNV_cancer2_LOSSatend <- anti_join(CNV_cancer2_LOSS3,CNV_cancer2_LOSS2) #52696

CNV_cancer2_LOSSatstart$match <- rep("start")
CNV_cancer2_LOSSatend$match <- rep("end")

CNV_cancer2_LOSStotal <- list(CNV_cancer2_LOSS1,CNV_cancer2_LOSSwithinCNV,CNV_cancer2_LOSSatstart,CNV_cancer2_LOSSatend) %>% reduce(full_join)
#742543
CNV_cancer_LOSStotal <- unique(full_join(CNV_cancer1_LOSStotal,CNV_cancer2_LOSStotal))

CNV_cancer <- CNV_cancer_LOSStotal[,c(1:3)]
CNV_cancer <- unique(CNV_cancer)#5783

write.table(CNV_cancer_LOSStotal, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_cancer_742543CNV_5783patients_recheck.txt", sep="\t")
write.table(CNV_cancer, file="/home/phalim/re_gecip/immune/phalim/SVCNV/List_of_5783cancerpatients.txt", sep="\t")

#CHECK OVERLAP rare
#overlap within the genes
CNV_rare1_LOSS1 <- inner_join(CNV.loss.rare1, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
CNV_rare1_LOSS1 <- unique(CNV_rare1_LOSS1)#44095
CNV_rare1_LOSS1$match <- rep("within")
#overlap at the start of the genes
CNV_rare1_LOSS2 <- inner_join(CNV.loss.rare1, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
CNV_rare1_LOSS2 <- unique(CNV_rare1_LOSS2)#713752
#overlap at the end of the genes
CNV_rare1_LOSS3 <- inner_join(CNV.loss.rare1, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
CNV_rare1_LOSS3 <- unique(CNV_rare1_LOSS3)#710792
#overlap counted twice CNVs larger than gene
CNV_rare1_LOSSwithinCNV <- semi_join(CNV_rare1_LOSS2,CNV_rare1_LOSS3)#652570
#CNVs at the start of the gene 
CNV_rare1_LOSSatstart <- anti_join(CNV_rare1_LOSS2,CNV_rare1_LOSS3) #61182
#CNVs at the end of the gene
CNV_rare1_LOSSatend <- anti_join(CNV_rare1_LOSS3,CNV_rare1_LOSS2) #58222
CNV_rare1_LOSSatstart$match <- rep("start")
CNV_rare1_LOSSatend$match <- rep("end")
CNV_rare1_LOSStotal <- list(CNV_rare1_LOSS1,CNV_rare1_LOSSwithinCNV,CNV_rare1_LOSSatstart,CNV_rare1_LOSSatend) %>% reduce(full_join)
#816069

#overlap within the genes
CNV_rare2_LOSS1 <- inner_join(CNV.loss.rare2, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
CNV_rare2_LOSS1 <- unique(CNV_rare2_LOSS1)#44802
CNV_rare2_LOSS1$match <- rep("within")
#overlap at the start of the genes
CNV_rare2_LOSS2 <- inner_join(CNV.loss.rare2, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
CNV_rare2_LOSS2 <- unique(CNV_rare2_LOSS2)#720544
#overlap at the end of the genes
CNV_rare2_LOSS3 <- inner_join(CNV.loss.rare2, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
CNV_rare2_LOSS3 <- unique(CNV_rare2_LOSS3)#718594
#overlap counted twice CNVs larger than gene
CNV_rare2_LOSSwithinCNV <- semi_join(CNV_rare2_LOSS2,CNV_rare2_LOSS3)
#CNVs at the start of the gene 
CNV_rare2_LOSSatstart <- anti_join(CNV_rare2_LOSS2,CNV_rare2_LOSS3)
#CNVs at the end of the gene
CNV_rare2_LOSSatend <- anti_join(CNV_rare2_LOSS3,CNV_rare2_LOSS2)
CNV_rare2_LOSSatstart$match <- rep("start")
CNV_rare2_LOSSatend$match <- rep("end")
CNV_rare2_LOSStotal <- list(CNV_rare2_LOSS1,CNV_rare2_LOSSwithinCNV,CNV_rare2_LOSSatstart,CNV_rare2_LOSSatend) %>% reduce(full_join)
#824611

#overlap within the genes
CNV_rare3_LOSS1 <- inner_join(CNV.loss.rare3, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
CNV_rare3_LOSS1 <- unique(CNV_rare3_LOSS1)
CNV_rare3_LOSS1$match <- rep("within")
#overlap at the start of the genes
CNV_rare3_LOSS2 <- inner_join(CNV.loss.rare3, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
CNV_rare3_LOSS2 <- unique(CNV_rare3_LOSS2)
#overlap at the end of the genes
CNV_rare3_LOSS3 <- inner_join(CNV.loss.rare3, GeneRangeDF, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
CNV_rare3_LOSS3 <- unique(CNV_rare3_LOSS3)
#overlap counted twice CNVs larger than gene
CNV_rare3_LOSSwithinCNV <- semi_join(CNV_rare3_LOSS2,CNV_rare3_LOSS3)
#CNVs at the start of the gene 
CNV_rare3_LOSSatstart <- anti_join(CNV_rare3_LOSS2,CNV_rare3_LOSS3)
#CNVs at the end of the gene
CNV_rare3_LOSSatend <- anti_join(CNV_rare3_LOSS3,CNV_rare3_LOSS2)
CNV_rare3_LOSSatstart$match <- rep("start")
CNV_rare3_LOSSatend$match <- rep("end")
CNV_rare3_LOSStotal <- list(CNV_rare3_LOSS1,CNV_rare3_LOSSwithinCNV,CNV_rare3_LOSSatstart,CNV_rare3_LOSSatend) %>% reduce(full_join)

CNV_rare_LOSStotal <- full_join(CNV_rare1_LOSStotal, full_join(CNV_rare2_LOSStotal, CNV_rare3_LOSStotal)) #1640680
CNV_rare <- CNV_rare_LOSStotal[,c(1:3)]
CNV_rare <- unique(CNV_rare)#12809

write.table(CNV_rare_LOSStotal, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_rare_1640680CNV_12809patients_recheck.txt", sep="\t")
write.table(CNV_rare, file="/home/phalim/re_gecip/immune/phalim/SVCNV/List_of_12809rarepatients.txt", sep="\t")

####################################
####frequency####
CNV_PIDp_LOSStotal <- read.table(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_226709CNV_487patients.txt", sep="\t")
CNV_cancer_LOSStotal <- read.table(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_cancer_742543CNV_5783patients.txt", sep="\t")
CNV_rare_LOSStotal <- read.table(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_rare_1640680CNV_12809patients.txt", sep="\t")

CNV_PIDp_LOSStotal <- unique(CNV_PIDp_LOSStotal[,c(1,6,31)])
CNV_cancer_LOSStotal <- unique(CNV_cancer_LOSStotal[,c(1,2,17)])
CNV_rare_LOSStotal <- unique(CNV_rare_LOSStotal[,c(1,2,17)])

Freq_PID_immune_loss <- group_by(CNV_PIDp_LOSStotal, CNV_PIDp_LOSStotal$gene.name) %>% summarize(count=n())
write.table(Freq_PID_immune_loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDpatients_1pergene.txt", sep= "\t")
#PID patients  with loss CNV = 486

Freq_control_cancer_loss <- group_by(CNV_cancer_LOSStotal, CNV_cancer_LOSStotal$gene.name) %>% summarize(count=n())
write.table(Freq_control_cancer_loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_cancercontrol_1pergene.txt", sep= "\t")
Freq_control_rare_loss <- group_by(CNV_rare_LOSStotal, CNV_rare_LOSStotal$gene.name) %>% summarize(count=n())
write.table(Freq_control_rare_loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_rarecontrol_1pergene.txt", sep= "\t")
#control with loss CNV --> Cancer = 5783; rare = 12809

names(Freq_PID_immune_loss)[1:2] <- c("gene.name", "PID.loss")
names(Freq_control_cancer_loss)[1:2] <- c("gene.name", "cancer.control.loss")
names(Freq_control_rare_loss)[1:2] <- c("gene.name", "rare.control.loss")

Freq_PIDvsrare_Loss <- left_join(Freq_PID_immune_loss, Freq_control_rare_loss)
Freq_PIDvscancer_Loss <- left_join(Freq_PID_immune_loss, Freq_control_cancer_loss)

write.table(Freq_PIDvscancer_Loss,file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsCancer_1pergene.txt", sep= "\t")
write.table(Freq_PIDvsrare_Loss,file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsRare_1pergene.txt", sep= "\t")

names(Freq_PID_immune_loss)[2] <- c("count") #2784 genes
names(Freq_control_cancer_loss)[2]  <- c("count") #10169
names(Freq_control_rare_loss)[2]  <- c("count") #16818

Freq_PID_immune_loss$count <- (Freq_PID_immune_loss$count/486)*100
Freq_PID_immune_loss$Group <- rep("Loss CNV in PID Patients")
Freq_control_cancer_loss$count <- (Freq_control_cancer_loss$count/5783)*100
Freq_control_cancer_loss$Group <- rep("Loss CNV in Unrelated Cancer Patients")
Freq_control_rare_loss$count <- (Freq_control_rare_loss$count/12809)*100
Freq_control_rare_loss$Group <- rep("Loss CNV in Unrelated Rare Disease Patients")
Casevscancer_freq_Loss <- full_join(Freq_PID_immune_loss, Freq_control_cancer_loss)
Casevsrare_freq_Loss <- full_join(Freq_PID_immune_loss, Freq_control_rare_loss)

ggplot(Casevscancer_freq_Loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group))+  
  labs(color= "Group",title="Comparison of Loss CNVs",x="Genes", y="Frequency(%)" ) + 
  theme(legend.position = c(0.85,0.92),plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)), 
        legend.title=element_text(size=20), legend.spacing=unit(1,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=15, lineheight=5), legend.key.size=unit(1,"cm"))
  
ggplot(Casevsrare_freq_Loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group))+  
  labs(color= "Group",title="Comparison of Loss CNVs" ,x="Genes", y="Frequency(%)" ) + 
  theme(legend.position = c(0.85,0.92),plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)), 
        legend.title=element_text(size=20), legend.spacing=unit(1,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=15, lineheight=5), legend.key.size=unit(1,"cm"))

Freq_control_cancer_loss <- Freq_control_cancer_loss [,c(1:2)]
Freq_control_rare_loss <- Freq_control_rare_loss [,c(1:2)]
Freq_PID_immune_loss <- Freq_PID_immune_loss[,c(1:2)]
CasevsCANCER_freq_Loss <- left_join(Freq_PID_immune_loss, Freq_control_cancer_loss, by="gene.name")
names(CasevsCANCER_freq_Loss)[2:3] <- c("loss.CNV.PIDpatient","loss.CNV.cancerpatient")
CasevsRARE_freq_Loss <- left_join(Freq_PID_immune_loss, Freq_control_rare_loss, by="gene.name")
names(CasevsRARE_freq_Loss)[2:3] <- c("loss.CNV.PIDpatient","loss.CNV.rarediseasepatient")

write.table(CasevsCANCER_freq_Loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/Overlap_frequencies_lossCNV_PIDvscancer.txt", sep = "\t")
write.table(CasevsRARE_freq_Loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/Overlap_frequencies_lossCNV_PIDvsrare.txt", sep = "\t")

ggplot(CasevsCANCER_freq_Loss, aes(x=loss.CNV.PIDpatient, y=loss.CNV.cancerpatient))+ geom_point(aes(col=gene.name))+  
  geom_abline(col = "#C42126",size = 1)+
  theme(legend.position ="none" ,plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)))+ 
  labs(color= "Genes", title="Comparison of Loss CNVs" ,x="Loss CNV in PID Patients", y="Loss CNV in Unrelated Cancer Patients" )

ggplot(CasevsRARE_freq_Loss, aes(x=loss.CNV.PIDpatient, y=loss.CNV.rarediseasepatient))+ geom_point(aes(col=gene.name))+  
  geom_abline(col = "#C42126",size = 1)+
  theme(legend.position ="none" ,plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)))+ 
  labs(color= "Genes", title="Comparison of Loss CNVs" ,x="Loss CNV in PID Patients", y="Loss CNV in Unrelated Rare Disease Patients" )

#####Control vs control plot####
Freq_control_cancer_loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_cancercontrol_1pergene.txt", sep= "\t")
Freq_control_rare_loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_rarecontrol_1pergene.txt", sep= "\t")
names(Freq_control_cancer_loss)[1] <- c("gene.name")
names(Freq_control_rare_loss)[1] <- c("gene.name")
Freq_control_cancer_loss$count <- (Freq_control_cancer_loss$count/5783)*100
Freq_control_cancer_loss$Group <- rep("Loss CNV in Unrelated Cancer Patients")
Freq_control_rare_loss$count <- (Freq_control_rare_loss$count/12809)*100
Freq_control_rare_loss$Group <- rep("Loss CNV in Unrelated Rare Disease Patients")
Cancervsrare_freq_Loss <- full_join(Freq_control_cancer_loss, Freq_control_rare_loss)

ggplot(Cancervsrare_freq_Loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group))+  
  labs(color= "Group",title="Comparison of Loss CNVs" ,x="Genes", y="Frequency(%)" ) + 
  theme(legend.position = c(0.85,0.92),plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)), 
        legend.title=element_text(size=20), legend.spacing=unit(1,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=15, lineheight=5), legend.key.size=unit(1,"cm"))

Freq_control_cancer_loss <- Freq_control_cancer_loss [,c(1:2)]
Freq_control_rare_loss <- Freq_control_rare_loss [,c(1:2)]
names(Freq_control_cancer_loss)[1:2] <- c("gene.name", "cancer.control.loss")
names(Freq_control_rare_loss)[1:2] <- c("gene.name", "rare.control.loss")
CancervsRare_freq_Loss <- full_join(Freq_control_cancer_loss, Freq_control_rare_loss, by="gene.name")
CancervsRare_freq_Loss <- CancervsRare_freq_Loss[complete.cases(CancervsRare_freq_Loss),]
write.table(CancervsRare_freq_Loss,file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_percentage_CancervsRare_1pergene.txt", sep= "\t")

ggplot(CancervsRare_freq_Loss, aes(x=cancer.control.loss, y=rare.control.loss))+ geom_point(aes(col=gene.name))+  
  geom_abline(col = "#C42126",size = 1)+
  theme(legend.position ="none" ,plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)))+ 
  labs(color= "Genes", title="Comparison of Loss CNVs" ,x="Loss CNV in Unrelated Cancer Patients", y="Loss CNV in Unrelated Rare Disease Patients" )
--------
#log10 graphs per Group
Freq_PID_immune_loss$count <- log10(Freq_PID_immune_loss$count)
Freq_control_cancer_loss$count <- log10(Freq_control_cancer_loss$count)
Freq_control_rare_loss$count <- log10(Freq_control_rare_loss$count)
Freq_PID_immune_loss$Group <- rep("Loss CNV in PID Patients")
Casevscancer_freq_Loss <- full_join(Freq_PID_immune_loss, Freq_control_cancer_loss)
Casevsrare_freq_Loss <- full_join(Freq_PID_immune_loss, Freq_control_rare_loss)
Freq_control_cancer_loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_cancercontrol_1pergene.txt", sep= "\t")
Freq_control_rare_loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_rarecontrol_1pergene.txt", sep= "\t")
names(Freq_control_cancer_loss)[1] <- c("gene.name")
names(Freq_control_rare_loss)[1] <- c("gene.name")
Freq_control_cancer_loss$count <- (Freq_control_cancer_loss$count/5783)*100
Freq_control_cancer_loss$Group <- rep("Loss CNV in Unrelated Cancer Patients")
Freq_control_rare_loss$count <- (Freq_control_rare_loss$count/12809)*100
Freq_control_rare_loss$Group <- rep("Loss CNV in Unrelated Rare Disease Patients")
Cancervsrare_freq_Loss <- full_join(Freq_control_cancer_loss, Freq_control_rare_loss)
Cancervsrare_freq_Loss <- CancervsRare_freq_Loss[complete.cases(CancervsRare_freq_Loss),]

ggplot(Casevscancer_freq_Loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group))+  
  labs(color= "Group",title="Comparison of Loss CNVs" ,x="Genes", y="Frequency(%)" ) + 
  theme(legend.position = c(0.85,0.92),plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)), 
        legend.title=element_text(size=20), legend.spacing=unit(1,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=15, lineheight=5), legend.key.size=unit(1,"cm"))
ggplot(Casevsrare_freq_Loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group))+  
  labs(color= "Group",title="Comparison of Loss CNVs" ,x="Genes", y="Frequency(%)" ) + 
  theme(legend.position = c(0.85,0.92),plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)), 
        legend.title=element_text(size=20), legend.spacing=unit(1,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=15, lineheight=5), legend.key.size=unit(1,"cm"))
ggplot(Cancervsrare_freq_Loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group))+  
  labs(color= "Group",title="Comparison of Loss CNVs" ,x="Genes", y="Frequency(%)" ) + 
  theme(legend.position = c(0.85,0.92),plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)), 
        legend.title=element_text(size=20), legend.spacing=unit(1,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=15, lineheight=5), legend.key.size=unit(1,"cm"))

#log10 graphs per gene
Freq_control_cancer_loss <- Freq_control_cancer_loss [,c(1:2)]
Freq_control_rare_loss <- Freq_control_rare_loss [,c(1:2)]
Freq_PID_immune_loss <- Freq_PID_immune_loss[,c(1:2)]
CasevsCANCER_freq_Loss <- left_join(Freq_PID_immune_loss, Freq_control_cancer_loss, by="gene.name")
names(CasevsCANCER_freq_Loss)[2:3] <- c("loss.CNV.PIDpatient","loss.CNV.cancerpatient")
CasevsRARE_freq_Loss <- left_join(Freq_PID_immune_loss, Freq_control_rare_loss, by="gene.name")
names(CasevsRARE_freq_Loss)[2:3] <- c("loss.CNV.PIDpatient","loss.CNV.rarediseasepatient")

write.table(CasevsCANCER_freq_Loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/Overlap_frequencies_lossCNV_PIDvscancer_log10.txt", sep = "\t")
write.table(CasevsRARE_freq_Loss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/Overlap_frequencies_lossCNV_PIDvsrare_log10.txt", sep = "\t")
write.table(CancervsRare_freq_Loss,file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_percentage_CancervsRare_log10.txt", sep= "\t")

ggplot(CasevsCANCER_freq_Loss, aes(x=loss.CNV.PIDpatient, y=loss.CNV.cancerpatient))+ geom_point(aes(col=gene.name))+  
  geom_abline(col = "#C42126",size = 1)+
  theme(legend.position ="none" ,plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)))+ 
  labs(color= "Genes", title="Comparison of Loss CNVs" ,x="Loss CNV in PID Patients", y="Loss CNV in Unrelated Cancer Disease Patients" )
ggplot(CasevsRARE_freq_Loss, aes(x=loss.CNV.PIDpatient, y=loss.CNV.rarediseasepatient))+ geom_point(aes(col=gene.name))+  
  geom_abline(col = "#C42126",size = 1)+
  theme(legend.position ="none" ,plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)))+ 
  labs(color= "Genes", title="Comparison of Loss CNVs" ,x="Loss CNV in PID Patients", y="Loss CNV in Unrelated Rare Disease Patients" )
ggplot(CancervsRare_freq_Loss, aes(x=count.x, y=count.y))+ geom_point(aes(col=gene.name))+  
  geom_abline(col = "#C42126",size = 1)+
  theme(legend.position ="none" ,plot.margin=margin(20,20,20,20),
        plot.title=element_text(size=35, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=25, margin = margin(30)),
        axis.title.y=element_text(size=25, margin = margin(0,30)))+ 
  labs(color= "Genes", title="Comparison of Loss CNVs" ,x="Loss CNV in Unrelated Cancer Patients", y="Loss CNV in Unrelated Rare Disease Patients" )

#######READ DEPTH BOXPLOT######## 25/05/2020
CNV_PIDp <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_226709CNV_487patients.txt", sep="\t")
CNV_rare <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_rare_1640680CNV_12809patients.txt", sep="\t")
CNV_cancer <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_cancer_742543CNV_5783patients.txt", sep="\t")

CNV_PIDp <- CNV_PIDp[,c(6,27)]
CNV_rare <- CNV_rare[,c(2,13)]
CNV_cancer <- CNV_cancer[,c(2,13)]

names(CNV_PIDp)[2] <- c("RD_PID")
names(CNV_rare)[2] <- c("RD_rare")
names(CNV_cancer)[2] <- c("RD_cancer")

RD_list <- list(CNV_PIDp, CNV_cancer, CNV_rare) %>% reduce(full_join)

op <- par(mar = c(6,6,6,3) + 0.1)
boxplot(RD_list$RD_PID,RD_list$RD_cancer,RD_list$RD_rare,
        main = "Read Depth Comparison", cex.main = 2.6,
        at = c(1,2,3),
        las = 2,
        col = c("red","green","blue"),
        border = "brown",
        axes = FALSE, ann = FALSE)
axis(1, at = 1:3, labels = c("PID", "Cancer", "Rare Disease"), cex.axis = 1.7)
axis(2, cex.axis = 1.7)
title(xlab = "Groups", cex.lab = 2.1,
      line = 3)
title(ylab = "Read Depth Units", cex.lab = 2.1,
      line = 3)
box()
par(op)

#####Freq per gene HISTOGRAM######
Freq_PIDvscancer_Loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsCancer.txt", sep= "\t")
Freq_PIDvsrare_Loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsRare.txt", sep= "\t")
Freq_PID_CANCER_RARE <- full_join(Freq_PIDvscancer_Loss,Freq_PIDvsrare_Loss)
names(Freq_PID_CANCER_RARE)[2:4] <- c("PID","Cancer","Rare Disease")
Freq_PID_CANCER_RARE$PID <- (Freq_PID_CANCER_RARE$PID/486*100)
Freq_PID_CANCER_RARE$Cancer <- (Freq_PID_CANCER_RARE$PID/5783*100)
Freq_PID_CANCER_RARE$`Rare Disease` <- (Freq_PID_CANCER_RARE$PID/12809*100)
dataset <- Freq_PID_CANCER_RARE %>% gather(Groups, Frequency, 2:4)

dataset %>%
  arrange(dataset$Frequency) %>%
  mutate(name=factor(gene.name, levels=rev(dataset$gene.name),ordered = TRUE))%>%
  ggplot(aes(x =gene.name, y =Frequency)) + geom_bar(aes(fill = Groups, y=..density..),position="identity") + 
  labs(title="Comparison of CNV Frequencies in All Genes" ,x="Genes", y="Frequency" )+
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,1))+
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=45),
        axis.title.x=element_text(size=55, margin = margin(50)),
        axis.title.y=element_text(size=55, margin = margin(0,40)), 
        plot.title=element_text(size=75, hjust=0.5, margin=margin(30,30,30,30)),
        legend.title=element_text(size=35), legend.spacing=unit(2,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=30, lineheight=10), legend.key.size=unit(1.5,"cm"),
        legend.position = c(0.94,0.92), plot.margin=margin(45,45,45,45))+
  facet_wrap(~Groups, nrow = 3, ncol = 1, scales = "fixed",
             shrink = TRUE, as.table = TRUE,
             switch = NULL, drop = TRUE, dir = "h", strip.position = "top")+
  theme(strip.text.x = element_text(size = 40, margin = margin(20,20,20,20)), panel.spacing = unit(1, "cm"))


#####Freq BOXPLOT######
Freq_PIDvscancer_Loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsCancer_1pergene.txt", sep= "\t")
Freq_PIDvsrare_Loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsRare_1pergene.txt", sep= "\t")
Freq_PID_CANCER_RARE <- full_join(Freq_PIDvscancer_Loss,Freq_PIDvsrare_Loss)
names(Freq_PID_CANCER_RARE)[2:4] <- c("PID","Cancer","Rare Disease")

op <- par(mar = c(6,12,6,3) + 0.5,las = 1)
boxplot(Freq_PID_CANCER_RARE$PID,Freq_PID_CANCER_RARE$Cancer,Freq_PID_CANCER_RARE$`Rare Disease`,
        main = "Frequency Comparison", cex.main = 4.1,
        at = c(1,2,3),
        col = c("red","green","blue"),
        border = "brown",
        axes = FALSE, ann = FALSE)
axis(1, at = 1:3, labels = c("PID", "Cancer", "Rare Disease"),cex.axis = 2.7, mgp= c(3,2,0))
axis(2, cex.axis = 2.7)
title(xlab = "Groups", cex.lab = 3.2,line = 5)
title(ylab = "Frequency", cex.lab = 3.2,line = 8)
box()
par(op)

#####Freq per sample HISTOGRAM######
CNV_PIDp_LOSStotal <- read.table(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_226709CNV_487patients.txt", sep="\t")
CNV_cancer_LOSStotal <- read.table(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_cancer_742543CNV_5783patients.txt", sep="\t")
CNV_rare_LOSStotal <- read.table(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_rare_1640680CNV_12809patients.txt", sep="\t")

CNV_PIDp_LOSStotal <- unique(CNV_PIDp_LOSStotal[,c(1,6,31)])
CNV_cancer_LOSStotal <- unique(CNV_cancer_LOSStotal[,c(1,2,17)])
CNV_rare_LOSStotal <- unique(CNV_rare_LOSStotal[,c(1,2,17)])

Freq_PID_immune_sample <- group_by(CNV_PIDp_LOSStotal, CNV_PIDp_LOSStotal$Plate.key) %>% summarize(count=n())
write.table(Freq_PID_immune_sample, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_CNV_persample_PID_1pergene.txt", sep= "\t")
Freq_control_cancer_sample <- group_by(CNV_cancer_LOSStotal, CNV_cancer_LOSStotal$Plate.key) %>% summarize(count=n())
write.table(Freq_control_cancer_sample, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_CNV_persample_cancer_1pergene.txt", sep= "\t")
Freq_control_rare_sample <- group_by(CNV_rare_LOSStotal, CNV_rare_LOSStotal$Plate.key) %>% summarize(count=n())
write.table(Freq_control_rare_sample, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_CNV_persample_rare_1pergene.txt", sep= "\t")
names(Freq_PID_immune_sample)[1:2] <- c("Plate.key","PID")
names(Freq_control_cancer_sample)[1:2] <- c("Plate.key","Cancer")
names(Freq_control_rare_sample)[1:2] <- c("Plate.key","Rare Disease")

Freq_PID_CANCER_RARE <- full_join(Freq_PID_immune_sample,full_join(Freq_control_cancer_sample,Freq_control_rare_sample))
dataset <- Freq_PID_CANCER_RARE %>% gather(Groups, Frequency, 2:4)
dataset <- dataset[complete.cases(dataset),]
#Freqgroup <-group_by(dataset, dataset$Groups) %>% summarize(count=n())
#names(Freqgroup)[1] <-"Groups"
#Dataset <- full_join(dataset, Freqgroup)
#Dataset$count <-Dataset$Frequency/Dataset$count*100 
#Dataset$count <- paste0(Dataset$count,"%")
Dat <- split.data.frame(dataset, dataset$Groups)
pid <- as.data.frame(Dat[[2]])
cancer <- as.data.frame(Dat[[1]])
rare <- as.data.frame(Dat[[3]])
Fixed <- rbind(pid,cancer,rare)
#combined plot + remove density
ggplot(Fixed, aes(x=Frequency,color=Groups, fill=Groups)) + 
  geom_histogram(aes(y=..density..),position="identity",alpha=0.7, bins=150)+
  labs(title = "Comparison of Genes with CNVs", x="Number of Genes with CNVs/Patient", 
       y="Number of Patient per Group (%)")+
  #geom_density(alpha=0.7)+ geom_vline(aes(xintercept = mean(Frequency),color=Groups),linetype="dashed")+
  scale_fill_manual(values =c("lightseagreen","purple1","orchid2"))+
  scale_color_manual(values =c("grey50","grey50","grey50"))+
    scale_y_continuous(labels = scales::percent_format(),limits = c(0,0.02))+
  scale_x_continuous(breaks = seq(0,600, by=50))+
  theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40),
        plot.title=element_text(size=75, hjust=0.5, margin=margin(20,20,20,20)),
        axis.title.x=element_text(size=50, margin = margin(30)),
        axis.title.y=element_text(size=50, margin = margin(0,30)), 
        legend.title=element_text(size=40), legend.spacing=unit(2,"cm"),
        legend.margin = margin(15,15,15,15),
        legend.text=element_text(size=33, lineheight=12), legend.key.size=unit(2,"cm"),
        legend.position = c(0.93,0.86), plot.margin=margin(50,45,45,45))
dev.off()
#grid.arranged plot
p <- ggplot(pid, aes(x=Frequency,color=Groups, fill=Groups)) + 
  geom_histogram(aes(y=..density..),position="identity",alpha=0.7, bins=150)+
  labs(title="Comparison of Genes with CNVs",x="Number of Genes with CNVs/Patient", y="Number of Patient per Group (%)")+
  geom_density(alpha=0.7)+ geom_vline(aes(xintercept = mean(Frequency),color=Groups),linetype="dashed")+
  scale_fill_manual(values =c("purple1"))+
  scale_color_manual(values =c("grey50","grey50","grey50"))+
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,0.02))+
  scale_x_continuous(breaks = seq(0,600, by=50))+
  theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40),
        plot.title=element_text(size=75, hjust=0.5, margin=margin(20,20,20,20)),
        axis.title.x=element_text(size=50, margin = margin(30)),
        axis.title.y=element_text(size=50, margin = margin(0,30)), 
        legend.title=element_text(size=40), legend.spacing=unit(2,"cm"),
        legend.margin = margin(15,15,15,15),
        legend.text=element_text(size=33, lineheight=12), legend.key.size=unit(2,"cm"),
        legend.position = c(0.95,0.86), plot.margin=margin(50,45,45,45))

c <- ggplot(cancer, aes(x=Frequency,color=Groups, fill=Groups)) + 
  geom_histogram(aes(y=..density..),position="identity",alpha=0.7, bins=150)+
  labs(x="Number of Genes with CNVs/Patient", y="Number of Patient per Group (%)")+
  geom_density(alpha=0.7)+ geom_vline(aes(xintercept = mean(Frequency),color=Groups),linetype="dashed")+
  scale_fill_manual(values =c("lightseagreen"))+
  scale_color_manual(values =c("grey50","grey50","grey50"))+
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,0.02))+
  scale_x_continuous(breaks = seq(0,600, by=50))+
  theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40),
        plot.title=element_text(size=75, hjust=0.5, margin=margin(20,20,20,20)),
        axis.title.x=element_text(size=50, margin = margin(30)),
        axis.title.y=element_text(size=50, margin = margin(0,30)), 
        legend.title=element_text(size=40), legend.spacing=unit(2,"cm"),
        legend.margin = margin(15,15,15,15),
        legend.text=element_text(size=33, lineheight=12), legend.key.size=unit(2,"cm"),
        legend.position = c(0.95,0.86), plot.margin=margin(50,45,45,45))
r <- ggplot(rare, aes(x=Frequency,color=Groups, fill=Groups)) + 
  geom_histogram(aes(y=..density..),position="identity",alpha=0.7, bins=150)+
  labs(x="Number of Genes with CNVs/Patient",y="Number of Patient per Group (%)")+
  geom_density(alpha=0.7)+ geom_vline(aes(xintercept = mean(Frequency),color=Groups),linetype="dashed")+
  scale_fill_manual(values =c("orchid2"))+
  scale_color_manual(values =c("grey50","grey50","grey50"))+
  scale_y_continuous(labels = scales::percent_format(),limits = c(0,0.02))+
  scale_x_continuous(breaks = seq(0,600, by=50))+
  theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40),
        plot.title=element_text(size=75, hjust=0.5, margin=margin(20,20,20,20)),
        axis.title.x=element_text(size=50, margin = margin(30)),
        axis.title.y=element_text(size=50, margin = margin(0,30)), 
        legend.title=element_text(size=40), legend.spacing=unit(2,"cm"),
        legend.margin = margin(15,15,15,15),
        legend.text=element_text(size=33, lineheight=12), legend.key.size=unit(2,"cm"),
        legend.position = c(0.93,0.86), plot.margin=margin(50,45,45,45))

grid.arrange(p,c,r,
             layout_matrix=
               rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                     c(NA,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                       2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                       2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                       2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                       2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                     c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                       3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                       3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                       3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                       3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)))

#Cases that is specific to PID patients
PIDvsCancer <- read.csv("SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsCancer_1pergene.txt", sep="\t")
PIDvsrare <- read.csv("SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsRare_1pergene.txt", sep="\t")
PIDvsCancer[is.na(PIDvsCancer)] <- 0
PIDvsrare[is.na(PIDvsrare)] <- 0
PCsplit <- split.data.frame(PIDvsCancer, PIDvsCancer$cancer.control.loss)
PRsplit <- split.data.frame(PIDvsrare,PIDvsrare$rare.control.loss)
PC_0 <- as.data.frame(PCsplit[[1]])
PR_0 <- as.data.frame(PRsplit[[1]])
PCR_0 <- full_join(PC_0, PR_0)
PCR_0 <- PCR_0[complete.cases(PCR_0),]
PIDcases <- read.csv("SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_226709CNV_487patients.txt", sep="\t")
PIDcases <- PIDcases[,c(1:3,6:8,20:22,26,27,29:31)]
PIDcaseonly <- full_join(PIDcases, PCR_0)
PIDcaseonly <- unique(PIDcaseonly[complete.cases(PIDcaseonly),]) 

library(data.table)
PIDcaseonly <- setDT(PIDcaseonly)[, Recruited.Disease.Or.Panel := paste0(as.character(Recruited.Disease.Or.Panel), collapse = ", "), by = Plate.key]
PIDcaseonly <- unique(PIDcaseonly[,c(1:15)]) #301/2 CNVs from 262 genes in list from 74 patients
write.table(PIDcaseonly, file="SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_302CNVs_caseonly_1pergene.txt", sep="\t")

# read depth vs length
CNV_PIDp <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_PID_226709CNV_487patients.txt", sep="\t")
CNV_rare <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_rare_1640680CNV_12809patients.txt", sep="\t")
CNV_cancer <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_Allgenes_comparison/CNV_Loss_cancer_742543CNV_5783patients.txt", sep="\t")

PIDpatientloss_n <- CNV_PIDp[,c(25,27)]
CNV.loss.cancer_n <- CNV_cancer[,c(11,13)]
CNV.loss.rare_n <- CNV_rare[,c(11,13)]

PIDpatientloss_n$Group <- rep("PID")
CNV.loss.cancer_n$Group <- rep("Cancer")
CNV.loss.rare_n$Group <- rep("Rare Disease")

LengthvsRD <- rbind(PIDpatientloss_n,CNV.loss.cancer_n,CNV.loss.rare_n)

ggplot(LengthvsRD, aes(x=dbl.1,color=Group, fill=Group)) + 
  geom_point(aes(y=length.loss),position="identity",size=3, alpha=0.5)+
  labs(title = "Comparison of Read Depth and CNV calls length", x="Read Depth", 
       y="Length of CNV calls")+
  #geom_density(alpha=0.7)+ geom_vline(aes(xintercept = mean(Frequency),color=Groups),linetype="dashed")+
  scale_fill_manual(values =c("lightseagreen","purple1","orchid2"))+
  scale_color_manual(values =c("lightseagreen","purple1","orchid2"))+
  scale_y_continuous(limits = c(0, 14000000))+
  scale_x_continuous(breaks = seq(10,40, by=1))+
  theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40),
        plot.title=element_text(size=75, hjust=0.5, margin=margin(20,20,20,20)),
        axis.title.x=element_text(size=50, margin = margin(30)),
        axis.title.y=element_text(size=50, margin = margin(0,30)), 
        legend.title=element_text(size=40), legend.spacing=unit(2,"cm"),
        legend.margin = margin(15,15,15,15),
        legend.text=element_text(size=33, lineheight=12), legend.key.size=unit(2,"cm"),
        legend.position = c(0.93,0.86), plot.margin=margin(50,45,45,45))

# length density
ggplot(LengthvsRD, aes(x=length.loss,color=Group, fill=Group)) + 
  geom_density(aes(y=..density..),position="identity",size=1, alpha=0.4)+
  labs(title = "Comparison of CN loss calls' length", x="CN loss length", 
       y="density")+
  #geom_density(alpha=0.7)+ geom_vline(aes(xintercept = mean(Frequency),color=Groups),linetype="dashed")+
  scale_fill_manual(values =c("lightseagreen","purple1","orchid2"))+
  scale_color_manual(values =c("lightseagreen","purple1","orchid2"))+
  scale_x_continuous(limits = c(0, 14000000))+
  #scale_x_continuous(breaks = seq(10,40, by=1))+
  theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40),
        plot.title=element_text(size=75, hjust=0.5, margin=margin(20,20,20,20)),
        axis.title.x=element_text(size=50, margin = margin(30)),
        axis.title.y=element_text(size=50, margin = margin(0,30)), 
        legend.title=element_text(size=40), legend.spacing=unit(2,"cm"),
        legend.margin = margin(15,15,15,15),
        legend.text=element_text(size=33, lineheight=12), legend.key.size=unit(2,"cm"),
        legend.position = c(0.94,0.86),plot.margin=margin(50,45,45,45))

#Freq per gene PID vs controls, control vs control + Sudmant
Freq_sudmantfull<- read.csv(file = "/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_1KGP_vs_case_control/Sudmant_CNVLoss_allgene_freq_1pergene_complete.csv", sep="\t")
Freq_sudmantfull <- Freq_sudmantfull[,c(1,7)]
names(Freq_sudmantfull)[1:2] <- c("gene.name", "sudmant.loss")
Freq_PIDvscancer_Loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsCancer_1pergene.txt", sep= "\t")
Freq_PIDvsrare_Loss <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_frequency_statistics/Count_allgenes_lossCNV_PIDvsRare_1pergene.txt", sep= "\t")
Freq_casevscontrols <- list(Freq_PIDvsrare_Loss,Freq_PIDvscancer_Loss,
                            Freq_sudmantfull)%>% reduce(full_join)
names(Freq_casevscontrols)[3:4] <- c("rare.loss","cancer.loss")
Freq_casevscontrols$PID.loss <- (Freq_casevscontrols$PID.loss/486*100)
Freq_casevscontrols$cancer.loss <- (Freq_casevscontrols$cancer.loss/5783*100)
Freq_casevscontrols$rare.loss <- (Freq_casevscontrols$rare.loss/12809*100)
Freq_casevscontrols$sudmant.loss <- (Freq_casevscontrols$sudmant.loss/2504*100)

ggplot(Freq_casevscontrols, aes(x=PID.loss, y=cancer.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in PID Group(%)", 
       y="CN Loss in Unrelated Cancer Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.position = "none")

ggplot(Freq_casevscontrols, aes(x=PID.loss, y=rare.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in PID Group(%)", 
       y="CN Loss in Unrelated Rare Disease Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.position = "none")

ggplot(Freq_casevscontrols, aes(x=PID.loss, y=sudmant.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in PID Group(%)", 
       y="CN Loss in Sudmant Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.position = "none")

ggplot(Freq_casevscontrols, aes(x=rare.loss, y=cancer.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in Unrelated Rare Disease Group(%)", 
       y="CN Loss in Unrelated Cancer Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.position = "none")

ggplot(Freq_casevscontrols, aes(x=rare.loss, y=sudmant.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in Unrelated Rare Disease Group(%)", 
       y="CN Loss in Sudmant Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.position = "none")

ggplot(Freq_casevscontrols, aes(x=cancer.loss, y=sudmant.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in Unrelated Cancer Group(%)", 
       y="CN Loss in Sudmant Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.position = "none")

PID <- Freq_casevscontrols[,c(1,2)]
PID$Group <- rep("PID")
names(PID)[2] <-"count"
RARE <- Freq_casevscontrols[,c(1,3)]
RARE$Group <- rep("Rare Disease")
names(RARE)[2] <-"count"
CANCER <- Freq_casevscontrols[,c(1,4)]
CANCER$Group <- rep("Cancer")
names(CANCER)[2] <-"count"
SUDMANT <- Freq_casevscontrols[,c(1,5)]
SUDMANT$Group <- rep("Sudmant")
names(SUDMANT)[2] <-"count"
Case_Control_freq_loss <- list(PID,RARE,CANCER,SUDMANT)%>% reduce(full_join)

ggplot(Case_Control_freq_loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group),size=2)+  
  labs(color= "Group", title="Comparison of CN Loss Proportion" ,x="Gene Name", y="Frequency (%)" ) + 
  theme(legend.position = c(0.94,0.82),plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_blank(), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)), 
        legend.title=element_text(size=17), legend.spacing=unit(0.8,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=13, lineheight=5), legend.key.size=unit(0.8,"cm"))
