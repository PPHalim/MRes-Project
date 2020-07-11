library(Rlabkey)
library(tidyverse, lib="~/re_gecip/immune/phalim/R_packages")
library(curl)
library(rjson)
library(RCurl)
library(dplyr)
library(purrr)
library(ggplot2)

#setting working directory
setwd("/home/phalim/re_gecip/immune/phalim/")

#open control data
ImmuneGene_CNV_cancertotal <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_2537cancersample.txt", sep= "\t")
ImmuneGene_CNV_raretotal <- read.csv(file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/List_MatchedCNVimmunegene_control_5786raresample.txt", sep= "\t")

ImmuneGene_CNV_cancertotal$Chromosome <- as.character(ImmuneGene_CNV_cancertotal$Chromosome)
ImmuneGene_CNV_raretotal$Chromosome <- as.character(ImmuneGene_CNV_raretotal$Chromosome)

ImmuneGene_CNV_cancertotal <- ImmuneGene_CNV_cancertotal[,c(1,2,3,6,15,14)]
ImmuneGene_CNV_raretotal <- ImmuneGene_CNV_raretotal[,c(1,2,3,6,15,14)]

#Separate control data into gain and loss
ImmuneGene_CNV_cancer_list <- split(ImmuneGene_CNV_cancertotal, ImmuneGene_CNV_cancertotal$condition)
ImmuneGene_CNV_rare_list <- split(ImmuneGene_CNV_raretotal, ImmuneGene_CNV_raretotal$condition)

ImmuneGene_CNV_cancer_loss <- as.data.frame(ImmuneGene_CNV_cancer_list[[2]])
ImmuneGene_CNV_rare_loss <- as.data.frame(ImmuneGene_CNV_rare_list[[2]])

#open immune patient data, remove patient disease info+unique to simplify
ImmuneGene_CNV_PID_all <- read.csv("/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_PID_immune_gene/List_MatchedCNV_immunegene.235patient.txt", sep= " ")
ImmuneGene_CNV_PID_list <- split(ImmuneGene_CNV_PID_all, ImmuneGene_CNV_PID_all$condition)
ImmuneGene_CNV_PID_loss <- as.data.frame(ImmuneGene_CNV_PID_list[[2]])
ImmuneGene_CNV_PID_loss <- unique(ImmuneGene_CNV_PID_loss[,c(1,2,3,6,13,12)])

ImmuneGene_CNV_PID_loss$Chromosome <- as.character(ImmuneGene_CNV_PID_loss$Chromosome)
Freq_sudmantPID <- read.csv(file = "/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_1KGP_vs_case_control/Sudmant_CNVLoss_PIDgene_freq_1pergene_complete.csv", sep="\t")
Freq_sudmantPID <- Freq_sudmantPID[,c(1,7)]
names(Freq_sudmantPID)[1:2] <- c("gene.name", "count")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#09/03/2020
#count frequencies of genes overlapped in immune patient: gain. loss, total
Freq_PID_immune_loss <- group_by(ImmuneGene_CNV_PID_loss, ImmuneGene_CNV_PID_loss$gene.name) %>% summarize(count=n())
Freq_PID_immune_loss$count <- (Freq_PID_immune_loss$count/164)*100
names(Freq_PID_immune_loss)[1] <- c("gene.name")  
Freq_PID_immune_loss$Group <- rep("PID")

#count frequencies of genes overlapped in control: rare, cancer, total (also gain, loss, total)
Freq_control_cancer_loss <- group_by(ImmuneGene_CNV_cancer_loss, ImmuneGene_CNV_cancer_loss$gene.name) %>% summarize(count=n())
names(Freq_control_cancer_loss)[1] <- c("gene.name")
Freq_control_cancer_loss$count <- (Freq_control_cancer_loss$count/1754)*100
Freq_control_cancer_loss$Group <- rep("Cancer")

Freq_control_rare_loss <- group_by(ImmuneGene_CNV_rare_loss, ImmuneGene_CNV_rare_loss$gene.name) %>% summarize(count=n())
names(Freq_control_rare_loss)[1] <- c("gene.name")
Freq_control_rare_loss$count <- (Freq_control_rare_loss$count/4105)*100
Freq_control_rare_loss$Group <- rep("Rare Disease")

Freq_sudmantPID$count <- (Freq_sudmantPID$count/2504)*100
Freq_sudmantPID$Group <- rep("Sudmant")

#combine data control freq -LOSS
Control_PID_freq_loss <- list(Freq_control_cancer_loss,Freq_control_rare_loss)%>% reduce(full_join)
Loss_Gene_CASE <- Freq_PID_immune_loss$gene.name
Control_PID_freq_loss <-Control_PID_freq_loss[Control_PID_freq_loss$gene.name %in% Loss_Gene_CASE ,]
Case_Control_freq_loss <- full_join(Freq_PID_immune_loss, Control_PID_freq_loss)

#18 GENES
write.table(Control_PID_freq_loss, file= "/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/Freq_percentage_loss_18genes_3controlgroups.txt", sep= "\t")
write.table(Case_Control_freq_loss, file= "/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/Freq_percentage_loss_18genes_4groups.txt", sep= "\t")

#another run to get the frequency (not percentage)
write.table(Case_Control_freq_loss, file= "/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_Control_unrelated/Freq_loss_18genes_4groups.txt", sep= "\t")
Case_Control_freq_loss[is.na(Case_Control_freq_loss)]<- "0"
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#02/07/2020

#Plot for loss CNV per category
ggplot(Case_Control_freq_loss, aes(x=gene.name, y=count))+ geom_point(aes(col=Group),size=2)+  
labs(color= "Group", title="Comparison of CN Loss Proportion" ,x="Gene Name", y="Frequency (%)" ) + 
  theme(legend.position = c(0.93,0.85),plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)), 
        legend.title=element_text(size=17), legend.spacing=unit(0.8,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=13, lineheight=5), legend.key.size=unit(0.8,"cm"))
  #saved

#:::::::::::::::::make plot for less common case only with total control only ::::::::::::::::::::::::::::::::::#
#02/07/2020

names(Freq_PID_immune_loss)[2] <-"PID.loss"
names(Freq_control_cancer_loss)[2]<-"cancer.loss"
names(Freq_control_rare_loss)[2] <-"rare.loss"
names(Freq_sudmantPID)[2] <- "sudmant.loss"
Freq_PID_immune_loss<-Freq_PID_immune_loss[,c(1,2)] 
Freq_control_cancer_loss<-Freq_control_cancer_loss[,c(1,2)] 
Freq_control_rare_loss<-Freq_control_rare_loss[,c(1,2)] 
Freq_sudmantPID<-Freq_sudmantPID[,c(1,2)]  
Freq_casevscontrols <-list(Freq_PID_immune_loss,Freq_control_cancer_loss,
                           Freq_control_rare_loss,Freq_sudmantPID)%>% reduce(full_join)
Freq_casevscontrols <- Freq_casevscontrols[Freq_casevscontrols$gene.name %in% Loss_Gene_CASE ,]

ggplot(Freq_casevscontrols, aes(x=PID.loss, y=cancer.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in PID Group(%)", 
       y="CN Loss in Unrelated Cancer Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.title=element_text(size=17), legend.spacing=unit(0.8,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=12, lineheight=5), legend.key.size=unit(0.8,"cm"))

ggplot(Freq_casevscontrols, aes(x=PID.loss, y=rare.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in PID Group(%)", 
       y="CN Loss in Unrelated Rare Disease Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.title=element_text(size=17), legend.spacing=unit(0.8,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=12, lineheight=5), legend.key.size=unit(0.8,"cm"))

ggplot(Freq_casevscontrols, aes(x=PID.loss, y=sudmant.loss))+ geom_point(aes(col=gene.name),size=2)+  
  geom_abline(col = "#C42126",size = 1)+
  labs(color= "Genes", title="Comparison of CN Loss per Gene" ,x="CN Loss in PID Group(%)", 
       y="CN Loss in Sudmant Group(%)" )+
  theme(plot.margin=margin(10,10,10,10),
        plot.title=element_text(size=25, hjust=0.5, margin=margin(10,10,10,10)),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20, margin = margin(20)),
        axis.title.y=element_text(size=20, margin = margin(0,20)),
        legend.title=element_text(size=17), legend.spacing=unit(0.8,"cm"),
        legend.margin = margin(10,10,10,10),
        legend.text=element_text(size=12, lineheight=5), legend.key.size=unit(0.8,"cm"))
