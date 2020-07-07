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

#gene information downloaded from ensembl.com
Gene_information <- read.csv("SVCNV/Geneinfo_ensembl.txt")
GeneRangeDF <- Gene_information[,c(1:4)]
names(GeneRangeDF)[1:4] <- c("Chromosome","start","end","gene.name")
GeneRangeDFa <- unique(GeneRangeDF)

#PID gene list
PIDgene1 <- read.csv("PanelApp_x389_PID_genesymbols.txt") #388 genes
PIDgene2 <- read.csv("PanelAlt_x87_plus_x33_PID_genesymbols.txt") #119 genes

names(PIDgene1)[1] <-"gene.name" 
names(PIDgene2)[1] <-"gene.name"
PIDgenefull <- full_join(PIDgene1, PIDgene2)#507 genes

#PID gene range
PIDgenerange <- inner_join(GeneRangeDFa, PIDgenefull, by="gene.name") #596 gene ranges
GeneRangeonly <- makeGRangesFromDataFrame(PIDgenerange, keep.extra.columns = TRUE)#must use library(GenomicRanges), do not work for this dataset

write.table(PIDgenerange, file = "SVCNV/CNV_in_PID_immune_gene/PID_genes_range.txt", sep = "\t")
#CNV in immune patient found from processing in SVCNVscript.R
ImmuneCNVdata <- read.csv("SVCNV/CNV_in_immune_domain/List_immune_patient_with_CNV.txt", sep="\t")

#take only the chromosome information needed for genomic range analysis
CNVRangeDF <- ImmuneCNVdata[,c(20:27,6)]
CNVRangeDF$Chromosome <- gsub("chr","",as.character(CNVRangeDF$Chromosome))
CNVRangeonly <- makeGRangesFromDataFrame(CNVRangeDF, keep.extra.columns = TRUE)

#overlap within the genes
ImmuneGeneCNV <- inner_join(CNVRangeDF, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y <= start.x) & (end.y >= end.x))
ImmuneGeneCNV <- unique(ImmuneGeneCNV)
ImmuneGeneCNV$match <- rep("within")
#overlap at the start of the genes
ImmuneGeneCNV2 <- inner_join(CNVRangeDF, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((start.y > start.x) & (start.y < end.x))
ImmuneGeneCNV2 <- unique(ImmuneGeneCNV2)#411

#overlap at the end of the genes
ImmuneGeneCNV3 <- inner_join(CNVRangeDF, PIDgenerange, by="Chromosome", suffix= c(".x",".y")) %>% filter((end.y < end.x) & (end.y > start.x))
ImmuneGeneCNV3 <- unique(ImmuneGeneCNV3)#436

#overlap counted twice in ImmuneGeneCNV2 and ImmuneGeneCNV3
ImmuneGenewithinCNV <- semi_join(ImmuneGeneCNV2,ImmuneGeneCNV3)#403
ImmuneGeneCNVatstart <- anti_join(ImmuneGeneCNV2,ImmuneGeneCNV3) #8
ImmuneGeneCNVatend <- anti_join(ImmuneGeneCNV3,ImmuneGeneCNV2) #33
ImmuneGeneCNVatstart$match <- rep("start")
ImmuneGeneCNVatend$match <- rep("end")



splitdeladdCNV <- split(ImmuneGenewithinCNV, ImmuneGenewithinCNV$dbl)
s0 <- as.data.frame(splitdeladdCNV[[1]])
s1 <- as.data.frame(splitdeladdCNV[[2]])
s3 <- as.data.frame(splitdeladdCNV[[3]])
s4 <- as.data.frame(splitdeladdCNV[[4]])
ImmunegeneCNVdel <- full_join(s0,s1, by= c( "Chromosome", "start.x", "end.x", "condition", "length.gain","length.loss","Plate.key","dbl","dbl.1","start.y","end.y","gene.name"))
ImmunegeneCNVdel$match <- rep("del")
ImmunegeneCNVadd <- full_join(s3,s4, by= c( "Chromosome", "start.x", "end.x", "condition", "length.gain","length.loss","Plate.key","dbl","dbl.1","start.y","end.y","gene.name"))
ImmunegeneCNVadd$match <- rep("add")

####check total match with immunegeneCNV####
ImmuneGeneCNVtotal <- full_join(ImmuneGeneCNV,full_join(ImmunegeneCNVadd, full_join(ImmunegeneCNVdel, full_join(ImmuneGeneCNVatstart, ImmuneGeneCNVatend))))

#rearrange ImmuneGeneCNVtotal collumn
ImmuneGeneCNVtotal <- ImmuneGeneCNVtotal[c(9,1,4,5,6,7,8,2,3,10,11,12,13)]
names(ImmuneGeneCNVtotal)[6:11] <- c("copy.number","read.depth","start","end","gene.start","gene.end")
#save
write.table(ImmuneGeneCNVtotal, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_PID_immune_gene/List_MatchedCNV_immunegene.235patient.txt", row.names = F)

#rearrange ImmuneGeneCNVtotal row
ImmuneGeneCNVtotal <- ImmuneGeneCNVtotal[order(ImmuneGeneCNVtotal$Plate.key, ImmuneGeneCNVtotal$Chromosome, ImmuneGeneCNVtotal$condition),c(1:13)]

#split data frame per Gene to make logic table (if needed)
splitpatient <- split.data.frame(ImmuneGeneCNVtotal, ImmuneGeneCNVtotal$gene.name, drop=TRUE)

p1 <- as.data.frame(splitpatient[[1]])
p10 <- as.data.frame(splitpatient[[10]])
p11 <- as.data.frame(splitpatient[[11]])
p12 <- as.data.frame(splitpatient[[12]])
p13 <- as.data.frame(splitpatient[[13]])
p14 <- as.data.frame(splitpatient[[14]])
p15 <- as.data.frame(splitpatient[[15]])
p16 <- as.data.frame(splitpatient[[16]])
p17 <- as.data.frame(splitpatient[[17]])
p18 <- as.data.frame(splitpatient[[18]])
p19 <- as.data.frame(splitpatient[[19]])
p2 <- as.data.frame(splitpatient[[2]])
p20 <- as.data.frame(splitpatient[[20]])
p21 <- as.data.frame(splitpatient[[21]])
p22 <- as.data.frame(splitpatient[[22]])
p23 <- as.data.frame(splitpatient[[23]])
p24 <- as.data.frame(splitpatient[[24]])
p25 <- as.data.frame(splitpatient[[25]])
p26 <- as.data.frame(splitpatient[[26]])
p27 <- as.data.frame(splitpatient[[27]])
p28 <- as.data.frame(splitpatient[[28]])
p29 <- as.data.frame(splitpatient[[29]])
p3 <- as.data.frame(splitpatient[[3]])
p30 <- as.data.frame(splitpatient[[30]])
p31 <- as.data.frame(splitpatient[[31]])
p32 <- as.data.frame(splitpatient[[32]])
p33 <- as.data.frame(splitpatient[[33]])
p34 <- as.data.frame(splitpatient[[34]])
p35 <- as.data.frame(splitpatient[[35]])
p36 <- as.data.frame(splitpatient[[36]])
p37 <- as.data.frame(splitpatient[[37]])
p38 <- as.data.frame(splitpatient[[38]])
p39 <- as.data.frame(splitpatient[[39]])
p4 <- as.data.frame(splitpatient[[4]])
p40 <- as.data.frame(splitpatient[[40]])
p41 <- as.data.frame(splitpatient[[41]])
p42 <- as.data.frame(splitpatient[[42]])
p43 <- as.data.frame(splitpatient[[43]])
p44 <- as.data.frame(splitpatient[[44]])
p45 <- as.data.frame(splitpatient[[45]])
p5 <- as.data.frame(splitpatient[[5]])
p6 <- as.data.frame(splitpatient[[6]])
p7 <- as.data.frame(splitpatient[[7]])
p8 <- as.data.frame(splitpatient[[8]])
p9 <- as.data.frame(splitpatient[[9]])

CNVcorr <-list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,
                p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,
                p33,p34,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45)%>% reduce(full_join, by = c("Plate.key","Chromosome", "length.gain", "length.loss", "condition", "copy.number","read.depth","start","end","gene.start","gene.end","match"))
#DO NOT USE library(GenomicRanges) etc because it will prevent the function above from working
  
 
CNVcorr <- CNVcorr[,c(1:11,13,12,14:57)]

#get A-Z list of gene that overlaps
ImmuneGeneCNVtotal1 <- ImmuneGeneCNVtotal[order(ImmuneGeneCNVtotal$gene.name), c(1:13)]
Genematched <- unique(ImmuneGeneCNVtotal1$gene.name)
#match it with logic table
names(CNVcorr)[13:57] <- Genematched
CNVcorr <- CNVcorr[order(CNVcorr$Plate.key,CNVcorr$Chromosome,CNVcorr$condition), c(1:57)]

#take important info only from CNV in immune domain patient
Pinfo <- ImmuneCNVdata[,c(1:3,6:12,15,16,19)]
CNVpatientlist <- unique(inner_join(Pinfo, ImmuneGeneCNVtotal, by= "Plate.key"))

write.table(CNVcorr, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_PID_immune_gene/List.immunegenes.CNV.235patients.txt", row.names = F)
write.table(CNVpatientlist, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_PID_immune_gene/List.immunegenes.CNV.235patients+diseaseinfo.txt", row.names = F)

CNVcond <- CNVpatientlist[order(CNVpatientlist$condition), c(1:25)]
splitcondition <- split(CNVcond, CNVcond$condition)

CNVgain <- as.data.frame(splitcondition[[1]])
CNVloss <- as.data.frame(splitcondition[[2]])

CNVgain <- CNVgain[order(CNVgain$Plate.key,CNVgain$Chromosome), c(1:25)]
CNVloss <- CNVloss[order(CNVloss$Plate.key,CNVloss$Chromosome), c(1:25)]

write.table(CNVgain, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_PID_immune_gene/List.immunegenes.CNV.235patients_gainonly.txt", row.names = F)#374 CNV+disease details
write.table(CNVloss, file="/home/phalim/re_gecip/immune/phalim/SVCNV/CNV_in_PID_immune_gene/List.immunegenes.CNV.235patients_lossonly.txt", row.names = F)#1129 CNV+disease details

#At this point, I get 5 tables: 
#1. CNVs of immune domain patient at immune genes                                                       :List_MatchedCNV_immunegene.235patient.txt
#2. CNVs of immune domain patient at immune genes with logic table (listing amount of overlap per gene) :List.immunegenes.CNV235patients_gene-CNV_logictable.xls
#3. CNVs of immune domain patient at immune genes with addition of patient information                  :List.immunegenes.CNV.235patients+diseaseinfo.txt
#4. CNVs of immune domain patient at immune genes with addition of patient information_only gain        :List.immunegenes.CNV.235patients_gainonly.txt
#5. CNVs of immune domain patient at immune genes with addition of patient information_only loss        :List.immunegenes.CNV.235patients_lossonly.txt
-------------------------------------------------------------------

#Trial with Genomic Range- do not work
ImmuneCNVdf <- read.csv("SVCNV/MatchedCNV.txt", sep=" ")
ImmuneGeneDF <- read.csv("SVCNV/Matchedgene.txt", sep=" ")


UniquePID_CNV <- unique(ImmuneCNVlocation)

#check dataframe
head(CNVRangeonly)

#############the ways that donot work --> too many chromosome with same name for multiple location####################
#for finding out the range that overlaps partially
##ranges <- merge(rangesA,rangesB,by="chrom",suffixes=c("A","B"))
##ranges[with(ranges, startB <= startA & stopB >= stopA),]

#for finding out number of overlaps with IRanges
##rangesA <- split(IRanges(rangesA$start, rangesA$stop), rangesA$chrom)
##rangesB <- split(IRanges(rangesB$start, rangesB$stop), rangesB$chrom) -->longer
##-which rangesB wholly contain at least one rangesA?
##ov <- countOverlaps(rangesB, rangesA, type="within")>0

#Range.CNV.Gene <- full_join(CNVRangeDF,GeneRangeDF, by="Chromosome", suffixes=c("CNV","Gene"))
#CNVRange1 <- split(IRanges(CNVRangeonly$start,CNVRangeonly$end), CNVRangeonly$Chromosome)
#GeneRange1 <- split(IRanges(GeneRangeonly$start,GeneRangeonly$end), GeneRangeonly$Chromosome)

#ImmuneCNVlocation <- countOverlaps(GeneRange1, CNVRange1, type="within")>0


##after overlap is found, list down the genes and find out the exon + intron start and end then crossref again with the data. 
##reffer again to GenomicRangesHOWTO page 18
#Genelocationlist <- list of genes that match immune CNV. 
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#library(GenomicFeatures)
#Genelocationlist_txs <- transcriptsBy(txdb, by="gene")[[Genelocationlist]]
#trak2_txs

#
#