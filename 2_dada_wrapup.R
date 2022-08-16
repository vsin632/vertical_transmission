### Both batches run together ####
#########
library(dada2)
library(phyloseq)
library(tidyverse)

sqtab<-readRDS("~/Desktop/Bin_stats/Both_together/seqtab.nochim.rds")
taxa_dictdb_cp<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_cp_2.rds")
dim(taxa_dictdb_cp)
dim(sqtab)
tx_tb<-tax_table(taxa_dictdb_cp)[-45017,]
otu_tb<-otu_table(sqtab, taxa_are_rows = F)

sdf<-read.table("~/Desktop/Bin_stats/Batch_2/batch_2_metadata.csv", sep=";", header = T)
rownames(sdf)<-sdf$SampleID

sam_df_1<-read.table("~/Desktop/Bin_stats/Batch_1/Batch_1_new_metadata.csv", header = T, sep = ";")
rownames(sam_df_1)<-sam_df_1$SampleID
sam_df_1$Batch<-1
colnames(sdf)<-c(colnames(sdf)[1:8], "Batch", "Caste_composite")

sdf<-rbind(sdf, sam_df_1)

sample_names(otu_tb)

ps3_mocks<-phyloseq(tax_table=tx_tb,
                    otu_table=otu_tb,
                    sample_data = sample_data(sdf))

ps3<-ps3_mocks
ps3c<-subset_taxa(ps3, taxa_sums(ps3)>25)
ps3c

# sample_data(ps3)$Batch_BGI<-c(rep("1", 45), rep("2", 22)) # Tree made here

ps3<-prune_samples(sample_data(ps3)$Pedigree!="Mock" &
                     !is.na(sample_data(ps3)$Pedigree), ps3)
sample_data(ps3)$Timepoint<-factor(sample_data(ps3)$Timepoint,
                                   levels = c("Mature colony", "Alate", "Fasting colony", "Fed colony"))
sdf<-sample_data(ps3)
sdf$Pedigree<-factor(sdf$Pedigree, levels = c("Mother", "Father", "Daughter"))
sdf$Caste<-factor(sdf$Caste)
sdf$Caste_composite<-as.character(sdf$Caste_composite)
sdf[sdf$Caste_composite=="Mature colony Queen", "Caste_composite"]<-"Mother colony Queen"
sdf[sdf$Caste_composite=="Mature colony Worker", "Caste_composite"]<-"Mother colony Worker"
sdf[sdf$Caste_composite=="Mature colony King", "Caste_composite"]<-"Mother colony King"
sdf[sdf$Caste_composite=="Alate Queen", "Caste_composite"]<-"Female Alate"
sdf[sdf$Caste_composite=="Alate King", "Caste_composite"]<-"Male Alate"
sdf[sdf$Caste_composite=="Fasting colony Larva", "Caste_composite"]<-"T1 colony Larva"
sdf[sdf$Caste_composite=="Fasting colony Worker", "Caste_composite"]<-"T1 colony Worker"
sdf[sdf$Caste_composite=="Fasting colony Queen", "Caste_composite"]<-"T1 colony Queen"
sdf[sdf$Caste_composite=="Fed colony Queen", "Caste_composite"]<-"T2 colony Queen"
sdf$Caste_2<-as.character(sdf$Caste_2)
sdf[sdf$Caste_composite=="Mother colony Queen" & is.na(sdf$Caste_2), "Caste_2"]<-"Section_C"

sdf$Caste_composite<-factor(sdf$Caste_composite, levels = c("Mother colony Queen", "Mother colony King", "Mother colony Worker",
                                                            "Female Alate", "Father colony Worker", "Male Alate", "T1 colony Larva",
                                                            "T1 colony Worker", "T1 colony Queen", "T2 colony Queen"))
sdf$SampleID<-factor(sdf$SampleID, levels = c("1", "13", "19", "7", "Blue-A", "Blue-B", "Green-A", "Green-B", "Yellow-A", "Yellow-B", "Red-A", "Red-B",
                                              "14", "2", "8", "16c", "4b", "110", "111", "112", "113", "114","115", "116", "117",
                                              "11b", "12c", "17", "18c", "5","51", "57c", "6c", "103", "104", "105", "101", "102",
                                              "25", "35", "42", "54", "23", "24", "26", "27", "33", "34b", "36", "37b", "43", "44",
                                              "46b", "47b", "52b", "53b", "55", "56b", "20c", "30b", "40c", "50c", "21c", "31c", "41c"))

sample_data(ps3)<-sdf
ps3_clean<-subset_samples(ps3,
                          sample_data(ps3)[,"SampleID"] !="4b")
ps3_clean<-subset_samples(ps3_clean,
                          sample_data(ps3_clean)[,"SampleID"] !="16c")
ps3_clean<-prune_taxa(taxa_sums(ps3_clean)>25, ps3_clean)
ps3_clean

### How many classified taxa? ####
tx_tb<-as.data.frame(tax_table(ps3_clean)@.Data)
class(tx_tb)
otu_tb<-as.data.frame(otu_table(ps3_clean)@.Data)

1-nrow(tx_tb[is.na(tx_tb$Genus),])/nrow(tx_tb)
1-nrow(tx_tb[is.na(tx_tb$Family),])/nrow(tx_tb)
1-nrow(tx_tb[is.na(tx_tb$Phylum),])/nrow(tx_tb)

g_tax<-rownames(tx_tb[is.na(tx_tb$Genus),])
1-sum(otu_tb[,g_tax])/sum(otu_tb)

f_tax<-rownames(tx_tb[is.na(tx_tb$Family),])
1-sum(otu_tb[,f_tax])/sum(otu_tb)

p_tax<-rownames(tx_tb[is.na(tx_tb$Phylum),])
1-sum(otu_tb[,p_tax])/sum(otu_tb)


length(unique(tax_table(ps3_clean)[,"Phylum"]))

