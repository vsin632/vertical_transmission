####################################################
######## Batch 2 16S vertical transmission #########
####################################################

library(dada2)
library(phyloseq)
library(tidyverse)
library(phangorn)
library(DECIPHER)
library(Biostrings)
library(ggplot2)
library(reshape2)

#########################
#### Import Raw Data ####
#########################

path <- "~/Path/to/file" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
## OWN DATA: The string manipulations may have to be modified if your filename format is different.


###########################
#### Plot Data Quality ####
###########################

plotQualityProfile(fnFs, TRUE)
plotQualityProfile(fnRs, TRUE)

## OWN DATA: Your reads must still overlap after truncation in order to merge them later!
## The tutorial is using 2x250 V4 sequence data, so the forward and reverse reads almost completely overlap:
## Our trimming can be completely guided by the quality scores.
## If you are using a less-overlapping primer set, like V1-V2 or V3-V4, your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Important step: truncLen is the amount of bases you leave wifthout cutting from F and R reads
# maxN: Maximum number of N bases, dada doesn't allow Ns!
# maxEE: reads with higher expected errors than X will be discarded. I usually relaz to 5 or 10 for reverse reads.
# trunQ: truncate reads at first instance of quality less than X
# type "?filterAndTrim" for more information.


################################################################
#### Remove Adapters and PhiX, and Filter Low Quality Reads ####
################################################################

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,260),   # May take a minute or two.
                     maxN=0, maxEE=c(1, 1), truncQ=1, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)                   # On Windows set multithread=FALSE

# saveRDS(out, "~/Desktop/Bin_stats/Batch_2/out.rds")
# out<-readRDS("~/Desktop/Bin_stats/Batch_2/out.rds")

head(out)

## OWN DATA: The standard filtering parameters are starting points, not set in stone.
## If you want to speed up downstream computation, consider tightening maxEE.
## If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails.
## Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.

## ITS DATA: NO truncLen as sizes vary so very much.

############################
#### Learn Error Rates #####
############################

set.seed(123)
errF <- learnErrors(filtFs, multithread=TRUE, verbose = T, nbases = 2e8, randomize = T)
plotErrors(errF, nominalQ = T)

set.seed(123)
errF <- learnErrors(filtRs, multithread=TRUE, verbose = T, nbases = 3e8, MAX_CONSIST = 100)
plotErrors(errR, nominalQ=TRUE)

# The black lines (error estimates from your data) should match the red lines (theoretical error values)


# for small, poor quality datasets use pool=true to increase accuracy of error correction at expense of computing time. ALWAYS ON SAMPLES FROM A UNIQUE BATCH!!!

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = T) 
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = T)
dadaFs[[1]]


##########################################
#### Merge Forward and Reverse Reads #####
##########################################

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 20) # Increase minOverlap above 12 (default) to reduce chimeras downstream
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
## OWN DATA: most reads should succesfuly merge, did you truncate away the overlap region?


seqtab <- makeSequenceTable(mergers)
dim(seqtab) # dim 1: samples, dim 2: taxa

###################################################
#### Inspect distribution of sequence lengths #####
###################################################

table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))

# Are any sequences far too long or short? Feel free to trim them away with the line below.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] # Change 250-256 to desired lengths. BEWARE Archae have atypical sizes! Do you want to loose them?


#########################################
#### Remove Chimera (Bimera) DeNovo #####
#########################################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

saveRDS(seqtab.nochim, "~/Path/to/file/seqtabnochim.rds")
saveRDS(seqtab, "~/Path/to/file/seqtab.rds")

sum(seqtab.nochim)/sum(seqtab) # Proportion of non-chimera reads. Anything above 0.8 is very good.
## If you get less, consider modifying upstream parameters in filterAndTrim to be more stringent.
## It's not uncommon to loose up to 40% in this step. But make sure you removed primers as these cause a lot of chimera (they have ambiguous bases)


############################################################
#### Track Reads Through Pipeline: How Well Did it Go? #####
############################################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) # Anything above 75% of reads kept in the nonchim column is grat.
colSums(track) # It's better to inspect sample by sample (line above), but looking at the total is also good.



##########################
#### Assign taxonomy #####
##########################


### Both batches run together ####
library(dada2)
library(phyloseq)
library(tidyverse)

sqtab<-readRDS("~/Desktop/Bin_stats/Both_together/seqtab.nochim.rds")

taxa_dictdb_2_test <- assignTaxonomy(sqtab[,1:5000],  "~/Desktop/Dict_db.fasta",
                                     multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test.rds")

taxa_dictdb_2_test_2 <- assignTaxonomy(sqtab[,5001:10000],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_2, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_2.rds")

taxa_dictdb_2_test_3 <- assignTaxonomy(sqtab[,10001:15000],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_3, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_3.rds")

taxa_dictdb_2_test_4 <- assignTaxonomy(sqtab[,15001:20000],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_4, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_4.rds")

taxa_dictdb_2_test_5 <- assignTaxonomy(sqtab[,20001:30000],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_5, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_5.rds")

taxa_dictdb_2_test_6 <- assignTaxonomy(sqtab[,30001:40000],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_6, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_6.rds")

taxa_dictdb_2_test_7 <- assignTaxonomy(sqtab[,40001:45016],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_7, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_7.rds")

taxa_dictdb_2_test_8 <- assignTaxonomy(sqtab[,45016:50016],  "~/Desktop/Dict_db.fasta",
                                       multithread=TRUE, taxLevels=c("Phylum", "Class","Order", "Family", "Genus", "Species"))
saveRDS(taxa_dictdb_2_test_8, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_8.rds")



taxa_dictdb_2_test<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test.rds")
taxa_dictdb_2_test_2<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_2.rds")
taxa_dictdb_2_test_3<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_3.rds")
taxa_dictdb_2_test_4<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_4.rds")
taxa_dictdb_2_test_5<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_5.rds")
taxa_dictdb_2_test_6<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_6.rds")
taxa_dictdb_2_test_7<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_7.rds")
taxa_dictdb_2_test_8<-readRDS("~/Desktop/Bin_stats/Both_together/taxa_dictdb_2_test_8.rds")

taxa_dict_b2<-rbind(taxa_dictdb_2_test,
                    taxa_dictdb_2_test_2,
                    taxa_dictdb_2_test_3,
                    taxa_dictdb_2_test_4,
                    taxa_dictdb_2_test_5,
                    taxa_dictdb_2_test_6,
                    taxa_dictdb_2_test_7,
                    taxa_dictdb_2_test_8[2:nrow(taxa_dictdb_2_test_8),]
)

reclass<-subset(taxa_dict_b2, is.na(taxa_dict_b2[,"Genus"]))
dim(reclass)
head(reclass)


taxa_silva_2 <- assignTaxonomy(reclass,  "~/Downloads/MiSeq_SOP/silva_nr_v132_train_set.fa.gz",
                               multithread=TRUE)
saveRDS(taxa_silva_2, "~/Desktop/Bin_stats/Both_together/taxa_silva_2.rds")

taxa_silva_2_cp<-taxa_silva_2[,2:6]
taxa_dictdb_cp<-taxa_dict_b2[,1:6]

for (r in rownames(taxa_silva_2_cp)){
  for (c in colnames(taxa_silva_2_cp)){
    taxa_dictdb_cp[r, c]<-taxa_silva_2_cp[r,c]
  }
}
saveRDS(taxa_dictdb_cp, "~/Desktop/Bin_stats/Both_together/taxa_dictdb_cp_2.rds")


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

# saveRDS(ps3, "~/Desktop/Bin_stats/ps3.rds")
# saveRDS(sdf, "~/Desktop/Bin_stats/sdf.rds")
# saveRDS(ps3_clean, "~/Desktop/Bin_stats/ps_3_clean_both.rds")

# write.csv2(sdf, "~/Desktop/ps_metadata.csv")

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









