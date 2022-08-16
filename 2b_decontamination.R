#### ddPCR data ####
library(tidyverse)
library(ggrepel)
library(phyloseq)
ddp<-read.csv2("~/Desktop/Bin_stats/ddPCR/ddpcr.csv")
ddp<-ddp[,1:13]
summary(ddp)

ddp$Corrected_16S<-ddp$Conc_16S*ddp$Dilution_ddpcr*ddp$Dilution_Justinn # Multiply by dilutions
ddp$Corrected_16S_2<-(ddp$Conc_16S*ddp$Dilution_ddpcr*ddp$Dilution_Justinn)/ddp$Nr_termites # And divide by number of termites in sample
ddp$Corrected_RPL_18<-ddp$Conc_RPL_18*ddp$Dilution_ddpcr*ddp$Dilution_Justinn # Multiply by dilutions
ddp$Prop_16S_RPL18<-ddp$Corrected_16S/ddp$Corrected_RPL_18 # Make the proportion

ggplot(ddp, aes(y=Conc_16S, x= Caste_composite, col = factor(Dilution_ddpcr)))+
  geom_point()+
  theme(axis.text.x = element_text (size = 15, angle = 90))

ggplot(ddp, aes(y=Corrected_16S, x= Caste_composite, col = factor(Dilution_ddpcr), label = SampleID))+
  geom_point()+
  theme(axis.text.x = element_text (size = 15, angle = 90))+
  scale_y_log10()+
  geom_text()+
  geom_hline(yintercept = 1e5, linetype = "dashed")+
  facet_grid(Batch~1)

ggplot(ddp, aes(y=Corrected_RPL_18, x= Caste_composite, col = factor(Dilution_ddpcr), label = SampleID))+
  geom_point()+
  theme(axis.text.x = element_text (size = 15, angle = 90))+
  scale_y_log10()+
  geom_text_repel()+
  geom_hline(yintercept = 1e5, linetype = "dashed")+
  facet_grid(Batch~1)

ggplot(ddp, aes(y=Prop_16S_RPL18, x= Caste_composite, col = Batch, label = SampleID))+
  geom_point()+
  theme(axis.text.x = element_text (size = 15, angle = 90))+
  scale_y_log10()+
  geom_text_repel()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_grid(Batch~1)

ddp2<-subset(ddp, !(Caste_composite %in% c("T2 colony Queen", "", "Control")))
ddp2$Caste_composite<-factor(ddp2$Caste_composite, levels = c(
  "Mother colony Queen", "Mother colony King", "Mother colony Worker", "Father colony Worker", "Female Alate",
  "Male Alate", "T1 colony Queen", "T1 colony Larva", "T1 colony Worker"
))
ggplot(ddp2, aes(y=Corrected_16S_2, x= Caste_composite, fill = Caste_composite))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text (size = 13, angle = 90),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14))+
  scale_y_continuous(limits = c(1, NA), trans = "log10", breaks = c(1, 10, 100, 1000, 10000, 100000))+
  ylab("16S rRNA gene copies / termite")+
  xlab(NULL)+
  guides(fill = F) # 7 - 6

##### DECONTAM #####
library("decontam")
sdf<-readRDS("~/Desktop/Bin_stats/sdf.rds") # Read metadata
sdf2<-as.data.frame(sdf@.Data) # Change format from phyloseq to data frame
colnames(sdf2)<-colnames(sdf) # Make sure column names are ok

sdf2$ddPCR_concentration<-NA # Make new column
for ( i in unique(sdf$SampleID)){
  sdf2[sdf2$SampleID == i, "ddPCR_concentration"]<-c(
    ddp[ddp$SampleID==i, "Prop_16S_RPL18"][1])
}

sdf2 # Check things are ok
sdf2$ddPCR_concentration


ps<-readRDS("~/Desktop/Bin_stats/ps3.rds") # read phyloseq
rownames(sdf2)<-sdf2$SampleID # Make rownames the sample names
sample_data(ps)<-sdf2 # Put sample data into phyloseq object

## Frequency driven
set.seed(123)
contamdf.freq <- isContaminant(ps, method="frequency", conc="ddPCR_concentration", normalize = F) # Calculate contaminants
head(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))
plot_frequency(ps, taxa_names(ps)[c(1,3)], conc="ddPCR_concentration") + 
  xlab("DNA Concentration (ddPCR signal)")
set.seed(100)
plot_frequency(ps, taxa_names(ps)[which(contamdf.freq$contaminant)], conc="ddPCR_concentration") +
  xlab("DNA Concentration (ddPCR signal)")+
  facet_grid(~taxa, labeller = labeller(taxa = tax_table(ps)[which(contamdf.freq$contaminant), "Genus"]))

contaminants_freq<-subset(contamdf.freq, contaminant == T) 

paste0(unique(tax_table(ps)[rownames(contaminants_freq), c( "Genus")])) # Which are they?


head(contaminants_freq)


contaminants_freq[, "Genus"]

ps_q<-subset_samples(ps, sample_data(ps)$Caste_composite== "Mother colony Queen")
ps_q_rel<-transform_sample_counts(ps_q, function(x) x/sum(x))

ps_k<-subset_samples(ps, sample_data(ps)$Caste_composite== "Mother colony King")
ps_k_rel<-transform_sample_counts(ps_k, function(x) x/sum(x))

ps_w<-subset_samples(ps, sample_data(ps)$Caste_composite== "Mother colony Worker")
ps_w_rel<-transform_sample_counts(ps_w, function(x) x/sum(x))


rownames(contaminants_freq)-> contaminants_freq$ASV
for (i in rownames(contaminants_freq)){
  contaminants_freq[i, c("Phylum", "Class", "Order", "Family", "Genus")]<-c(tax_table(ps)[i, c("Phylum", "Class", "Order", "Family", "Genus")])
  focal_ps_q<-subset_taxa(ps_q_rel, taxa_names(ps_q_rel) == i)
  contaminants_freq[i, "Ab_in_queens"]<-mean(sample_sums(focal_ps_q))
  contaminants_freq[i, "sd_in_queens"]<-sd(sample_sums(focal_ps_q))
  
  focal_ps_k<-subset_taxa(ps_k_rel, taxa_names(ps_k_rel) == i)
  contaminants_freq[i, "Ab_in_kings"]<-mean(sample_sums(focal_ps_k))
  contaminants_freq[i, "sd_in_kings"]<-sd(sample_sums(focal_ps_k))
  
  focal_ps_w<-subset_taxa(ps_w_rel, taxa_names(ps_w_rel) == i)
  contaminants_freq[i, "Ab_in_workers"]<-mean(sample_sums(focal_ps_w))
  contaminants_freq[i, "sd_in_workers"]<-sd(sample_sums(focal_ps_w))
}

rownames(contaminants_freq)<-NULL

head(contaminants_freq, n = 100)
head(subset(contaminants_freq, Ab_in_queens > 0), n = 100)


rownames(contaminants_freq)<-NULL
contaminants_freq
write.csv2(contaminants_freq, "~/Desktop/Bin_stats/Sup_data_paper/Table_s9.csv")


######## How much of the microbiome do they make up? ###
ps_rel<-transform_sample_counts(ps, function(x) x/sum(x))
ps_contam_all<-subset_taxa(ps_rel, taxa_names(ps_rel) %in% contaminants_freq$ASV)
plot_bar(ps_contam_all, x = "SampleID")+
  theme(axis.text.y = element_text(size = 12))

mean(sample_sums(ps_contam_all))


