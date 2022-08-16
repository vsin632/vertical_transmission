#### New steps classification clean ####
library(tidyverse)

##### Transmission Steps ######
ps_3_clean<-ps3_clean
ps_3_clean
ps_3_rel<-transform_sample_counts(ps_3_clean, function(x) x/sum(x))

ps_Y<-prune_samples(sample_data(ps_3_clean)$Colony=="Yellow", ps_3_clean)
ps_R<-prune_samples(sample_data(ps_3_clean)$Colony=="Red", ps_3_clean)
ps_G<-prune_samples(sample_data(ps_3_clean)$Colony=="Green", ps_3_clean)
ps_B<-prune_samples(sample_data(ps_3_clean)$Colony=="Blue", ps_3_clean)
ps_Father<-prune_samples(sample_data(ps_3_clean)$Colony=="Father", ps_3_clean)

ps_Y_A<-prune_samples(sample_data(ps_Y)$Timepoint=="Alate", ps_Y)
ps_Y_A1<-prune_samples(c(TRUE, FALSE), ps_Y_A)
ps_Y_A2<-prune_samples(c(FALSE, TRUE), ps_Y_A)
ps_R_A<-prune_samples(sample_data(ps_R)$Timepoint=="Alate", ps_R)
ps_R_A1<-prune_samples(c(TRUE, FALSE), ps_R_A)
ps_R_A2<-prune_samples(c(FALSE, TRUE), ps_R_A)
ps_G_A<-prune_samples(sample_data(ps_G)$Timepoint=="Alate", ps_G)
ps_G_A1<-prune_samples(c(TRUE, FALSE), ps_G_A)
ps_G_A2<-prune_samples(c(FALSE, TRUE), ps_G_A)
ps_B_A<-prune_samples(sample_data(ps_B)$Timepoint=="Alate", ps_B)
ps_B_A1<-prune_samples(c(TRUE, FALSE), ps_B_A)
ps_B_A2<-prune_samples(c(FALSE, TRUE), ps_B_A)

ps_F_A<-prune_samples(sample_data(ps_Father)$Timepoint=="Alate", ps_Father)
ps_F_A1<-prune_samples(c(TRUE, FALSE), ps_F_A)
ps_F_A2<-prune_samples(c(FALSE, TRUE), ps_F_A)

ps_Y_D<-prune_samples(sample_data(ps_Y)$Pedigree=="Daughter", ps_Y)
ps_R_D<-prune_samples(sample_data(ps_R)$Pedigree=="Daughter", ps_R)
ps_G_D<-prune_samples(sample_data(ps_G)$Pedigree=="Daughter", ps_G)
ps_B_D<-prune_samples(sample_data(ps_B)$Pedigree=="Daughter", ps_B)

ps_Y_M<-prune_samples(sample_data(ps_Y)$Pedigree=="Mother" &
                        sample_data(ps_Y)$Timepoint!="Alate", ps_Y)
ps_R_M<-prune_samples(sample_data(ps_R)$Pedigree=="Mother"&
                        sample_data(ps_R)$Timepoint!="Alate", ps_R)
ps_G_M<-prune_samples(sample_data(ps_G)$Pedigree=="Mother"&
                        sample_data(ps_G)$Timepoint!="Alate", ps_G)
ps_B_M<-prune_samples(sample_data(ps_B)$Pedigree=="Mother"&
                        sample_data(ps_B)$Timepoint!="Alate", ps_B)

ps_Y_M1<-prune_samples(sample_data(ps_Y_M)$Caste_composite=="Mother colony Worker", ps_Y_M)
ps_R_M1<-prune_samples(sample_data(ps_R_M)$Caste_composite=="Mother colony Worker", ps_R_M)
ps_G_M1<-prune_samples(sample_data(ps_G_M)$Caste_composite=="Mother colony Worker", ps_G_M)
ps_B_M1<-prune_samples(sample_data(ps_B_M)$Caste_composite=="Mother colony Worker", ps_B_M)
ps_Y_M1<-prune_taxa(taxa_sums(ps_Y_M1)>0, ps_Y_M1)
ps_R_M1<-prune_taxa(taxa_sums(ps_R_M1)>0, ps_R_M1)
ps_G_M1<-prune_taxa(taxa_sums(ps_G_M1)>0, ps_G_M1)
ps_B_M1<-prune_taxa(taxa_sums(ps_B_M1)>0, ps_B_M1)


ps_F_M<-prune_samples(sample_data(ps_Father)$Pedigree=="Father"
                      &
                        sample_data(ps_Father)$Timepoint!="Alate", ps_Father)

##
ps_Y_A1<-prune_taxa(taxa_sums(otu_table(ps_Y_A1))>0, ps_Y_A1)
ps_R_A1<-prune_taxa(taxa_sums(otu_table(ps_R_A1))>0, ps_R_A1)
ps_G_A1<-prune_taxa(taxa_sums(otu_table(ps_G_A1))>0, ps_G_A1)
ps_B_A1<-prune_taxa(taxa_sums(otu_table(ps_B_A1))>0, ps_B_A1)
ps_Y_A2<-prune_taxa(taxa_sums(otu_table(ps_Y_A2))>0, ps_Y_A2)
ps_R_A2<-prune_taxa(taxa_sums(otu_table(ps_R_A2))>0, ps_R_A2)
ps_G_A2<-prune_taxa(taxa_sums(otu_table(ps_G_A2))>0, ps_G_A2)
ps_B_A2<-prune_taxa(taxa_sums(otu_table(ps_B_A2))>0, ps_B_A2)

ps_F_A1<-prune_taxa(taxa_sums(otu_table(ps_F_A1))>0, ps_F_A1)
ps_F_A2<-prune_taxa(taxa_sums(otu_table(ps_F_A2))>0, ps_F_A2)

ps_Y_D<-prune_taxa(taxa_sums(otu_table(ps_Y_D))>0, ps_Y_D)
ps_R_D<-prune_taxa(taxa_sums(otu_table(ps_R_D))>0, ps_R_D)
ps_G_D<-prune_taxa(taxa_sums(otu_table(ps_G_D))>0, ps_G_D)
ps_B_D<-prune_taxa(taxa_sums(otu_table(ps_B_D))>0, ps_B_D)
ps_Y_M<-prune_taxa(taxa_sums(otu_table(ps_Y_M))>0, ps_Y_M)
ps_R_M<-prune_taxa(taxa_sums(otu_table(ps_R_M))>0, ps_R_M)
ps_G_M<-prune_taxa(taxa_sums(otu_table(ps_G_M))>0, ps_G_M)
ps_B_M<-prune_taxa(taxa_sums(otu_table(ps_B_M))>0, ps_B_M)
ps_Y_A<-prune_taxa(taxa_sums(otu_table(ps_Y_A))>0, ps_Y_A) 
ps_R_A<-prune_taxa(taxa_sums(otu_table(ps_R_A))>0, ps_R_A)
ps_G_A<-prune_taxa(taxa_sums(otu_table(ps_G_A))>0, ps_G_A)
ps_B_A<-prune_taxa(taxa_sums(otu_table(ps_B_A))>0, ps_B_A)

ps_F_M<-prune_taxa(taxa_sums(otu_table(ps_F_M))>0, ps_F_M)
ps_F_A<-prune_taxa(taxa_sums(otu_table(ps_F_A))>0, ps_F_A)

##
tax_Y_A1<-taxa_names(ps_Y_A1)
tax_R_A1<-taxa_names(ps_R_A1)
tax_G_A1<-taxa_names(ps_G_A1)
tax_B_A1<-taxa_names(ps_B_A1)
tax_Y_A2<-taxa_names(ps_Y_A2)
tax_R_A2<-taxa_names(ps_R_A2)
tax_G_A2<-taxa_names(ps_G_A2)
tax_B_A2<-taxa_names(ps_B_A2)
tax_F_A1<-taxa_names(ps_F_A1)
tax_F_A2<-taxa_names(ps_F_A2)

tax_Y_D<-taxa_names(ps_Y_D)
tax_R_D<-taxa_names(ps_R_D)
tax_G_D<-taxa_names(ps_G_D)
tax_B_D<-taxa_names(ps_B_D)
tax_Y_M<-taxa_names(ps_Y_M)
tax_R_M<-taxa_names(ps_R_M)
tax_G_M<-taxa_names(ps_G_M)
tax_B_M<-taxa_names(ps_B_M)
tax_F_M<-taxa_names(ps_F_M)

tax_Y_A<-taxa_names(ps_Y_A)
tax_R_A<-taxa_names(ps_R_A)
tax_G_A<-taxa_names(ps_G_A)
tax_B_A<-taxa_names(ps_B_A)
tax_F_A<-taxa_names(ps_F_A)

##### Step 1 #### Mother colony (FULL) vs Both alates

# yellow
tax_Y_s1<-subset(tax_Y_M, tax_Y_M %in% tax_Y_A1)
tax_Y_s1<-subset(tax_Y_s1, tax_Y_s1 %in% tax_Y_A2)
tax_R_s1<-subset(tax_R_M, tax_R_M %in% tax_R_A1)
tax_R_s1<-subset(tax_R_s1, tax_R_s1 %in% tax_R_A2)
tax_G_s1<-subset(tax_G_M, tax_G_M %in% tax_G_A1)
tax_G_s1<-subset(tax_G_s1, tax_G_s1 %in% tax_G_A2)
tax_B_s1<-subset(tax_B_M, tax_B_M %in% tax_B_A1)
tax_B_s1<-subset(tax_B_s1, tax_B_s1 %in% tax_B_A2)

tax_F_s1<-subset(tax_F_M, tax_F_M %in% tax_F_A1)
tax_F_s1<-subset(tax_F_s1, tax_F_s1 %in% tax_F_A2)


#### Step 2 ####
ps_Y_D1<-prune_samples(sample_data(ps_Y_D)$Timepoint=="Fasting colony", ps_Y_D)
ps_Y_D1<-prune_taxa(taxa_sums(otu_table(ps_Y_D1))>0, ps_Y_D1)
ps_G_D1<-prune_samples(sample_data(ps_G_D)$Timepoint=="Fasting colony", ps_G_D)
ps_G_D1<-prune_taxa(taxa_sums(otu_table(ps_G_D1))>0, ps_G_D1)
ps_R_D1<-prune_samples(sample_data(ps_R_D)$Timepoint=="Fasting colony", ps_R_D)
ps_R_D1<-prune_taxa(taxa_sums(otu_table(ps_R_D1))>0, ps_R_D1)
ps_B_D1<-prune_samples(sample_data(ps_B_D)$Timepoint=="Fasting colony", ps_B_D)
ps_B_D1<-prune_taxa(taxa_sums(otu_table(ps_B_D1))>0, ps_B_D1)

ps_G_D2<-prune_samples(sample_data(ps_G_D)$Timepoint!="Fasting colony", ps_G_D)
ps_G_D2<-prune_taxa(taxa_sums(otu_table(ps_G_D2))>0, ps_G_D2)
ps_R_D2<-prune_samples(sample_data(ps_R_D)$Timepoint!="Fasting colony", ps_R_D)
ps_R_D2<-prune_taxa(taxa_sums(otu_table(ps_R_D2))>0, ps_R_D2)
ps_B_D2<-prune_samples(sample_data(ps_B_D)$Timepoint!="Fasting colony", ps_B_D)
ps_B_D2<-prune_taxa(taxa_sums(otu_table(ps_B_D2))>0, ps_B_D2)

tax_Y_D1<-taxa_names(ps_Y_D1)
tax_G_D1<-taxa_names(ps_G_D1)
tax_R_D1<-taxa_names(ps_R_D1)
tax_B_D1<-taxa_names(ps_B_D1)

tax_G_D2<-taxa_names(ps_G_D2)
tax_R_D2<-taxa_names(ps_R_D2)
tax_B_D2<-taxa_names(ps_B_D2)

#### From  mother-alate-daughter colony
tax_Y_s2<-subset(tax_Y_s1, tax_Y_s1 %in% tax_Y_D1)
tax_R_s2<-subset(tax_R_s1, tax_R_s1 %in% tax_R_D1)
tax_G_s2<-subset(tax_G_s1, tax_G_s1 %in% tax_G_D1)
tax_B_s2<-subset(tax_B_s1, tax_B_s1 %in% tax_B_D1)

tax_Y_s2F<-subset(tax_F_s1, tax_F_s1 %in% tax_Y_D1)
tax_R_s2F<-subset(tax_F_s1, tax_F_s1 %in% tax_R_D1)
tax_G_s2F<-subset(tax_F_s1, tax_F_s1 %in% tax_G_D1)
tax_B_s2F<-subset(tax_F_s1, tax_F_s1 %in% tax_B_D1)

length(tax_Y_s1)
length(tax_R_s1)
length(tax_G_s1)
length(tax_B_s1)
length(tax_F_s1)
length(unique(c(tax_Y_s1, tax_R_s1, tax_G_s1, tax_B_s1, tax_F_s1)))

#### non-cumulative step 2 ####
tax_Y_s2_nc<-subset(tax_Y_A1, tax_Y_A1 %in% tax_Y_D1)
tax_Y_s2_nc<-subset(tax_Y_s2_nc, tax_Y_s2_nc %in% tax_Y_A2)
tax_R_s2_nc<-subset(tax_R_A1, tax_R_A1 %in% tax_R_D1)
tax_R_s2_nc<-subset(tax_R_s2_nc, tax_R_s2_nc %in% tax_R_A2)
tax_G_s2_nc<-subset(tax_G_A1, tax_G_A1 %in% tax_G_D1)
tax_G_s2_nc<-subset(tax_G_s2_nc, tax_G_s2_nc %in% tax_G_A2)
tax_B_s2_nc<-subset(tax_B_A1, tax_B_A1 %in% tax_B_D1)
tax_B_s2_nc<-subset(tax_B_s2_nc, tax_B_s2_nc %in% tax_B_A2)

length(tax_Y_s2_nc)
length(tax_R_s2_nc)
length(tax_G_s2_nc)
length(tax_B_s2_nc)

tax_Y_s2_nc_F<-subset(tax_F_A1, tax_F_A1 %in% tax_Y_D1)
tax_Y_s2_nc_F<-subset(tax_Y_s2_nc_F, tax_Y_s2_nc_F %in% tax_F_A2)
tax_R_s2_nc_F<-subset(tax_F_A1, tax_F_A1 %in% tax_R_D1)
tax_R_s2_nc_F<-subset(tax_R_s2_nc_F, tax_R_s2_nc_F %in% tax_F_A2)
tax_G_s2_nc_F<-subset(tax_F_A1, tax_F_A1 %in% tax_G_D1)
tax_G_s2_nc_F<-subset(tax_G_s2_nc_F, tax_G_s2_nc_F %in% tax_F_A2)
tax_B_s2_nc_F<-subset(tax_F_A1, tax_F_A1 %in% tax_B_D1)
tax_B_s2_nc_F<-subset(tax_B_s2_nc_F, tax_B_s2_nc_F %in% tax_F_A2)

length(tax_Y_s2_nc_F)
length(tax_R_s2_nc_F)
length(tax_G_s2_nc_F)
length(tax_B_s2_nc_F)


### APPENDIX VERONICAA

#### From  mother-alate-daughter (not alates to queen)
tax_Y_s2b<-subset(tax_Y_s1, tax_Y_s1 %in% tax_Y_D)
tax_R_s2b<-subset(tax_R_s1, tax_R_s1 %in% tax_R_D)
tax_G_s2b<-subset(tax_G_s1, tax_G_s1 %in% tax_G_D)
tax_B_s2b<-subset(tax_B_s1, tax_B_s1 %in% tax_B_D)

tax_Y_s2Fb<-subset(tax_F_s1, tax_F_s1 %in% tax_Y_D)
tax_R_s2Fb<-subset(tax_F_s1, tax_F_s1 %in% tax_R_D)
tax_G_s2Fb<-subset(tax_F_s1, tax_F_s1 %in% tax_G_D)
tax_B_s2Fb<-subset(tax_F_s1, tax_F_s1 %in% tax_B_D)

length(tax_Y_s2b)
length(tax_R_s2b)
length(tax_G_s2b)
length(tax_B_s2b)
length(tax_Y_s2Fb)
length(tax_R_s2Fb)
length(tax_G_s2Fb)
length(tax_B_s2Fb)
length(unique(c(tax_Y_s2b, tax_R_s2b, tax_G_s2b, tax_B_s2b,
                tax_Y_s2Fb, tax_R_s2Fb, tax_G_s2Fb, tax_B_s2Fb
)))


s2_Asv<-unique(c(tax_Y_s2b, tax_R_s2b, tax_G_s2b, tax_B_s2b,
                 tax_Y_s2Fb, tax_R_s2Fb, tax_G_s2Fb, tax_B_s2Fb
))


#### Queen-specific VT taxa
ps_Y_MQ<-prune_samples(sample_data(ps_Y_M)$Caste_composite=="Mother colony Queen", ps_Y_M)
ps_Y_MQ<-subset_taxa(ps_Y_MQ, taxa_sums(ps_Y_MQ)>0)
ps_R_MQ<-prune_samples(sample_data(ps_R_M)$Caste_composite=="Mother colony Queen", ps_R_M)
ps_R_MQ<-subset_taxa(ps_R_MQ, taxa_sums(ps_R_MQ)>0)
ps_G_MQ<-prune_samples(sample_data(ps_G_M)$Caste_composite=="Mother colony Queen", ps_G_M)
ps_G_MQ<-subset_taxa(ps_G_MQ, taxa_sums(ps_G_MQ)>0)
ps_B_MQ<-prune_samples(sample_data(ps_B_M)$Caste_composite=="Mother colony Queen", ps_B_M)
ps_B_MQ<-subset_taxa(ps_B_MQ, taxa_sums(ps_B_MQ)>0)

tax_Y_s1Q<-subset(tax_Y_s1, tax_Y_s1 %in% taxa_names(ps_Y_MQ))
tax_Y_s2Q<-subset(tax_Y_s2b, tax_Y_s2b %in% taxa_names(ps_Y_MQ))
tax_R_s1Q<-subset(tax_R_s1, tax_R_s1 %in% taxa_names(ps_R_MQ))
tax_R_s2Q<-subset(tax_R_s2b, tax_R_s2b %in% taxa_names(ps_R_MQ))
tax_G_s1Q<-subset(tax_G_s1, tax_G_s1 %in% taxa_names(ps_G_MQ))
tax_G_s2Q<-subset(tax_G_s2b, tax_G_s2b %in% taxa_names(ps_G_MQ))
tax_B_s1Q<-subset(tax_B_s1, tax_B_s1 %in% taxa_names(ps_B_MQ))
tax_B_s2Q<-subset(tax_B_s2b, tax_B_s2b %in% taxa_names(ps_B_MQ))

length(tax_Y_s1Q)
length(tax_Y_s2Q)
length(tax_R_s1Q)
length(tax_R_s2Q)
length(tax_G_s1Q)
length(tax_G_s2Q)
length(tax_B_s1Q)
length(tax_B_s2Q)

tax_all_1_Q<-list(tax_Y_s1Q, tax_R_s1Q, tax_G_s1Q, tax_B_s1Q)
tax_all_2_Q<-list(tax_Y_s2Q, tax_R_s2Q, tax_G_s2Q, tax_B_s2Q)

s1q_gen<-list()
s2q_gen<-list()
for(i in 1:length(tax_all_1_Q)){
  s1q_gen[[i]]<-c(unique(tax_table(ps_3_clean)[tax_all_1_Q[[i]],"Genus"]))
  s2q_gen[[i]]<-c(unique(tax_table(ps_3_clean)[tax_all_2_Q[[i]],"Genus"]))
}

Queen_s1_core<-subset(s1q_gen[[1]], s1q_gen[[1]] %in% s1q_gen[[2]] &
                        s1q_gen[[1]] %in% s1q_gen[[3]] &
                        s1q_gen[[1]] %in% s1q_gen[[4]])

Queen_s2_core<-subset(s2q_gen[[1]], s2q_gen[[1]] %in% s2q_gen[[2]] &
                        s2q_gen[[1]] %in% s2q_gen[[3]] &
                        s2q_gen[[1]] %in% s2q_gen[[4]])
Queen_s2_core

#### Linear model on ASV abundances ####
ASVs<-c(length(tax_Y_s1), length(tax_R_s1), length(tax_G_s1), length(tax_B_s1), length(tax_F_s1),
        length(tax_Y_s2b), length(tax_R_s2b), length(tax_G_s2b), length(tax_B_s2b),
        length(tax_Y_s2Fb), length(tax_R_s2Fb), length(tax_G_s2Fb), length(tax_B_s2Fb))

ASVs_nc<-c(length(tax_Y_s1), length(tax_R_s1), length(tax_G_s1), length(tax_B_s1), length(tax_F_s1),
           length(tax_Y_s2_nc), length(tax_R_s2_nc), length(tax_G_s2_nc), length(tax_B_s2_nc),
           length(tax_Y_s2_nc_F), length(tax_R_s2_nc_F), length(tax_G_s2_nc_F), length(tax_B_s2_nc_F))

Steps<-c(1, 1, 1, 1, 1,
         2, 2, 2, 2,
         2, 2, 2, 2)

Lineages<-c("Y", "R", "G", "B", "NA",
            "Y", "R", "G", "B",
            "Y", "R", "G", "B")

Batches<-c("No", "No", "No", "No", "Yes", 
           "Yes", "Yes", "Yes", "Yes",
           "No", "No", "No", "No")

Sexes<-c("M", "M", "M", "M", "F",
         "M", "M", "M", "M",
         "F",  "F",  "F",  "F")

Names<-c("Y_s1", "R_s1", "G_s1", "B_s1", "F_s1",
         "Y_s2", "R_s2", "G_s2", "B_s2", "Y_s2F", "R_s2F", "G_s2F", "B_s2F")


Names2<-c("Y_s1", "R_s1", "G_s1", "B_s1", "F_s1",
          "Y_s2b", "R_s2b", "G_s2b", "B_s2b", "Y_s2Fb", "R_s2Fb", "G_s2Fb", "B_s2Fb")

df<-data.frame(ASVs = ASVs, non_cum_ASVs = ASVs_nc, Step = Steps, Lineage = Lineages, Batch = Batches, Sex = Sexes,
               Name = Names, Name2 = Names2)
head(df)

m1<-lm(data= df, non_cum_ASVs ~ Batch+ Step + Sex)
summary(aov(m1))
aov(m1)
m1
shapiro.test(residuals(m1))
lmtest::bptest(m1)

df001<-subset(df, Step == 1)
m2<-lm(data= df001, non_cum_ASVs ~ Batch + Sex)
summary(aov(m2))
aov(m2)
m2
shapiro.test(residuals(m2))
lmtest::bptest(m2)

df$Step
df %>% group_by(Step) %>%
  summarise(mean = mean(ASVs))


df$Name<-factor(df$Name, levels = Names)

ggplot(df, aes(y=ASVs, x= Name))+geom_bar(stat = "identity")

core_asvs<-length(intersect(intersect(intersect(intersect(tax_Y_M, tax_B_M), tax_G_M), tax_R_M), tax_F_M))

for (i in df$Name2){
  taxa<-get(paste0("tax_", i))
  df[df$Name2 == i, "Genera"]<-length(unique(paste0(tax_table(ps_3_clean)[taxa,"Phylum"], ":_",
                                                    tax_table(ps_3_clean)[taxa,"Class"], ":_",
                                                    tax_table(ps_3_clean)[taxa,"Order"], ":_",
                                                    tax_table(ps_3_clean)[taxa,"Family"], ":_",
                                                    tax_table(ps_3_clean)[taxa,"Genus"])))
}


## Add core genera
tax_Y_M<-unique(paste0(tax_table(ps_Y_M)[,"Phylum"], ":_",
                       tax_table(ps_Y_M)[,"Class"], ":_",
                       tax_table(ps_Y_M)[,"Order"], ":_",
                       tax_table(ps_Y_M)[,"Family"], ":_",
                       tax_table(ps_Y_M)[,"Genus"]))

tax_R_M<-unique(paste0(tax_table(ps_R_M)[,"Phylum"], ":_",
                       tax_table(ps_R_M)[,"Class"], ":_",
                       tax_table(ps_R_M)[,"Order"], ":_",
                       tax_table(ps_R_M)[,"Family"], ":_",
                       tax_table(ps_R_M)[,"Genus"]))

tax_G_M<-unique(paste0(tax_table(ps_G_M)[,"Phylum"], ":_",
                       tax_table(ps_G_M)[,"Class"], ":_",
                       tax_table(ps_G_M)[,"Order"], ":_",
                       tax_table(ps_G_M)[,"Family"], ":_",
                       tax_table(ps_G_M)[,"Genus"]))

tax_B_M<-unique(paste0(tax_table(ps_B_M)[,"Phylum"], ":_",
                       tax_table(ps_B_M)[,"Class"], ":_",
                       tax_table(ps_B_M)[,"Order"], ":_",
                       tax_table(ps_B_M)[,"Family"], ":_",
                       tax_table(ps_B_M)[,"Genus"]))

tax_F_M<-unique(paste0(tax_table(ps_F_M)[,"Phylum"], ":_",
                       tax_table(ps_F_M)[,"Class"], ":_",
                       tax_table(ps_F_M)[,"Order"], ":_",
                       tax_table(ps_F_M)[,"Family"], ":_",
                       tax_table(ps_F_M)[,"Genus"]))

core<-length(intersect(intersect(intersect(intersect(tax_Y_M, tax_B_M), tax_G_M), tax_R_M), tax_F_M))

df$Core_genera<-core
df$Core_ASVs<-core_asvs

ggplot(df)+
  geom_bar(stat = "identity", aes(y=Core_ASVs, x= Name), fill = "orange")+
  geom_bar(stat = "identity", aes(y=ASVs, x= Name), fill = "dark red")+
  geom_bar(stat = "identity", aes(y=Core_genera*10, x= Name), fill = "dark green")+
  geom_bar(stat = "identity", aes(y=Genera*10, x= Name), fill = "blue")+
  ylim(0, NA)+
  ylab("ASVs")


df$Core_genera_10x<-df$Core_genera*10
df$Genera_10x<-df$Genera*10
library(reshape2)
df2<-melt(df, id.vars = c("Name", "Name2", "Sex", "Batch", "Lineage", "Step", "non_cum_ASVs"))

df2$variable<-factor(df2$variable, levels = c("Core_ASVs","ASVs","Core_genera_10x","Genera_10x", 
                                              "Genera", "Core_genera"))

for (i in 1:nrow(df2)){
  if(df2[i, "Name"]=="Y_s1"){
    df2[i, "Name_new"]<-"Yellow_step_1"}
  if(df2[i, "Name"]=="R_s1"){
    df2[i, "Name_new"]<-"Red_step_1"}
  if(df2[i, "Name"]=="G_s1"){
    df2[i, "Name_new"]<-"Green_step_1"}
  if(df2[i, "Name"]=="B_s1"){
    df2[i, "Name_new"]<-"Blue_step_1"}
  if(df2[i, "Name"]=="F_s1"){
    df2[i, "Name_new"]<-"Father_step_1"}
  
  if(df2[i, "Name"]=="Y_s2"){
    df2[i, "Name_new"]<-"Yellow_step_2"}
  if(df2[i, "Name"]=="R_s2"){
    df2[i, "Name_new"]<-"Red_step_2"}
  if(df2[i, "Name"]=="G_s2"){
    df2[i, "Name_new"]<-"Green_step_2"}
  if(df2[i, "Name"]=="B_s2"){
    df2[i, "Name_new"]<-"Blue_step_2"}
  if(df2[i, "Name"]=="Y_s2F"){
    df2[i, "Name_new"]<-"Yellow_step_2_father"}
  if(df2[i, "Name"]=="R_s2F"){
    df2[i, "Name_new"]<-"Red_step_2_father"}
  if(df2[i, "Name"]=="G_s2F"){
    df2[i, "Name_new"]<-"Green_step_2_father"}
  if(df2[i, "Name"]=="B_s2F"){
    df2[i, "Name_new"]<-"Blue_step_2_father"}
}

df3<-subset(df2, !(variable %in% c("Genera", "Core_genera")))
df3$Name_new<-factor(df3$Name_new, levels = c("Yellow_step_1", "Red_step_1", "Blue_step_1", "Green_step_1", "Father_step_1",
                                              "Yellow_step_2", "Red_step_2", "Blue_step_2", "Green_step_2",
                                              "Yellow_step_2_father", "Red_step_2_father", "Blue_step_2_father", "Green_step_2_father"
))

ggplot(df3[order(df3$variable), ],
       aes(y=value, fill = variable, x = Name_new, order = variable))+
  geom_bar(stat = "identity", position = "identity")+
  theme_classic()+
  scale_y_continuous("ASV count",
                     sec.axis = sec_axis(~ . / 10, name = "Genus count"))+
  theme(axis.text.x = element_text(size = 13, angle = 90),
        axis.text.y = element_text(size =13),
        axis.title = element_text(size = 14))




##################
#### Save VTT ####
##################

taxa_Y_s1<-unique(paste0(tax_table(ps_3_clean)[tax_Y_s1, "Phylum"], ":_", tax_table(ps_3_clean)[tax_Y_s1, "Class"], ":_",
                         tax_table(ps_3_clean)[tax_Y_s1, "Order"], ":_", tax_table(ps_3_clean)[tax_Y_s1, "Family"], ":_",
                         tax_table(ps_3_clean)[tax_Y_s1, "Genus"]))

taxa_R_s1<-unique(paste0(tax_table(ps_3_clean)[tax_R_s1, "Phylum"], ":_", tax_table(ps_3_clean)[tax_R_s1, "Class"], ":_",
                         tax_table(ps_3_clean)[tax_R_s1, "Order"], ":_", tax_table(ps_3_clean)[tax_R_s1, "Family"], ":_",
                         tax_table(ps_3_clean)[tax_R_s1, "Genus"]))

taxa_B_s1<-unique(paste0(tax_table(ps_3_clean)[tax_B_s1, "Phylum"], ":_", tax_table(ps_3_clean)[tax_B_s1, "Class"], ":_",
                         tax_table(ps_3_clean)[tax_B_s1, "Order"], ":_", tax_table(ps_3_clean)[tax_B_s1, "Family"], ":_",
                         tax_table(ps_3_clean)[tax_B_s1, "Genus"]))

taxa_G_s1<-unique(paste0(tax_table(ps_3_clean)[tax_G_s1, "Phylum"], ":_", tax_table(ps_3_clean)[tax_G_s1, "Class"], ":_",
                         tax_table(ps_3_clean)[tax_G_s1, "Order"], ":_", tax_table(ps_3_clean)[tax_G_s1, "Family"], ":_",
                         tax_table(ps_3_clean)[tax_G_s1, "Genus"]))

taxa_F_s1<-unique(paste0(tax_table(ps_3_clean)[tax_F_s1, "Phylum"], ":_", tax_table(ps_3_clean)[tax_F_s1, "Class"], ":_",
                         tax_table(ps_3_clean)[tax_F_s1, "Order"], ":_", tax_table(ps_3_clean)[tax_F_s1, "Family"], ":_",
                         tax_table(ps_3_clean)[tax_F_s1, "Genus"]))

##
taxa_Y_s2b<-unique(paste0(tax_table(ps_3_clean)[tax_Y_s2b, "Phylum"], ":_", tax_table(ps_3_clean)[tax_Y_s2b, "Class"], ":_",
                          tax_table(ps_3_clean)[tax_Y_s2b, "Order"], ":_", tax_table(ps_3_clean)[tax_Y_s2b, "Family"], ":_",
                          tax_table(ps_3_clean)[tax_Y_s2b, "Genus"]))
taxa_R_s2b<-unique(paste0(tax_table(ps_3_clean)[tax_R_s2b, "Phylum"], ":_", tax_table(ps_3_clean)[tax_R_s2b, "Class"], ":_",
                          tax_table(ps_3_clean)[tax_R_s2b, "Order"], ":_", tax_table(ps_3_clean)[tax_R_s2b, "Family"], ":_",
                          tax_table(ps_3_clean)[tax_R_s2b, "Genus"]))
taxa_G_s2b<-unique(paste0(tax_table(ps_3_clean)[tax_G_s2b, "Phylum"], ":_", tax_table(ps_3_clean)[tax_G_s2b, "Class"], ":_",
                          tax_table(ps_3_clean)[tax_G_s2b, "Order"], ":_", tax_table(ps_3_clean)[tax_G_s2b, "Family"], ":_",
                          tax_table(ps_3_clean)[tax_G_s2b, "Genus"]))
taxa_B_s2b<-unique(paste0(tax_table(ps_3_clean)[tax_B_s2b, "Phylum"], ":_", tax_table(ps_3_clean)[tax_B_s2b, "Class"], ":_",
                          tax_table(ps_3_clean)[tax_B_s2b, "Order"], ":_", tax_table(ps_3_clean)[tax_B_s2b, "Family"], ":_",
                          tax_table(ps_3_clean)[tax_B_s2b, "Genus"]))
taxa_Y_s2Fb<-unique(paste0(tax_table(ps_3_clean)[tax_Y_s2Fb, "Phylum"], ":_", tax_table(ps_3_clean)[tax_Y_s2Fb, "Class"], ":_",
                           tax_table(ps_3_clean)[tax_Y_s2Fb, "Order"], ":_", tax_table(ps_3_clean)[tax_Y_s2Fb, "Family"], ":_",
                           tax_table(ps_3_clean)[tax_Y_s2Fb, "Genus"]))
taxa_R_s2Fb<-unique(paste0(tax_table(ps_3_clean)[tax_R_s2Fb, "Phylum"], ":_", tax_table(ps_3_clean)[tax_R_s2Fb, "Class"], ":_",
                           tax_table(ps_3_clean)[tax_R_s2Fb, "Order"], ":_", tax_table(ps_3_clean)[tax_R_s2Fb, "Family"], ":_",
                           tax_table(ps_3_clean)[tax_R_s2Fb, "Genus"]))
taxa_G_s2Fb<-unique(paste0(tax_table(ps_3_clean)[tax_G_s2Fb, "Phylum"], ":_", tax_table(ps_3_clean)[tax_G_s2Fb, "Class"], ":_",
                           tax_table(ps_3_clean)[tax_G_s2Fb, "Order"], ":_", tax_table(ps_3_clean)[tax_G_s2Fb, "Family"], ":_",
                           tax_table(ps_3_clean)[tax_G_s2Fb, "Genus"]))
taxa_B_s2Fb<-unique(paste0(tax_table(ps_3_clean)[tax_B_s2Fb, "Phylum"], ":_", tax_table(ps_3_clean)[tax_B_s2Fb, "Class"], ":_",
                           tax_table(ps_3_clean)[tax_B_s2Fb, "Order"], ":_", tax_table(ps_3_clean)[tax_B_s2Fb, "Family"], ":_",
                           tax_table(ps_3_clean)[tax_B_s2Fb, "Genus"]))

genera_vtt_list<-list(taxa_Y_s1, taxa_R_s1, taxa_B_s1, taxa_G_s1, taxa_F_s1, 
                      taxa_Y_s2b, taxa_R_s2b, taxa_B_s2b, taxa_G_s2b,
                      taxa_Y_s2Fb, taxa_R_s2Fb, taxa_B_s2Fb, taxa_G_s2Fb
)
names(genera_vtt_list)<-c("taxa_Y_s1", "taxa_R_s1", "taxa_B_s1", "taxa_G_s1", "taxa_F_s1", 
                          "taxa_Y_s2b", "taxa_R_s2b", "taxa_B_s2b", "taxa_G_s2b",
                          "taxa_Y_s2Fb", "taxa_R_s2Fb", "taxa_B_s2Fb", "taxa_G_s2Fb"
)

summary(genera_vtt_list)
unique(unlist(genera_vtt_list))
library(rlist)
# list.save(genera_vtt_list, file = "~/Desktop/Bin_stats/VTT_genera_list.yaml")

#### Make sup table with VT ####

taxa<-unique(c(taxa_names(ps_B_M), taxa_names(ps_R_M), taxa_names(ps_G_M), taxa_names(ps_Y_M)))
data_fr<-data.frame(matrix(nrow=length(taxa), ncol = 19))
colnames(data_fr)<-c("ASV", "Phylum", "Class", "Order", "Family", "Genus",
                     "Yellow_step_1", "Red_step_1", "Blue_step_1", "Green_step_1", "Father_step_1",
                     "Yellow_step_2", "Red_step_2", "Blue_step_2", "Green_step_2",
                     "Yellow_step_2_father", "Red_step_2_father", "Blue_step_2_father", "Green_step_2_father")

data_fr$ASV<-taxa
data_fr$Phylum<-tax_table(ps_3_clean)[data_fr$ASV, "Phylum"]
data_fr$Class<-tax_table(ps_3_clean)[data_fr$ASV, "Class"]
data_fr$Order<-tax_table(ps_3_clean)[data_fr$ASV, "Order"]
data_fr$Family<-tax_table(ps_3_clean)[data_fr$ASV, "Family"]
data_fr$Genus<-tax_table(ps_3_clean)[data_fr$ASV, "Genus"]

for (i in 1:nrow(data_fr)){
  if(data_fr[i, "ASV"] %in% tax_Y_s1){
    data_fr[i, "Yellow_step_1"]<-"Yes"
  }else(data_fr[i, "Yellow_step_1"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_R_s1){
    data_fr[i, "Red_step_1"]<-"Yes"
  }else(data_fr[i, "Red_step_1"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_G_s1){
    data_fr[i, "Green_step_1"]<-"Yes"
  }else(data_fr[i, "Green_step_1"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_B_s1){
    data_fr[i, "Blue_step_1"]<-"Yes"
  }else(data_fr[i, "Blue_step_1"]<-"No")
  
  if(data_fr[i, "ASV"] %in% tax_Y_s2b){
    data_fr[i, "Yellow_step_2"]<-"Yes"
  }else(data_fr[i, "Yellow_step_2"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_R_s2b){
    data_fr[i, "Red_step_2"]<-"Yes"
  }else(data_fr[i, "Red_step_2"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_G_s2b){
    data_fr[i, "Green_step_2"]<-"Yes"
  }else(data_fr[i, "Green_step_2"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_B_s2b){
    data_fr[i, "Blue_step_2"]<-"Yes"
  }else(data_fr[i, "Blue_step_2"]<-"No")
  
  if(data_fr[i, "ASV"] %in% tax_Y_s2Fb){
    data_fr[i, "Yellow_step_2_father"]<-"Yes"
  }else(data_fr[i, "Yellow_step_2_father"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_R_s2Fb){
    data_fr[i, "Red_step_2_father"]<-"Yes"
  }else(data_fr[i, "Red_step_2_father"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_G_s2Fb){
    data_fr[i, "Green_step_2_father"]<-"Yes"
  }else(data_fr[i, "Green_step_2_father"]<-"No")
  if(data_fr[i, "ASV"] %in% tax_B_s2Fb){
    data_fr[i, "Blue_step_2_father"]<-"Yes"
  }else(data_fr[i, "Blue_step_2_father"]<-"No")
  
  if (i %% 1000 == 0){
    print(i)
  }
}

# write.csv2(data_fr, "~/Desktop/Bin_stats/Sup_data_paper/table_s5.csv")



summary(genera_vtt_list)
step_2 <- genera_vtt_list[6:13]
summary(step_2)
step_2_gen<-subset(step_2[["taxa_Y_s2b"]], step_2[["taxa_Y_s2b"]] %in% step_2[["taxa_R_s2b"]] | step_2[["taxa_Y_s2b"]] %in% step_2[["taxa_R_s2Fb"]])
step_2_gen<-subset(step_2_gen, step_2_gen %in% step_2[["taxa_G_s2b"]] | step_2_gen %in% step_2[["taxa_G_s2Fb"]])
step_2_gen<-subset(step_2_gen, step_2_gen %in% step_2[["taxa_B_s2b"]] | step_2_gen %in% step_2[["taxa_B_s2Fb"]])
length(step_2_gen)

step_2_gen_father<-subset(step_2[["taxa_Y_s2Fb"]], step_2[["taxa_Y_s2Fb"]] %in% step_2[["taxa_R_s2b"]] | step_2[["taxa_Y_s2Fb"]] %in% step_2[["taxa_R_s2Fb"]])
step_2_gen_father<-subset(step_2_gen_father, step_2_gen_father %in% step_2[["taxa_G_s2b"]] | step_2_gen_father %in% step_2[["taxa_G_s2Fb"]])
step_2_gen_father<-subset(step_2_gen_father, step_2_gen_father %in% step_2[["taxa_B_s2b"]] | step_2_gen_father %in% step_2[["taxa_B_s2Fb"]])
length(step_2_gen_father)

length(unique(c(step_2_gen, step_2_gen_father)))/
  length(unique(unlist(genera_vtt_list)))

