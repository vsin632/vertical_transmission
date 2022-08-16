#### Mock analysis anew ####

ps3_mocks

ps_3_clean.prop <- transform_sample_counts(ps_3_clean, function(x) x/sum(x))

sample_data(ps3_mocks)
ps3_mocks<-prune_taxa(taxa_sums(ps3_mocks)>25, ps3_mocks)
ps_mock<-prune_samples(c("121", "mock"), ps3_mocks)
ps_mock<-prune_taxa(taxa_sums(ps_mock)>0, ps_mock)
ps_mock

tax_table(ps_mock)[,"Genus"]<-paste0(tax_table(ps_mock)[,"Family"], ":_", tax_table(ps_mock)[,"Genus"])
plot_bar(ps_mock, fill = "Genus")

# taxa_names(ps_mock)<-paste0("ASV_", 1:length(taxa_names(ps_mock)), ":_", tax_table(ps_mock)[,"Genus"])

ps_mock_gen<-tax_glom(ps_mock, "Genus")
plot_bar(ps_mock_gen, fill = "Genus")
ps_mock_gen_prop<-transform_sample_counts(ps_mock_gen, function(x) x/sum(x))
plot_bar(ps_mock_gen_prop, fill = "Genus")

df<-data.frame(Genus = tax_table(ps_mock)[,"Genus"], Ab = taxa_sums(ps_mock))
df

sum(rownames(df) %in% rownames(df))
sum(rownames(df) %in% tax_Y_s2F)


ps3<-subset_taxa(ps_3_clean.prop, taxa_names(ps_3_clean.prop) %in% taxa_names(ps_mock))
plot_bar(ps3, fill = "Genus", x = "SampleID")+facet_grid(Caste_composite~Colony, scale = "free_x", space = "free")

ps3g<-tax_glom(ps3, "Genus")
plot_bar(ps3g, fill = "Genus", x = "SampleID")+facet_grid(Caste_composite~Colony, scale = "free_x", space = "free")


### Study intermock variation. ###

tax_table(ps_mock)[tax_table(ps_mock)[, "Genus"] == "Enterococcaceae:_Enterococcus_3", "Genus"]<-"Enterococcaceae:_Enterococcus"
tax_table(ps_mock)[tax_table(ps_mock)[, "Genus"] == "Lactobacillaceae:_Lactobacillus_2", "Genus"]<-"Lactobacillaceae:_Lactobacillus"
tax_table(ps_mock)[tax_table(ps_mock)[, "Genus"] == "Lactobacillaceae:_Lactobacillus_2", "Genus"]<-"Lactobacillaceae:_Lactobacillus"


ps_mock
ps_mock_rel<-transform_sample_counts(ps_mock, function(x) x/sum(x))
ps_mock_rel_g<-tax_glom(ps_mock_rel, "Genus")

a<-c()

ntaxa(ps_mock_rel_g)
for (i in taxa_names(ps_mock_rel_g)){
  focal_ps<-subset_taxa(ps_mock_rel_g, taxa_names(ps_mock_rel_g)==i)
  # print(sample_sums(focal_ps))
  print(((sample_sums(focal_ps)[1]-sample_sums(focal_ps)[2])/sample_sums(focal_ps)[2])*100)
  a<-c(a, ((sample_sums(focal_ps)[1]-sample_sums(focal_ps)[2])/sample_sums(focal_ps)[2])*100)
  
}
a
a<-subset(a, is.finite(a) & a>-100)
a
mean(abs(a))
