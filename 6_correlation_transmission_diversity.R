#### Workers abundance predicts vertical transmission ####
ps3_clean

Lineages<-c("Yellow", "Red", "Green", "Blue")


all_tax<-list(tax_Y_s1, tax_B_s1, tax_R_s1, tax_G_s1, tax_F_s1,
              tax_Y_s2b, tax_B_s2b, tax_R_s2b, tax_G_s2b,
              tax_Y_s2Fb, tax_B_s2Fb, tax_R_s2Fb, tax_G_s2Fb)

cols<-c("Yellow", "Blue", "Red", "Green", "Father",
        "Yellow", "Blue", "Red", "Green",
        "Yellow", "Blue", "Red", "Green")

Sex<-c("M", "M", "M", "M", "F",
       "M", "M", "M", "M",
       "F", "F", "F", "F")

Step<-c("1", "1", "1", "1", "1",
        "2", "2", "2", "2",
        "2", "2", "2", "2")

all_taxa<-list()
all_sum<-list()
for (i in 1:length(all_tax)){ # Get the strain taxonomic classification out and count it.
  all_taxa[[i]]<-tax_table(ps_3_clean)[all_tax[[i]]]
  all_sum[[i]]<-all_taxa[[i]][,c("Phylum", "Class", "Order", "Family", "Genus")] %>%
    plyr::count()
  print((i/0.13))
}

for (xy in 1:length(all_sum)){ # Paste all tax ranks into single vector called "Taxon" to compare and count properly
  for (i in 1:nrow(all_sum[[xy]])){
    all_sum[[xy]][i, "Taxon"]<-paste0(all_sum[[xy]][i, "x.Phylum"], ":_", all_sum[[xy]][i, "x.Class"], ":_", 
                                      all_sum[[xy]][i, "x.Order"], ":_",
                                      all_sum[[xy]][i, "x.Family"], ":_", all_sum[[xy]][i, "x.Genus"])
  }
}

taxa_all_s1<-rbind(all_sum[[1]],
                   all_sum[[2]],
                   all_sum[[3]],
                   all_sum[[4]],
                   all_sum[[5]],
                   all_sum[[5]],
                   all_sum[[5]],
                   all_sum[[5]]
)

taxa_all_s2<-rbind(all_sum[[6]],
                   all_sum[[7]],
                   all_sum[[8]],
                   all_sum[[9]],
                   all_sum[[10]],
                   all_sum[[11]],
                   all_sum[[12]],
                   all_sum[[13]]
)
unique_taxa<-unique(c(taxa_all_s1$Taxon, taxa_all_s2$Taxon))

dff<-expand.grid(unique_taxa, Lineages)
colnames(dff)<-c("Taxon", "Lineage")
dff$n_strains_dataset<-NA
dff$rel_ab_data_pct<-NA
dff$Genus<-NA
dff$Family<-NA
dff$Order<-NA
dff$Class<-NA
dff$Phylum<-NA

for (xy in Lineages){
  ps<-subset_samples(ps_3_rel, sample_data(ps_3_rel)[,"Caste_composite"]=="T1 colony Worker")
  ps<-subset_samples(ps, sample_data(ps)[,"Colony"]==xy)
  ps<-prune_taxa(taxa_sums(ps)>0, ps)
  
  for (i in unique(dff$Taxon)){i
    
    # ps2<-prune_taxa
    
    #
    tryCatch({
      t1<-str_split(i, ":_")[[1]][1]
      t2<-str_split(i, ":_")[[1]][2]
      t3<-str_split(i, ":_")[[1]][3]
      t4<-str_split(i, ":_")[[1]][4]
      t5<-str_split(i, ":_")[[1]][5]
      if(t5=="NA"){
        if(t4=="NA"){
          if (t3=="NA"){
            if(t2=="NA"){
              if(t1=="NA"){
                ps_int<-subset_taxa(ps,
                                    is.na(tax_table(ps)[,1]) &
                                      is.na(tax_table(ps)[,2]) &
                                      is.na(tax_table(ps)[,3]) &
                                      is.na(tax_table(ps)[,4]) &
                                      is.na(tax_table(ps)[,5])) 
              }else{
                ps_int<-subset_taxa(ps,
                                    tax_table(ps)[,1] == t1 &
                                      is.na(tax_table(ps)[,2]) &
                                      is.na(tax_table(ps)[,3]) &
                                      is.na(tax_table(ps)[,4]) &
                                      is.na(tax_table(ps)[,5]))
              }}
            else{
              ps_int<-subset_taxa(ps,
                                  tax_table(ps)[,1] == t1 &
                                    tax_table(ps)[,2] == t2 &
                                    is.na(tax_table(ps)[,3]) &
                                    is.na(tax_table(ps)[,4]) &
                                    is.na(tax_table(ps)[,5]))
            }}
          else{
            ps_int<-subset_taxa(ps,
                                tax_table(ps)[,1] == t1 &
                                  tax_table(ps)[,2] == t2 &
                                  tax_table(ps)[,3] == t3 &
                                  is.na(tax_table(ps)[,4]) &
                                  is.na(tax_table(ps)[,5]))
          }}else{
            ps_int<-subset_taxa(ps,
                                tax_table(ps)[,1] == t1 &
                                  tax_table(ps)[,2] == t2 &
                                  tax_table(ps)[,3] == t3 &
                                  tax_table(ps)[,4] == t4 &
                                  is.na(tax_table(ps)[,5] ))
          }}else{
            ps_int<-subset_taxa(ps,
                                tax_table(ps)[,1] == t1 &
                                  tax_table(ps)[,2] == t2 &
                                  tax_table(ps)[,3] == t3 &
                                  tax_table(ps)[,4] == t4 &
                                  tax_table(ps)[,5] == t5 
            )}
      dff[dff$Taxon==i, "Genus"]<-t5
      dff[dff$Taxon==i, "Family"]<-t4
      dff[dff$Taxon==i, "Class"]<-t3
      dff[dff$Taxon==i, "Order"]<-t2
      dff[dff$Taxon==i, "Phylum"]<-t1
      dff[dff$Taxon==i, "n_strains_dataset"]<-ntaxa(ps_int)
      dff[dff$Taxon==i, "rel_ab_data_pct"]<-(sum(taxa_sums(ps_int))/sum(taxa_sums(ps)))*100
      
      print((which(i==dff$Taxon)/length(dff$Taxon)*100))
    }, error = function(e){print(e)})
    #
  }
  
}
head(dff)


## add whether those taxa are VT in thoat lineage ##

head(taxa_Y_s1)

for (i in unique(dff$Taxon)){
  dff[dff$Taxon==i & dff$Lineage == "Yellow", "Transmission"]<-
    dff[dff$Taxon==i & dff$Lineage == "Yellow", "Taxon"] %in% taxa_Y_s1
  dff[dff$Taxon==i & dff$Lineage == "Red", "Transmission"]<-
    dff[dff$Taxon==i & dff$Lineage == "Red", "Taxon"] %in% taxa_R_s1
  dff[dff$Taxon==i & dff$Lineage == "Green", "Transmission"]<-
    dff[dff$Taxon==i & dff$Lineage == "Green", "Taxon"] %in% taxa_G_s1
  dff[dff$Taxon==i & dff$Lineage == "Blue", "Transmission"]<-
    dff[dff$Taxon==i & dff$Lineage == "Blue", "Taxon"] %in% taxa_B_s1
}

summary(dff)

ggplot(dff, aes(y=rel_ab_data_pct, x = Transmission, col = Lineage))+
  geom_jitter()+
  scale_y_log10()

ggplot(dff, aes(y=rel_ab_data_pct, x = Transmission, fill = Lineage))+
  geom_boxplot()+
  scale_y_log10()

ggplot(dff, aes(y=n_strains_dataset, x = Transmission, fill = Lineage))+
  geom_boxplot()+
  scale_y_log10()


kruskal.test(data = dff, Transmission~n_strains_dataset)
kruskal.test(data = dff, Transmission~rel_ab_data_pct)


