#### Venn diagrams Figure 2 ####


library(UpSetR)
library(ComplexUpset)
Lin_Y<-list(
  "Yellow_step_1"=tax_Y_s1,
  "Father_step_1"=tax_F_s1,
  "Yellow_step_2"=tax_Y_s2b,
  "Yellow_step_2_father"=tax_Y_s2Fb
)

Lin_R<-list(
  "Red_step_1"=tax_R_s1,
  "Father_step_1"=tax_F_s1,
  "Red_step_2"=tax_R_s2b,
  "Red_step_2_father"=tax_R_s2Fb
)

Lin_G<-list(
  "Green_step_1"=tax_G_s1,
  "Father_step_1"=tax_F_s1,
  "Green_step_2"=tax_G_s2b,
  "Green_step_2_father"=tax_G_s2Fb
)

Lin_B<-list(
  "Blue_step_1"=tax_B_s1,
  "Father_step_1"=tax_F_s1,
  "Blue_step_2"=tax_B_s2b,
  "Blue_step_2_father"=tax_B_s2Fb
)

step_1<-list(
  "Yellow_step_1"=tax_Y_s1,
  "Red_step_1"=tax_R_s1,
  "Green_step_1"=tax_G_s1,
  "Blue_step_1"=tax_B_s1,
  "Father_step_1"=tax_F_s1
)

step_2b<-list(
  "Yellow_step_2_father"=tax_Y_s2Fb,
  "Red_step_2_father"=tax_R_s2Fb,
  "Green_step_2_father"=tax_G_s2Fb,
  "Blue_step_2_father"=tax_B_s2Fb,
  "Yellow_step_2"=tax_Y_s2b,
  "Red_step_2"=tax_R_s2b,
  "Green_step_2"=tax_G_s2b,
  "Blue_step_2"=tax_B_s2b
)

both_steps<-c(step_1, step_2b)
names(both_steps)<-c(names(step_1), names(step_2b))

# three_steps<-c(step_1, step_2, step_3)
# names(three_steps)<-c(names(step_1), names(step_2), names(step_3))

pb<-UpSetR::upset(data = fromList(Lin_B), nintersects = NA, nsets = 8, keep.order = T, point.size = 4, text.scale = 2.2,
                  sets = names(Lin_B))
py<-UpSetR::upset(data = fromList(Lin_Y), nintersects = NA, nsets = 8, keep.order = T, point.size = 4, text.scale = 2.2,
                  sets = names(Lin_Y))
pr<-UpSetR::upset(data = fromList(Lin_R), nintersects = NA, nsets = 8, keep.order = T, point.size = 4, text.scale = 2.2,
                  sets = names(Lin_R))
pg<-UpSetR::upset(data = fromList(Lin_G), nintersects = NA, nsets = 8, keep.order = T, point.size = 4, text.scale = 2.2,
                  sets = names(Lin_G))

p1<-UpSetR::upset(data = fromList(step_1), nintersects = NA, nsets = 8, keep.order = T, point.size = 4, text.scale = 2.2,
                  sets = names(step_1))
p2b<-UpSetR::upset(data = fromList(step_2b), nintersects = NA, nsets = 8, keep.order = T, point.size = 4, text.scale = 2.2,
                   sets = names(step_2b))

pfull<-UpSetR::upset(data = fromList(both_steps), nintersects = NA, nsets = 13, keep.order = T, point.size = 4, text.scale = 2.2,
                     sets = names(both_steps))
p1
p2b
pfull

pb 
py
pr
pg


# Play with pfull

setting_ints_1<-list(list("Red_step_1"), list("Yellow_step_1"), 
                     list("Green_step_1"), list("Blue_step_1"),
                     list("Father_step_1"),
                     list("Red_step_1", "Yellow_step_1", "Green_step_1", "Blue_step_1", "Father_step_1")
)

setting_ints_1_sup<-list(list("Red_step_1"), list("Yellow_step_1"), 
                         list("Green_step_1"), list("Blue_step_1"),
                         list("Father_step_1"),
                         list("Red_step_1", "Yellow_step_1", "Green_step_1", "Blue_step_1", "Father_step_1")
)

pclean_1<-UpSetR::upset(data = fromList(both_steps), nintersects = 16, nsets = 4, keep.order = T, point.size = 4, text.scale = 2.2,
                        sets = names(both_steps), intersections = setting_ints_1,
                        main.bar.color	= "black",
                        queries = list(list(query = intersects, params = list("Blue_step_1"), color = list("deepskyblue"), active = T, query.name = "Blue-specific taxa"),
                                       list(query = intersects, params = list("Yellow_step_1"), color = list("yellow3"), active = T, query.name = "Yellow-specific taxa"),
                                       list(query = intersects, params = list("Red_step_1"), color = list("brown1"), active = T, query.name = "Red-specific taxa"),
                                       list(query = intersects, params = list("Green_step_1"), color = list("green3"), active = T, query.name = "Green-specific taxa"),
                                       list(query = intersects, params = list("Father_step_1"), color = list("darkgrey"), active = T, query.name = "Blue-male-specific taxa")
                                       # list(query = intersects, params = list("Yellow_step_2_father"), color = list("orange"), active = T, query.name = "Yellow-male-specific taxa"),
                                       # list(query = intersects, params = list("Red_step_2_father"), color = list("dark red"), active = T, query.name = "Red-male-specific taxa"),
                                       # list(query = intersects, params = list("Green_step_2_father"), color = list("dark green"), active = T, query.name = "Green-male-specific taxa")
                        ))
pclean_1

pclean_1_sup<-UpSetR::upset(data = fromList(step_1), #nintersects = 16,
                            nsets = 4, keep.order = T, point.size = 4, text.scale = 2.2,
                            sets = names(step_1), #intersections = both_steps,
                            main.bar.color	= "black", show.numbers = "no", set_size.show = "no",
                            queries = list(list(query = intersects, params = list("Blue_step_1"), color = list("deepskyblue"), active = T, query.name = "Blue-specific taxa"),
                                           list(query = intersects, params = list("Yellow_step_1"), color = list("yellow3"), active = T, query.name = "Yellow-specific taxa"),
                                           list(query = intersects, params = list("Red_step_1"), color = list("brown1"), active = T, query.name = "Red-specific taxa"),
                                           list(query = intersects, params = list("Green_step_1"), color = list("green3"), active = T, query.name = "Green-specific taxa"),
                                           list(query = intersects, params = list("Father_step_1"), color = list("darkgrey"), active = T, query.name = "Blue-male-specific taxa")
                                           # list(query = intersects, params = list("Yellow_step_2_father"), color = list("orange"), active = T, query.name = "Yellow-male-specific taxa"),
                                           # list(query = intersects, params = list("Red_step_2_father"), color = list("dark red"), active = T, query.name = "Red-male-specific taxa"),
                                           # list(query = intersects, params = list("Green_step_2_father"), color = list("dark green"), active = T, query.name = "Green-male-specific taxa")
                            ))
pclean_1_sup # S_fig_venn_step1


setting_ints_2<-list(list("Red_step_2"), list("Yellow_step_2"), 
                     list("Green_step_2"), list("Blue_step_2"),
                     list("Red_step_2_father"), list("Yellow_step_2_father"), 
                     list("Green_step_2_father"), list("Blue_step_2_father"),
                     list("Red_step_2", "Yellow_step_2", "Green_step_2", "Blue_step_2",
                          "Red_step_2_father", "Yellow_step_2_father", "Green_step_2_father", "Blue_step_2_father")
)


pclean_2<-UpSetR::upset(data = fromList(step_2b), nintersects = NA, nsets = 4, keep.order = T, point.size = 4, text.scale = 2.2,
                        sets = names(step_2b), group.by = "degree", main.bar.color	= "black",
                        intersections = setting_ints_2, 
                        query.legend = "bottom", 
                        queries = list(list(query = intersects, params = list("Blue_step_2"), color = list("deepskyblue"), active = T, query.name = "Blue-specific taxa"),
                                       list(query = intersects, params = list("Yellow_step_2"), color = list("yellow3"), active = T, query.name = "Yellow-specific taxa"),
                                       list(query = intersects, params = list("Red_step_2"), color = list("brown1"), active = T, query.name = "Red-specific taxa"),
                                       list(query = intersects, params = list("Green_step_2"), color = list("green3"), active = T, query.name = "Green-specific taxa"),
                                       list(query = intersects, params = list("Blue_step_2_father"), color = list("dark blue"), active = T, query.name = "Blue-male-specific taxa"),
                                       list(query = intersects, params = list("Yellow_step_2_father"), color = list("darkorange3"), active = T, query.name = "Yellow-male-specific taxa"),
                                       list(query = intersects, params = list("Red_step_2_father"), color = list("brown4"), active = T, query.name = "Red-male-specific taxa"),
                                       list(query = intersects, params = list("Green_step_2_father"), color = list("dark green"), active = T, query.name = "Green-male-specific taxa")
                        ))
pclean_2

pclean_2_sup<-UpSetR::upset(data = fromList(step_2b), #nintersects = 16,
                            nsets = 4, keep.order = T, point.size = 4, text.scale = 2.2,
                            sets = names(step_2b), #intersections = both_steps,
                            main.bar.color	= "black", show.numbers = "no", set_size.show = "no",
                            queries = list(list(query = intersects, params = list("Blue_step_2"), color = list("deepskyblue"), active = T, query.name = "Blue-specific taxa"),
                                           list(query = intersects, params = list("Yellow_step_2"), color = list("yellow3"), active = T, query.name = "Yellow-specific taxa"),
                                           list(query = intersects, params = list("Red_step_2"), color = list("brown1"), active = T, query.name = "Red-specific taxa"),
                                           list(query = intersects, params = list("Green_step_2"), color = list("green3"), active = T, query.name = "Green-specific taxa"),
                                           list(query = intersects, params = list("Blue_step_2_father"), color = list("dark blue"), active = T, query.name = "Blue-male-specific taxa"),
                                           list(query = intersects, params = list("Yellow_step_2_father"), color = list("darkorange3"), active = T, query.name = "Yellow-male-specific taxa"),
                                           list(query = intersects, params = list("Red_step_2_father"), color = list("brown4"), active = T, query.name = "Red-male-specific taxa"),
                                           list(query = intersects, params = list("Green_step_2_father"), color = list("dark green"), active = T, query.name = "Green-male-specific taxa")
                            ))
pclean_2_sup

df<-data.frame(matrix(ncol=7, nrow = ntaxa(ps_3_clean)))
colnames(df)<-c("Taxon", "Genus", "Family", "Order", "Class", "Phylum", "Abundance")
df
df$Taxon<-taxa_names(ps_3_clean)

df$B_s1<-df$Taxon %in% tax_B_s1
df$R_s1<-df$Taxon %in% tax_R_s1
df$Y_s1<-df$Taxon %in% tax_Y_s1
df$G_s1<-df$Taxon %in% tax_G_s1
df$F_s1<-df$Taxon %in% tax_F_s1

df$B_s2<-df$Taxon %in% tax_B_s2b 
df$R_s2<-df$Taxon %in% tax_R_s2b 
df$Y_s2<-df$Taxon %in% tax_Y_s2b 
df$G_s2<-df$Taxon %in% tax_G_s2b 

df$B_s2F<-df$Taxon %in% tax_B_s2Fb 
df$R_s2F<-df$Taxon %in% tax_B_s2Fb 
df$Y_s2F<-df$Taxon %in% tax_Y_s2Fb 
df$G_s2F<-df$Taxon %in% tax_G_s2Fb 

df$Y_M<-df$Taxon %in% tax_Y_M
df$R_M<-df$Taxon %in% tax_R_M
df$G_M<-df$Taxon %in% tax_G_M
df$B_M<-df$Taxon %in% tax_B_M

df$Y_D<-df$Taxon %in% tax_Y_D
df$R_D<-df$Taxon %in% tax_R_D
df$G_D<-df$Taxon %in% tax_G_D
df$B_D<-df$Taxon %in% tax_B_D

df$Y_A1<-df$Taxon %in% tax_Y_A1
df$Y_A2<-df$Taxon %in% tax_Y_A2
df$R_A1<-df$Taxon %in% tax_R_A1
df$R_A2<-df$Taxon %in% tax_R_A2
df$G_A1<-df$Taxon %in% tax_G_A1
df$G_A2<-df$Taxon %in% tax_G_A2
df$B_A1<-df$Taxon %in% tax_B_A1
df$B_A2<-df$Taxon %in% tax_B_A2

interesting_cols<-colnames(df)[8:length(colnames(df))]
interesting_cols

for (i in taxa_names(ps_3_clean)){
  df[df$Taxon==i, "Genus"]<-tax_table(ps_3_clean)[i, "Genus"]
  df[df$Taxon==i, "Family"]<-tax_table(ps_3_clean)[i, "Family"]
  df[df$Taxon==i, "Order"]<-tax_table(ps_3_clean)[i, "Order"]
  df[df$Taxon==i, "Class"]<-tax_table(ps_3_clean)[i, "Class"]
  df[df$Taxon==i, "Phylum"]<-tax_table(ps_3_clean)[i, "Phylum"]
  df[df$Taxon==i, "Full_tax"]<-paste0(
    tax_table(ps_3_clean)[i, "Phylum"], ":_",
    tax_table(ps_3_clean)[i, "Class"], ":_",
    tax_table(ps_3_clean)[i, "Order"], ":_",
    tax_table(ps_3_clean)[i, "Family"], ":_",
    tax_table(ps_3_clean)[i, "Genus"]
  )
  print(which(taxa_names(ps_3_clean) %in% i))
}

df0<-head(df)
df0



for ( i in interesting_cols){
  df[df[,i]==T,i]<-1
}

head(df)

## Abundance missing




## Looking for specificity amongst genera

df1<-subset(df, !is.na(df$Phylum))

library(reshape2)
interesting_cols2<-interesting_cols[1:8]
interesting_cols3<-interesting_cols[1:4]
df2<-melt(df1, measure.vars = interesting_cols3)

df3<-subset(df2, value == 1)

ggplot(df3, aes(fill = Phylum, x = Genus))+
  geom_bar(stat="count", position = "dodge")+
  facet_grid(variable~Phylum, scale="free_x", space="free_x")+
  theme_classic()+
  scale_y_log10()+
  theme(#legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 9))



#### Tease specific taxa apart ####
ubi_taxa<-subset(tax_Y_s1,
                 tax_Y_s1 %in% tax_B_s1 &
                   tax_Y_s1 %in% tax_G_s1 &
                   tax_Y_s1 %in% tax_R_s1 &
                   tax_Y_s1 %in% tax_F_s1
)

yel_taxa<-subset(tax_Y_s1,
                 !(tax_Y_s1 %in% tax_B_s1) &
                   !(tax_Y_s1 %in% tax_G_s1) &
                   !(tax_Y_s1 %in% tax_R_s1)
)

red_taxa<-subset(tax_R_s1,
                 !(tax_R_s1 %in% tax_B_s1) &
                   !(tax_R_s1 %in% tax_G_s1) &
                   !(tax_R_s1 %in% tax_Y_s1)
)

green_taxa<-subset(tax_G_s1,
                   !(tax_G_s1 %in% tax_B_s1) &
                     !(tax_G_s1 %in% tax_Y_s1) &
                     !(tax_G_s1 %in% tax_R_s1)
)

blue_taxa<-subset(tax_B_s1,
                  !(tax_B_s1 %in% tax_Y_s1) &
                    !(tax_B_s1 %in% tax_G_s1) &
                    !(tax_B_s1 %in% tax_R_s1)
)
##

df$Ubiquitous<-df$Taxon %in% ubi_taxa
df$Yellow_specific<-df$Taxon %in% yel_taxa
df$Red_specific<-df$Taxon %in% red_taxa
df$Green_specific<-df$Taxon %in% green_taxa
df$Blue_specific<-df$Taxon %in% blue_taxa
df$Yellow_all<-df$Taxon %in% tax_Y_s2b | df$Taxon %in% tax_Y_s2Fb
df$Red_all<-df$Taxon %in% tax_R_s2b | df$Taxon %in% tax_R_s2Fb
df$Green_all<-df$Taxon %in% tax_G_s2b | df$Taxon %in% tax_G_s2Fb
df$Blue_all<-df$Taxon %in% tax_B_s2b | df$Taxon %in% tax_B_s2Fb

df1<-subset(df, !is.na(df$Phylum))
df2<-melt(df1, measure.vars = c("Ubiquitous", "Yellow_specific", "Red_specific", "Green_specific", "Blue_specific",
                                "Yellow_all", "Red_all", "Green_all", "Blue_all"))

df3<-subset(df2, value == T)
dim(df3)

head(df2, 100)

df3 %>%
  # group_by() %>%
  group_by(variable, Genus, Family, Order, Class, Phylum) %>%
  summarise(count = n()) -> df4

df4$count<-df4$count+1

for (i in 1:nrow(df4)){
  if(is.na(df4[i, "Genus"])){
    if(is.na(df4[i, "Family"])){
      if(is.na(df4[i, "Order"])){
        if(is.na(df4[i, "Class"])){
          df4[i, "Genus"]<-paste0(df4[i, "Phylum"], "_Unclassified")
        }else{
          df4[i, "Genus"]<-paste0(df4[i, "Class"], "_Unclassified")
        }
      }else{
        df4[i, "Genus"]<-paste0(df4[i, "Order"], "_Unclassified")
      }
    }else{
      df4[i, "Genus"]<-paste0(df4[i, "Family"], "_Unclassified")
    }
  }
}

df4[df4$Genus=="NA_Unclassified",]
colnames(df4)
unique(df4$Genus)

df4$variable<-factor(df4$variable, levels = c("Ubiquitous", "Yellow_specific", "Red_specific", "Green_specific", "Blue_specific",
                                              "Yellow_all", "Red_all", "Green_all", "Blue_all"
))

# saveRDS(df4, "~/Desktop/Bin_stats/df4_venn_plots_all.rds") # including the all category
# saveRDS(df4, "~/Desktop/Bin_stats/df4_venn_plots.rds") # original in Veronica's Figure 2
# saveRDS(df4, "~/Desktop/Bin_stats/df4_venn_plots_good_filter.rds")

# df4<-readRDS("~/Desktop/Bin_stats/df4_venn_plots_all.rds")
df4b<-subset(df4, variable %in% c("Yellow_all", "Red_all", "Green_all", "Blue_all"))
core_genera<-unique(df4$Genus)


col_df<-subset(df4, variable %in% unique(df4$variable)[2:5])
levels(col_df$variable)<-levels(df4$variable)

ggplot(df4, aes(fill = variable, x = Genus, y = count))+
  geom_bar(stat= "identity")+
  facet_grid(variable~Phylum, scale="free_x", space="free_x")+
  theme_classic()+
  scale_y_log10()+
  scale_fill_manual(values = c("black", "yellow3", "brown1", "green3", "deepskyblue",
                               "yellow3", "brown1", "green3", "deepskyblue"))+
  theme(legend.position = "none",
        strip.text.x = element_text(angle = 90, size = 11),
        # strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+#element_text(angle = 90, size = 8, hjust = 1, vjust = 0.7))+
  geom_hline(yintercept = 1)

sort(unique(df4$Phylum))



## Save Table SX ##

# write.csv2(df4, "~/Desktop/Bin_stats/Sup_data_paper/table_s8.csv")


## Counting phyla ##
phyla <- c("Firmicutes", "Proteobacteria", "Bacteroidetes", "Actinobacteria")
for (i in phyla){
  df<-subset(df4b, Phylum == i)
  n_asvs<-sum(df$count)
  n_genera<-length(unique(df$Genus))
  print(i)
  print(n_asvs)
  print(n_genera)
}

ps_p<-tax_glom(ps3_clean, "Phylum")

plot_bar(ps_p, fill = "Phylum")

## Subsetted plots
valid_genera<-df4 %>%
  group_by(Genus, Family, Order, Class, Phylum) %>%
  transmute(count = sum(count)) %>%
  filter(count > 10)

df5<-subset(df4, df4$Genus %in% unique(valid_genera$Genus))
ggplot(df5, aes(fill = variable, x = Genus, y = count))+
  geom_bar(stat= "identity")+
  facet_grid(variable~Phylum, scale="free_x", space="free_x")+
  theme_classic()+
  scale_y_log10()+
  scale_fill_manual(values = c("black", "yellow3", "brown1", "green3", "deepskyblue",
                               "yellow3", "brown1", "green3", "deepskyblue"))+
  theme(legend.position = "none",
        strip.text.x = element_text(angle = 90, size = 11),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.7))+
  geom_hline(yintercept = 1)








