#### Strain & genus hyper test ####
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggrepel)


###############################
### Get strains in transmission
###############################

all_tax<-list(tax_Y_s1, tax_B_s1, tax_R_s1, tax_G_s1, tax_F_s1,
              tax_Y_s2b, tax_B_s2b, tax_R_s2b, tax_G_s2b,
              tax_Y_s2Fb, tax_B_s2Fb, tax_R_s2Fb, tax_G_s2Fb
)

names(all_tax)<-c("Y_s1", "B_s1", "R_s1", "G_s1", "F_s1",
                  "Y_s2", "B_s2", "R_s2", "G_s2",
                  "Y_s2F", "B_s2F", "R_s2F", "G_s2F"
)

cols<-c("Yellow", "Blue", "Red", "Green", "Father",
        "Yellow", "Blue", "Red", "Green",
        "Yellow", "Blue", "Red", "Green"
)

Sex<-c("M", "M", "M", "M", "F",
       "M", "M", "M", "M",
       "F", "F", "F", "F"
)

Step<-c("1", "1", "1", "1", "1",
        "2", "2", "2", "2",
        "2", "2", "2", "2"
)

Batch<-c("No", "No", "No", "No", "Yes", 
         "Yes", "Yes", "Yes", "Yes",
         "No", "No", "No", "No")


all_suport<-list(ps_Y_M1, ps_B_M1, ps_R_M1, ps_G_M1, ps_F_M,
                 ps_Y_M1, ps_B_M1, ps_R_M1, ps_G_M1,
                 ps_F_M, ps_F_M, ps_F_M, ps_F_M
)

all_2<-list()
for (i in 1:length(all_tax)){ # Count all tax classification for each strain.
  tax<-tax_table(ps_3_clean)[all_tax[[i]]]
  taxa<-tax[, c("Phylum", "Class", "Order", "Family", "Genus")] %>% plyr::count()
  all_2[[i]]<-taxa #Sum_X_sX
  print(paste0(round(((i/length(all_tax))*100),2
  ), " %"))
}

all_support_2<-list() 
for (i in 1:length(all_suport)){ # Count all tax classs for each strain in the support phyloseqs. 
  ps<-all_suport[[i]]
  tax<-tax_table(ps)
  taxa<-tax[, c("Phylum", "Class", "Order", "Family","Genus")] %>% plyr::count()
  taxa$Taxon<-paste0(taxa[i, "x.Phylum"], ":_", taxa[i, "x.Class"], ":_", # Create a taxon column to compare and count strains properly
                     taxa[i, "x.Order"], ":_", taxa[i, "x.Family"], ":_", taxa[i, "x.Genus"])
  all_support_2[[i]]<-taxa # Sum_X_sX
  print(i)
}


fish_tests<-data.frame(matrix(nrow=length(all_2), ncol=1))
colnames(fish_tests)<-"p_val"

for (i in 1:length(all_2)){
  Sum_M<-tax_table(all_suport[[i]])[, c("Phylum", "Class", "Order", "Family","Genus")] %>% plyr::count()
  m<-matrix(ncol=2, nrow=2)
  m[1,1]<-ntaxa(all_suport[[i]])-sum(all_2[[i]][,"freq"]) # Strains not transmitted
  m[1,2]<-sum(all_2[[i]][,"freq"]) # Strains transmitted
  m[2,1]<-nrow(Sum_M)-nrow(all_2[[i]]) # Genera not transmitted
  m[2,2]<-nrow(all_2[[i]]) # Genera transmitted
  m2<-fisher.test(m)
  fish_tests[i, "p_val"]<-m2$p.value
  fish_tests[i, "odds"]<-c(m2$estimate)
  fish_tests[i, "lwr"]<-c(m2$conf.int[1])
  fish_tests[i, "uppr"]<-c(m2$conf.int[2])
  print(i)
}
fish_tests$fish_sig<-fish_tests$p_val<0.05
fish_tests$Group<-names(all_tax)

ggplot(fish_tests, aes(y=p_val, x=odds, label=Group, col=fish_sig))+
  geom_point()+scale_y_log10()+geom_text_repel(size=4)

max(fish_tests$p_val)

##### Run by simulation #####
predict_gen<-data.frame(matrix(ncol=2, nrow=length(all_2)))
colnames(predict_gen)<-c("Random_genera", "VT_genera")

List_of_lists<-list()

set.seed(123)
for(xy in 1:length(all_2)){
  Sum_M<-tax_table(all_suport[[xy]])@.Data[, c("Phylum", "Class", "Order", "Family","Genus")] 
  Sum_M<-data.frame(Sum_M)
  Sum_M$Taxon<-NA
  for (i in 1:nrow(Sum_M)){
    Sum_M[i, "Taxon"]<-paste0(unique(Sum_M[i, "Phylum"]), ":_", unique(Sum_M[i, "Class"]), ":_",
                              unique(Sum_M[i, "Order"]), ":_", unique(Sum_M[i, "Family"]), ":_", unique(Sum_M[i, "Genus"]))
  }
  for (i in 1:nrow(all_2[[xy]])){
    all_2[[xy]][i, "Taxon"]<-paste0(all_2[[xy]][i, "x.Phylum"], ":_", all_2[[xy]][i, "x.Class"], ":_", 
                                    all_2[[xy]][i, "x.Order"], ":_",
                                    all_2[[xy]][i, "x.Family"], ":_", all_2[[xy]][i, "x.Genus"])
  }
  # vec<-Sum_M$Taxon
  vec<-rownames(Sum_M)
  ab_vector<-taxa_sums(all_suport[[xy]]) # Get the strain abundances to weigh
  List_of_lists[[xy]]<-list()
  for (i in 1:1000){
    List_of_lists[[xy]][[i]]<-sample(vec, sum(all_2[[xy]][,"freq"]), replace = F, prob=ab_vector) # Random sampling
  }
  
  rand<-length(unique(Sum_M[List_of_lists[[xy]][[i]], "Taxon"])) # How many genera is there? 
  biol<-length(unique(all_2[[xy]][, "Taxon"])) # How many were there in the first place
  predict_gen[xy, "Random_genera"]<-rand
  predict_gen[xy, "VT_genera"]<-biol
  print(paste0(round(((xy/length(all_2))*100),2
  ), " %"))
}
# saveRDS(List_of_lists, "~/Desktop/Bin_stats/Random_models/List_of_lists.rds")
# saveRDS(predict_gen, "~/Desktop/Bin_stats/Random_models/predict_gen.rds")
List_of_lists<-readRDS("~/Desktop/Bin_stats/Random_models/List_of_lists.rds")
predict_gen<-readRDS("~/Desktop/Bin_stats/Random_models/predict_gen.rds")

summary(List_of_lists)
summary(List_of_lists[[1]])
head(List_of_lists[[1]][[1]])
summary(predict_gen)

Simulations<-c(1:100)
names(all_tax)
new_pred_df<-expand.grid(Simulations, names(all_tax))
colnames(new_pred_df)<-c("Simulation", "Group")

for (xy in 1:length(all_2)){
  Sum_M<-tax_table(all_suport[[xy]])@.Data[, c("Phylum", "Class", "Order", "Family","Genus")] 
  Sum_M<-data.frame(Sum_M)
  Sum_M$Taxon<-NA
  for (i2 in 1:nrow(Sum_M)){
    Sum_M[i2, "Taxon"]<-paste0(unique(Sum_M[i2, "Phylum"]), ":_", unique(Sum_M[i2, "Class"]), ":_",
                               unique(Sum_M[i2, "Order"]), ":_", unique(Sum_M[i2, "Family"]), ":_", unique(Sum_M[i2, "Genus"]))
  }
  for (i in Simulations){
    new_pred_df[new_pred_df$Simulation==i & new_pred_df$Group == names(all_tax)[[xy]], "Genera"]<-
      length(unique(Sum_M[List_of_lists[[xy]][[i]], "Taxon"]))
  }
  print(paste0(round(((xy/length(all_2))*100),2
  ), " %"))
}

new_pred_df$Random_or_Biol<-"Random_abund"

new_our_df<-expand.grid(1, names(all_tax))
colnames(new_our_df)<-c("Simulation", "Group")
new_our_df$Random_or_Biol<-"Biological"


for (i in 1:length(all_2)){
  new_our_df[new_our_df$Group==names(all_tax)[[i]], "Genera"]<-
    predict_gen[i, "VT_genera"]
}

df_pred_biol<-rbind(new_our_df, new_pred_df)
df_pred_biol<-subset(df_pred_biol, Group %in% c("Y_s1", "B_s1", "R_s1", "G_s1",
                                                "Y_s2", "B_s2", "R_s2", "G_s2"))
df_pred_biol$Group<-as.character(df_pred_biol$Group)
df_pred_biol[df_pred_biol$Group=="Y_s1", "Group"]<-"Yellow_step_1"
df_pred_biol[df_pred_biol$Group=="R_s1", "Group"]<-"Red_step_1"
df_pred_biol[df_pred_biol$Group=="G_s1", "Group"]<-"Green_step_1"
df_pred_biol[df_pred_biol$Group=="B_s1", "Group"]<-"Blue_step_1"

df_pred_biol[df_pred_biol$Group=="Y_s2", "Group"]<-"Yellow_step_2"
df_pred_biol[df_pred_biol$Group=="R_s2", "Group"]<-"Red_step_2"
df_pred_biol[df_pred_biol$Group=="G_s2", "Group"]<-"Green_step_2"
df_pred_biol[df_pred_biol$Group=="B_s2", "Group"]<-"Blue_step_2"

df_pred_biol$Group<-factor(df_pred_biol$Group, levels = c(
  "Yellow_step_1", "Red_step_1", "Green_step_1", "Blue_step_1",
  "Yellow_step_2", "Red_step_2", "Green_step_2", "Blue_step_2"
))

ggplot(df_pred_biol, aes(y = Genera, x = Group, fill = Random_or_Biol))+
  geom_point(aes(y = Genera, x = Group, colour = Random_or_Biol, size = 3))+
  geom_violin()+
  # scale_colour_manual(values = c("red", NA))+
  ylim(0, NA)+
  theme_classic()+
  theme(axis.text.x = element_text(size=13, angle = 90),
        axis.text.y = element_text(size=13),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        axis.title = element_text(size=14),
        title = element_text(size=15)
  )+ylab("Transmitted genera")+xlab("Steps")+
  # ggtitle("Both batches: # Transmitted genera")+
  scale_fill_manual(values = c(NA, "#b49242"))+
  scale_color_manual(values = c("#b348c7", "white"))+
  guides(size = F, fill = F, color = F) # 4 - 5


##### Are they more similar than expected by chance?? #####

overlaps<-list()
Simulations<-c(1:100)
names(all_tax)


for (xy in 1:length(all_tax)){
  Sum_M<-tax_table(all_suport[[xy]])@.Data[, c("Phylum", "Class", "Order", "Family","Genus")] 
  Sum_M<-data.frame(Sum_M)
  Sum_M$Taxon<-NA
  for (i2 in 1:nrow(Sum_M)){
    Sum_M[i2, "Taxon"]<-paste0(unique(Sum_M[i2, "Phylum"]), ":_", unique(Sum_M[i2, "Class"]), ":_",
                               unique(Sum_M[i2, "Order"]), ":_", unique(Sum_M[i2, "Family"]), ":_", unique(Sum_M[i2, "Genus"]))
  }
  overlaps[[xy]]<-list()
  for (i in Simulations){
    overlaps[[xy]][[i]]<-unique(Sum_M[List_of_lists[[xy]][[i]], "Taxon"])
  }
  print(xy)
}

summary(overlaps)
summary(overlaps[[13]])
summary(overlaps[[1]][[1]])

names(overlaps)<-names(all_tax)

over_df<-expand_grid(Simulations, names(all_tax))
colnames(over_df)<-c("Simulation", "Group")

for (xy in 1:length(all_2)){
  # over_df[[xy]]<-list()
  for (i in Simulations){
    if (xy == 1){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["B_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["R_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["G_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    if (xy == 2){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["Y_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["R_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["G_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    if (xy == 3){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["B_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["Y_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["G_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    if (xy == 4){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["B_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["R_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["Y_s1"]][[i]])/
        length(overlaps[[xy]][[i]])
    }

    if (xy == 6){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["B_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["R_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["G_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    if (xy == 7){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["Y_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["R_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["G_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    if (xy == 8){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["B_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["Y_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["G_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    if (xy == 9){
      other_1<-sum(overlaps[[xy]][[i]] %in% overlaps[["B_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_2<-sum(overlaps[[xy]][[i]] %in% overlaps[["R_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
      other_3<-sum(overlaps[[xy]][[i]] %in% overlaps[["Y_s2"]][[i]])/
        length(overlaps[[xy]][[i]])
    }
    over_df[over_df$Simulation==i & over_df$Group==names(all_tax)[xy],"Value"]<-mean(other_1, other_2, other_3)
  }
}

summary(over_df)
over_df$Group<-factor(over_df$Group, levels = names(all_tax))
over_df$Random_or_Biological<-"Random_abund"

over_our_df<-expand.grid(1, names(all_tax))
colnames(over_our_df)<-c("Simulation", "Group")
over_our_df$Random_or_Biological<-"Biological"

similar_gen_tru<-list()
for(xy in 1:length(all_2)){
  for (i in 1:nrow(all_2[[xy]])){
    all_2[[xy]][i, "Taxon"]<-paste0(all_2[[xy]][i, "x.Phylum"], ":_", all_2[[xy]][i, "x.Class"], ":_", 
                                    all_2[[xy]][i, "x.Order"], ":_",
                                    all_2[[xy]][i, "x.Family"], ":_", all_2[[xy]][i, "x.Genus"])
  }
  similar_gen_tru[[xy]]<-all_2[[xy]][,"Taxon"]
  print(paste0(round(((xy/length(all_2))*100),2
  ), " %"))
}

similar_gen_tru_u<-list()
for (i in 1:length(similar_gen_tru)){
  similar_gen_tru_u[[i]]<-unique(similar_gen_tru[[i]])
}

names(similar_gen_tru_u)<-names(all_tax)
for (xy in 1:length(similar_gen_tru_u)){
  if (xy == 1){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["B_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["R_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["G_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
  }
  if (xy == 2){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["Y_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["R_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["G_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
  }
  if (xy == 3){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["B_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["Y_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["G_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
  }
  if (xy == 4){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["B_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["R_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["Y_s1"]]) / 
      length(similar_gen_tru_u[[xy]])
  }

  if (xy == 6){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["B_s2"]])/length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["R_s2"]]) /length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["G_s2"]]) / length(similar_gen_tru_u[[xy]])
  }
  if (xy == 7){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["Y_s2"]])/length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["R_s2"]]) /length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["G_s2"]]) / length(similar_gen_tru_u[[xy]])
  }
  if (xy == 8){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["B_s2"]])/length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["Y_s2"]]) /length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["G_s2"]]) / length(similar_gen_tru_u[[xy]])
  }
  if (xy == 9){
    val1<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["B_s2"]])/length(similar_gen_tru_u[[xy]])
    val2<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["R_s2"]]) /length(similar_gen_tru_u[[xy]])
    val3<-sum(similar_gen_tru_u[[xy]] %in% similar_gen_tru_u[["Y_s2"]]) / length(similar_gen_tru_u[[xy]])
  }
  val<-mean(val1, val2, val3)
  over_our_df[over_our_df$Group==names(all_tax)[xy] , "Value"]<-
    val
}

colnames(over_our_df)
colnames(over_df)

overlaps_df<-rbind(over_our_df, over_df)
overlaps_df<-subset(overlaps_df, Group %in% c("Y_s1", "B_s1", "R_s1", "G_s1",
                                              "Y_s2", "B_s2", "R_s2", "G_s2"))
overlaps_df$Group<-as.character(overlaps_df$Group)
overlaps_df$Value<-overlaps_df$Value*100
ggplot(overlaps_df, aes(y = Value, x = Group, fill = Random_or_Biological))+
  geom_violin()+
  geom_point(data= subset(overlaps_df, Random_or_Biological == "Biological"),
             aes(y = Value, x = Group, colour = Random_or_Biological, size = 3))+
  scale_colour_manual(values = c("red", NA))+
  ylim(0, NA)+
  theme_classic()+
  theme(axis.text.x = element_text(size=13, angle = 90),
        axis.text.y = element_text(size=13),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        axis.title = element_text(size=14),
        title = element_text(size=15)
  )+ylab("")+xlab("Steps")+
  scale_fill_manual(values = c(NA, "#b49242"))+
  scale_color_manual(values = c("#b348c7", NA))+
  guides(size = F, fill = F, color = F) # 4 - 5

for (i in 1:13){
  a<-sum(
    overlaps_df[overlaps_df$Group==names(all_tax)[i] & 
                  overlaps_df$Random_or_Biological=="Random_abund", "Value"] >
      overlaps_df[overlaps_df$Group==names(all_tax)[i] & 
                    overlaps_df$Random_or_Biological=="Biological", "Value"]
  )
  print(names(all_tax)[i])
  print(1-(a/100))
}




new_pred_df<-expand.grid(Simulations, names(all_tax))
colnames(new_pred_df)<-c("Simulation", "Group")

for (xy in 1:length(all_2)){
  for (i in Simulations){
    new_pred_df[new_pred_df$Simulation==i & new_pred_df$Group == names(all_tax)[[xy]], "Genera"]<-
      length(unique(Sum_M[List_of_lists[[xy]][[i]], "Taxon"]))
  }
  print(paste0(round(((xy/length(all_2))*100),2
  ), " %"))
}





