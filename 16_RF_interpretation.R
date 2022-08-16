#### Sergio SHAP values ####
library(tidyverse)
library(phyloseq)
models<-list()

# path<-"~/Desktop/Bin_stats/Complete_shap/"
path<-"~/Desktop/Bin_stats/Sup_data_paper/Files_s1/"
model_names<-list.files(path)
model_names_shap<-subset(model_names, str_detect(model_names, "_shap.tsv"))
model_names_data<-subset(model_names, str_detect(model_names, "_data.tsv"))

models_shap<-list()
models_data<-list()

for ( i in 1:length(model_names_data)){
  tb_data<-read.table(paste0(path, model_names_data[i]), header = F)
  tb_shap<-read.table(paste0(path, model_names_shap[i]), header = T)
  models_shap[[i]]<-tb_shap
  colnames(tb_data)<-colnames(tb_shap)
  models_data[[i]]<-tb_data
}

## For binary variables
bin_vars<-colnames(models_shap[[i]])
num_vars<-c("Abundance_in_origin", "Netw_degree", "Netw_degree_pos", "Netw_degree_neg")
bin_vars<-subset(bin_vars, !(bin_vars %in% num_vars))
bin_vars
vars<-c(bin_vars, num_vars)

summary_list<-list()
df<-data.frame(matrix(ncol = 2, nrow = length(vars)))
colnames(df)<-c("Variable", "Value")
df$Variable<-vars


for (i in 1:length(model_names_data)){
  for (col in vars){
    if(col %in% bin_vars){
      rows_interest<-models_shap[[i]][,col]==0
      mn_1<-mean(models_data[[i]][rows_interest, col])
      rows_interest_2<-models_shap[[i]][,col]==1
      mn_2<-mean(models_data[[i]][rows_interest_2, col])
      df[df$Variable==col, "Value"]<-mn_2-mn_1
      df[df$Variable==col, "Mean_1"]<-mn_2
      df[df$Variable==col, "Mean_2"]<-mn_1  
    }
    if(col %in% num_vars){
      median_test<-median(models_shap[[i]][, col])
      rows_interest<-models_shap[[i]][,col]<=median_test
      mn_1<-mean(models_data[[i]][rows_interest, col])
      rows_interest_2<-models_shap[[i]][,col]>median_test
      mn_2<-mean(models_data[[i]][rows_interest_2, col])
      df[df$Variable==col, "Value"]<-mn_2-mn_1
      df[df$Variable==col, "Mean_1"]<-mn_2
      df[df$Variable==col, "Mean_2"]<-mn_1
    }
    if(col == "Step"){
      rows_interest<-models_shap[[i]][,col]==1
      mn_1<-mean(models_data[[i]][rows_interest, col])
      rows_interest_2<-models_shap[[i]][,col]==2
      mn_2<-mean(models_data[[i]][rows_interest_2, col])
      df[df$Variable==col, "Value"]<-mn_2-mn_1    
      df[df$Variable==col, "Mean_1"]<-mn_2
      df[df$Variable==col, "Mean_2"]<-mn_1
    }
  }
  df$Model<-model_names_data[i]
  summary_list[[i]]<-df
  print(i)
}
print(names(summary_list)<-model_names_data)
summary(summary_list)


summary_list_1<-subset(summary_list, c(T, F, F, F, F, F, F, F, F, F, F, F))
summary_list_1<-summary_list[[1]]
summary(summary_list_1)
summary_list_2<-subset(summary_list, c(T, F, T, F, F, F, F, F, F, F, F, F))
summary_list_2<-summary_list[c(1,3)]
summary(summary_list_2)
summary_list_3<-subset(summary_list, c(T, F, T, F, T, F, T, F, T, F, T, F))
summary_list_3<-summary_list[seq(2, 12, 2)]

summary_list_raw<-subset(summary_list, !str_detect(names(summary_list), "new"))
summary_list_raw<-summary_list[seq(1, 11, 2)]
summary_list_norm<-subset(summary_list, str_detect(names(summary_list), "new"))
summary_list_norm<-summary_list[seq(2, 12, 2)]

data_raw<-do.call(rbind, summary_list_raw)
ggplot(data_raw, aes(y=Value, x=reorder(Variable, Value)))+
  geom_bar(stat = "identity")+facet_grid(~Model)+
  coord_flip()+
  geom_point(aes(y=Mean_1, x=reorder(Variable, Value)), col = "red")+
  geom_point(aes(y=Mean_2, x=reorder(Variable, Value)), col = "blue")

data_norm<-do.call(rbind, summary_list_norm)
ggplot(data_norm, aes(y=Value, x=reorder(Variable, Value)))+
  geom_bar(stat = "identity")+facet_grid(~Model)+
  coord_flip()

data_norm_2<-subset(data_raw, Variable== "Alphaproteobacteria")
ggplot(data_norm_2, aes(y=Value, x=reorder(Variable, Value)))+
  geom_bar(stat = "identity")+facet_grid(~Model)+
  coord_flip()


#####################
###### Normalize ####
#####################

ratio_all<-rbind(summary_list_norm[[1]], summary_list_norm[[3]], summary_list_norm[[5]])
ratio_all$Ratio<-c(rbind(summary_list_norm[[1]]$Value-summary_list_norm[[2]]$Value),
                   rbind(summary_list_norm[[3]]$Value-summary_list_norm[[4]]$Value),
                   rbind(summary_list_norm[[5]]$Value-summary_list_norm[[6]]$Value))
ratio_all$Value<-c(rbind(summary_list_norm[[1]]$Value),
                   rbind(summary_list_norm[[3]]$Value),
                   rbind(summary_list_norm[[5]]$Value))
ratio_all$Value_random<-c(rbind(0+summary_list_norm[[2]]$Value),
                          rbind(0+summary_list_norm[[4]]$Value),
                          rbind(0+summary_list_norm[[6]]$Value))



threshold<-0.015

ratio_2<-subset(ratio_all, Ratio > threshold |
                  Ratio < -threshold )
vdata<-ratio_2 %>%
  group_by(Model) %>%
  summarise(mean_biol = mean(Value), 
            mean_random = mean(Value_random))

ggplot(ratio_2)+
  geom_bar(aes(y=Ratio, x=reorder(Variable, Ratio)), stat = "identity")+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text (size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text( size = 14))+
  geom_point(aes(y=Value,  x=reorder(Variable, Ratio)), col = "red")+
  geom_point(aes(y=Value_random,  x=reorder(Variable, Ratio)), col = "blue")+
  facet_grid(~Model)

write.csv(ratio_all, "~/Desktop/Table_s11.csv")



### Highlight DA taxa ####
Results_s1<-readRDS( "~/Desktop/Bin_stats/Aldex_s1.rds")

Results_s1$Taxon<-paste0(Results_s1$Phylum, "._", Results_s1$Class, "._", Results_s1$Order, "._", Results_s1$Family, "._", Results_s1$Genus)

ratio_all[ratio_all$Variable %in% Results_s1[Results_s1$Alate_estimate>0, "Taxon"] &
            ratio_all$Model == "F_step1_new_data.tsv", "DA"]<- "Enriched_in_alates"
ratio_all[ratio_all$Variable %in% Results_s1[Results_s1$Alate_estimate<0, "Taxon"] &
            ratio_all$Model == "F_step1_new_data.tsv", "DA"]<- "Decreased_in_alates"

ggplot(ratio_all, aes(y=Value, x=reorder(Variable, Value), fill = DA))+
  geom_bar(stat = "identity")+facet_grid(~Model)+
  coord_flip()

##

Results_s2<-readRDS( "~/Desktop/Bin_stats/Aldex_s2.rds")
Results_s2$Taxon<-paste0(Results_s2$Phylum, "._", Results_s2$Class, "._", Results_s2$Order, "._", Results_s2$Family, "._", Results_s2$Genus)

ratio_all[ratio_all$Variable %in% Results_s2[Results_s2$Daughter_estimate>0, "Taxon"] &
            ratio_all$Model == "J_step2_new_data.tsv", "DA"]<- "Enriched_in_daughters"
ratio_all[ratio_all$Variable %in% Results_s2[Results_s2$Daughter_estimate<0, "Taxon"] &
            ratio_all$Model == "J_step2_new_data.tsv", "DA"]<- "Decreased_in_daughters"


threshold<-0.025

ratio_2<-ratio_all

for (i in 1:nrow(ratio_2)){
  if(ratio_2[i, "Model"] == "B_Complete_new_data.tsv"){
    ratio_2[i, "Model"]<-"Both_steps"
  }
  if(ratio_2[i, "Model"] == "F_step1_new_data.tsv"){
    ratio_2[i, "Model"]<-"Step_1"
  }
  if(ratio_2[i, "Model"] == "J_step2_new_data.tsv"){
    ratio_2[i, "Model"]<-"Step_2"
  }
  length_now<-length(str_split(ratio_2[i, "Variable"], "\\._")[[1]])
  if(length_now>1){
    ratio_2[i, "Variable"]<-
      paste0(str_split(ratio_2[i, "Variable"], "\\._")[[1]][(length_now-1)], ":_",
             str_split(ratio_2[i, "Variable"], "\\._")[[1]][length_now])  
  }else{
    ratio_2[i, "Variable"]<-
      paste0(str_split(ratio_2[i, "Variable"], "\\._")[[1]][(length_now)])
  }
}


ggplot(ratio_2, aes(y=Ratio, x=reorder(Variable, Ratio), fill = DA))+
  geom_bar(stat = "identity")+facet_grid(~Model)+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text (size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)
  )+
  geom_point(aes(y=Value,  x=reorder(Variable, Ratio)), color = "red")+
  geom_point(aes(y=Value_random,  x=reorder(Variable, Ratio)), color = "blue")+
  xlab(NULL)

subset(ratio_2, DA == "Enriched_in_daughters")$Variable
##



###########################################
########### Split by categories ###########
###########################################
Variables<-unique(ratio_2$Variable)

Phyla <- c(unique(tax_table(ps_3_clean)[,"Phylum"]))
Classes <- c(unique(tax_table(ps_3_clean)[,"Class"]))
Orders <- c(unique(tax_table(ps_3_clean)[,"Order"]))
Families <- c(unique(tax_table(ps_3_clean)[,"Family"]))
Genera <- unique(paste0(c(tax_table(ps_3_clean)[,"Family"]), ":_", c(tax_table(ps_3_clean)[,"Genus"])))
Nets <- c("Netw_degree_pos", "Netw_degree_neg", "Netw_degree")
Colonies<- c("Yellow", "Red", "Green", "Blue")
Steps<-c("Step")
Abundance<- c("Abundance_in_origin")

Others<-c(Nets, Colonies, Steps, Abundance)

Variables_up<-subset(Variables, !(Variables %in% Phyla))
Variables_up<-subset(Variables_up, !(Variables_up %in% Classes))
Variables_up<-subset(Variables_up, !(Variables_up %in% Orders))
Variables_up<-subset(Variables_up, !(Variables_up %in% Families))
Variables_up<-subset(Variables_up, !(Variables_up %in% Genera))
Variables_up<-subset(Variables_up, !(Variables_up %in% Nets))
Variables_up<-subset(Variables_up, !(Variables_up %in% Colonies))
Variables_up<-subset(Variables_up, !(Variables_up %in% Steps))
Variables_up

sort(Genera)[130:160]
Variables_up

ratio_2$Model<-as.character(ratio_2$Model)
ratio_2[ratio_2$Model=="B_Complete_new_data.tsv", "Model"]<-"Both steps"
ratio_2[ratio_2$Model=="F_step1_new_data.tsv", "Model"]<-"Step 1"
ratio_2[ratio_2$Model=="J_step2_new_data.tsv", "Model"]<-"Step 2"


groups<-list(Phyla, Classes, Orders, Families, Genera, Others)
p<-list()
for (i in 1:length(groups)){
  df010<-subset(ratio_2, Variable %in% groups[[i]])
  p[[i]]<-ggplot( df010, aes(y=Ratio, x=reorder(Variable, Ratio)))+
    geom_bar(stat = "identity", fill = NA, col = "black")+facet_grid(~Model)+
    coord_flip()+
    theme_classic()+
    theme(axis.text.y = element_text (size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 90),
          axis.title.x = element_blank())+
    # ylim(c(-0.3, 0.5))+
    geom_point(aes(y=Value,  x=reorder(Variable, Ratio)), color = "red", size = 3)+
    geom_point(aes(y=Value_random,  x=reorder(Variable, Ratio)), color = "blue", size = 3)+
    xlab(NULL)
  print(groups[i])
  print(mean(abs(df010$Value)))
}
library(gridExtra)
do.call(grid.arrange, p)
summary(p)

df010<-subset(ratio_2, Variable %in% Genera)
df010[is.na(df010$DA), "DA"]<-"Not enriched"
df010[df010$DA=="Enriched_in_daughters","DA"]<-NA
df010$DA<-as.character(df010$DA)
df010[is.na(df010$DA), "DA"]<-"Not enriched"
p[[5]]<-
  ggplot( df010, aes(y=Ratio, x=reorder(Variable, Ratio), fill = DA))+
  geom_bar(stat = "identity", col = "black")+facet_grid(~Model)+
  scale_fill_manual(values = c("#b348c7", NA))+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_text (size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x = element_blank())+
  # ylim(c(-0.3, 0.5))+
  geom_point(aes(y=Value,  x=reorder(Variable, Ratio)), color = "red", size = 3)+
  geom_point(aes(y=Value_random,  x=reorder(Variable, Ratio)), color = "blue", size = 3)+
  xlab(NULL)

p[[1]] # 8 - 4
p[[2]] # 8 - 6
p[[3]] # 8 - 6
p[[4]] # 8 - 7
p[[5]] # 13  -8
p[[6]] # 8 - 4


max(ratio_2$Value)
min(ratio_2$Value)

max(ratio_2$Value_random)
min(ratio_2$Value_random)

max(ratio_2$Ratio)
min(ratio_2$Ratio)

### How importnt are taxa?

df_taxa<-subset(ratio_2, Variable %in% c(Phyla, Classes, Orders, Families, Genera))
df_taxa$Variable
mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value"]))/
  mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value_random"]))

mean(abs(df_taxa[df_taxa$Model=="Step_1", "Value"]))/
  mean(abs(df_taxa[df_taxa$Model=="Step_1", "Value_random"]))

mean(abs(df_taxa[df_taxa$Model=="Step_2", "Value"]))/
  mean(abs(df_taxa[df_taxa$Model=="Step_2", "Value_random"]))

tax_levels<-list(Phyla, Classes, Orders, Families, Genera)
vec1<-c()
vec2<-c()
vec3<-c()
for (i in 1:length(tax_levels)){
  df_taxa<-subset(ratio_2, Variable %in% tax_levels[[i]])
  print(mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value"]))/
          mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value_random"])))
  vec1[i]<-(mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value"]))-
              mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value_random"])))
  vec2[i]<-mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value"]))
  vec3[i]<-mean(abs(df_taxa[df_taxa$Model=="Both_steps", "Value_random"]))
}


dfp<-data.frame(Taxonomic_rank = c("Phylum", "Class", "Order", "Family", "Genus"),
                Predictive_power_ratio = vec1,
                Predictive_power_biological = vec2,
                Predictive_power_random = vec3
)
dfp$Taxonomic_rank<-factor(dfp$Taxonomic_rank, levels =  c("Genus", "Family", "Order", "Class", "Phylum"))
ggplot()+
  geom_bar(data = dfp, aes(y = Predictive_power_ratio, x = Taxonomic_rank), fill = "white", color = "black", stat = "identity")+
  geom_point(data = dfp, aes(y= Predictive_power_biological, x = Taxonomic_rank), size = 4, color = "red")+
  geom_point(data = dfp, aes(y= Predictive_power_random, x = Taxonomic_rank), size = 4, color = "blue")+
  ylim(0, NA)+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))+
  xlab("Taxonomic rank")+
  ylab("Predictive power")# 4-4


