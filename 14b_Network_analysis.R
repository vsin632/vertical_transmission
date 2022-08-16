### Networks both batches ####

ps_3_clean

ps_D<-subset_samples(ps_3_clean, sample_data(ps_3_clean)[,"Caste_composite"]=="T1 colony Worker")
ps_D<-prune_taxa(taxa_sums(ps_D)>200, ps_D)
ps_D

vt_cols<-c("#B6A971FF", "#00204DFF")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

library(tidyverse)
library(microbiome)
library(microbiomeutilities)
library(phyloseq)
library(ggplot2)
library(network)
library(igraph)
library(intergraph)
library(GGally)
library(ggnet)
library(ggnetwork)
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling

net_beta<-readRDS("~/Desktop/Bin_stats/both_nets/net_D_200_beta.rds")
net_refit<-readRDS("~/Desktop/Bin_stats/both_nets/net_D_200_refit.rds")

rownames(net_beta)<-colnames(net_beta)<-rownames(net_refit)<-colnames(net_refit)<-taxa_names(ps_D)

vsize <- log2(apply(otu_table(ps_D), 2, mean)) # add log abundance as properties of vertex/nodes.
net_ig <- graph.adjacency(net_beta, mode='undirected', weighted = TRUE)
coords.fdr = layout_with_fr(net_ig)
E(net_ig)[weight > 0]$color<-"steelblue" #now color the edges based on their values positive is steelblue
E(net_ig)[weight < 0]$color<-"orange"  #now color the edges based on their values
E(net_ig)[weight>0]
E(net_ig)[weight<0]

netw <- asNetwork(net_ig)
network::set.edge.attribute(netw, "color", ifelse(netw %e% "weight" > 0, "grey", "orange"))

netw %v% "vertex.names"
netw %v% "VT" <- netw %v% "vertex.names" %in% c(tax_Y_s2b, tax_R_s2b, tax_G_s2b, tax_B_s2b,
                                                tax_Y_s2Fb, tax_R_s2Fb, tax_G_s2Fb, tax_B_s2Fb)
netw %v% "nodesize" <- vsize-3
netw %v% "degree_2" <- netw %v% "degree"


stl.mb <- degree.distribution(net_ig)
plot(0:(length(stl.mb)-1), stl.mb, type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")


ig.mod <- graph.adjacency(net_refit, mode='undirected', add.rownames = TRUE)
modules =cluster_fast_greedy(ig.mod)
print(modules)
modularity(modules)
V(ig.mod)$color=modules$membership
plot(ig.mod, col = modules, vertex.size = 4, vertex.label = NA)

net.mod<-asNetwork(ig.mod)

net.mod %v% "membership" <- modules$membership


#############################
#### Taxon characteristics

hubs<-igraph::degree(ig.mod, v = V(ig.mod))
author<-authority_score(ig.mod, scale = F)$vector
betwenness<-betweenness(net.mod, gmode = "graph")
edge_btw<-edge_betweenness(ig.mod, directed = F)
centr_btw<-centr_betw(ig.mod, directed = F)$res
degree<-degree(net.mod, gmode = "graph")
closeness<-igraph::closeness(ig.mod, mode = "all")
eig_centr<-eigen_centrality(ig.mod, directed = F, scale = F)$vector
centr_clo<-centr_clo(ig.mod, mode = "all", normalized = F)$res
core <- coreness(ig.mod, mode="all")
df<-data.frame(Degree = hubs, Authroity = author, Betweeness = betwenness, Centrality = centr_btw,
               Closeness = closeness, Eig_centrality = eig_centr, Clo_centrality = centr_clo,
               Coreness = core, Abundance = vsize)

############################
#### Network attributes ####
############################

ass<-assortativity(ig.mod, V(ig.mod), directed = F)
ass_nom<-assortativity_nominal(ig.mod, directed = F, types = V(ig.mod))
ass_degr<-assortativity_degree(ig.mod, directed = F)

V(net_ig)$Genus <- paste0(tax_table(ps_D)[V(net_ig)$name,"Family"], ":_",
                          tax_table(ps_D)[V(net_ig)$name,"Genus"]
)

V(net_ig)$Genus[V(net_ig)$Genus=="NA:_NA"]<-"AD3:_Unclassified"


V(net_ig)$random_VT <- sample(x = c(0, 1), replace = T, size = ntaxa(ps_D))
V(net_ig)$VT_s2 <- V(net_ig)$name %in% c(tax_Y_s2b, tax_R_s2b, tax_G_s2b, tax_B_s2b,
                                         tax_Y_s2Fb, tax_R_s2Fb, tax_G_s2Fb, tax_B_s2Fb)
V(net_ig)$VT_s1 <- V(net_ig)$name %in% c(tax_Y_s1, tax_R_s1, tax_G_s1, tax_B_s1, tax_F_s1)

random_strains<-readRDS("~/Desktop/Bin_stats/Final_Results.rds")

random_strains_prob_s2<-c(random_strains[[6]], random_strains[[7]],
                          random_strains[[8]], random_strains[[9]],
                          random_strains[[10]], random_strains[[11]],
                          random_strains[[12]], random_strains[[13]])
random_strains_prob_s1<-c(random_strains[[1]], random_strains[[2]],
                          random_strains[[3]], random_strains[[4]],
                          random_strains[[5]])

V(net_ig)$VT_prob_random_s1 <- V(net_ig)$name %in% random_strains_prob_s1
V(net_ig)$VT_prob_random_s2 <- V(net_ig)$name %in% random_strains_prob_s2

assortativity(net_ig, V(net_ig)$random_VT, directed=F)
assortativity(net_ig, V(net_ig)$VT_s1, directed=F)
assortativity(net_ig, V(net_ig)$VT_s2, directed=F)
assortativity(net_ig, V(net_ig)$VT_prob_random_s1, directed=F)
assortativity(net_ig, V(net_ig)$VT_prob_random_s2, directed=F)

##################################
#### Subset pos and neg edges ####
##################################
summary(net_ig)
V(net_ig)$degree<-vsize

pos_net<-subgraph.edges(net_ig, eids = E(net_ig)[weight>0], delete.vertices = F)
summary(pos_net)

neg_net<-subgraph.edges(net_ig, eids = E(net_ig)[weight<0], delete.vertices = F)
summary(neg_net)

pos_net_ig<-asNetwork(pos_net)
neg_net_ig<-asNetwork(neg_net)

hubs_pos<-igraph::degree(pos_net, v = V(pos_net))
author_pos<-authority_score(pos_net, scale = F)$vector
betwenness_pos<-betweenness(pos_net_ig, gmode = "graph")
edge_btw_pos<-edge_betweenness(pos_net, directed = F)
centr_btw_pos<-centr_betw(pos_net, directed = F)$res
degree_pos<-degree(pos_net_ig, gmode = "graph")
closeness_pos<-igraph::closeness(pos_net, mode = "all") 
eig_centr_pos<-eigen_centrality(pos_net, directed = F, scale = F)$vector
centr_clo_pos<-centr_clo(pos_net, mode = "all", normalized = F)$res
core_pos <- coreness(pos_net, mode="all")
df_pos00<-data.frame(Degree = hubs_pos, Authroity = author_pos, Betweeness = betwenness_pos, Centrality = centr_btw_pos,
               Closeness = closeness_pos, Eig_centrality_pos = eig_centr_pos, Clo_centrality = centr_clo_pos,
               Coreness = core_pos, Abundance = vsize)


hubs_neg<-igraph::degree(neg_net, v = V(neg_net))
author_neg<-authority_score(neg_net, scale = F)$vector
betwenness_neg<-betweenness(neg_net_ig, gmode = "graph")
centr_btw_neg<-centr_betw(neg_net, directed = F)$res
degree_neg<-degree(neg_net_ig, gmode = "graph")
closeness_neg<-igraph::closeness(neg_net, mode = "all")
eig_centr_neg<-eigen_centrality(neg_net, directed = F, scale = F)$vector
centr_clo_neg<-centr_clo(neg_net, mode = "all", normalized = F)$res
core_neg <- coreness(neg_net, mode="all")
df_neg00<-data.frame(Degree = hubs_neg, Authroity = author_neg, Betweeness = betwenness_neg, Centrality = centr_btw_neg,
               Eig_centrality = eig_centr_neg, 
               Coreness = core_neg, Abundance = vsize)

##############################
#### Replication #############
##############################
List_of_lists<-readRDS("~/Desktop/Bin_stats/Random_models/List_of_lists.rds")

reps<-length(List_of_lists[[1]])
Methods<-c("Random", "Random_abund", "Biological")
Steps <- c("1", "2")
Associations <- c("Positive", "Negative")

df_as<-expand.grid(1:reps, Methods, Steps, Associations)
df_as$Assortativity<-NA
colnames(df_as)<-c("Replicate", "Method", "Step", "Association", "Assortativity")
df_as<-subset(df_as, Method != "Random" | Step != "2")
df_as

table(df_as$Method, df_as$Step, df_as$Association)

random_strain_list_s1<-list()
random_strain_list_s2<-list()
strains<-list()
for (i in 1:reps){
  random_strain_list_s1[[i]]<-c(List_of_lists[[1]][[i]], List_of_lists[[2]][[i]],
                                List_of_lists[[3]][[i]], List_of_lists[[4]][[i]],
                                List_of_lists[[5]][[i]])
  random_strain_list_s2[[i]]<-c(List_of_lists[[6]][[i]], List_of_lists[[7]][[i]],
                                List_of_lists[[8]][[i]], List_of_lists[[9]][[i]],
                                List_of_lists[[10]][[i]], List_of_lists[[11]][[i]],
                                List_of_lists[[12]][[i]], List_of_lists[[13]][[i]])
  
  V(pos_net)$random_strain_list_s1<- V(pos_net)$name %in% random_strain_list_s1[[i]]
  df_as[df_as$Replicate == i & df_as$Association == "Positive" & df_as$Step=="1" & df_as$Method == "Random_abund", "Assortativity"]<-
    assortativity(pos_net, V(pos_net)$random_strain_list_s1, directed=F)
  
  V(neg_net)$random_strain_list_s1<- V(neg_net)$name %in% random_strain_list_s1[[i]]
  df_as[df_as$Replicate == i & df_as$Association == "Negative" & df_as$Step=="1" & df_as$Method == "Random_abund", "Assortativity"]<-
    assortativity(neg_net, V(neg_net)$random_strain_list_s1, directed=F)
  
  V(pos_net)$random_strain_list_s2<- V(pos_net)$name %in% random_strain_list_s2[[i]]
  df_as[df_as$Replicate == i & df_as$Association == "Positive" & df_as$Step=="2" & df_as$Method == "Random_abund", "Assortativity"]<-
    assortativity(pos_net, V(pos_net)$random_strain_list_s2, directed=F)
  
  V(neg_net)$random_strain_list_s2<- V(neg_net)$name %in% random_strain_list_s2[[i]]
  df_as[df_as$Replicate == i & df_as$Association == "Negative" & df_as$Step=="2" & df_as$Method == "Random_abund", "Assortativity"]<-
    assortativity(neg_net, V(neg_net)$random_strain_list_s2, directed=F)
  
  ## Add purely random sampling
  
  strains[[i]]<-sample(x = c(0, 1), replace = T, size = ntaxa(ps_D))
  V(pos_net)$random_VT<-V(neg_net)$random_VT<-strains[[i]]
  df_as[df_as$Rep == i & df_as$Method == "Random", "Assortativity"]<-assortativity(neg_net, V(neg_net)$random_VT, directed=F)
  
  if (i %% 100 == 0){
    print(paste0(round(i*100/reps, 4), " %"))
  }
}

biol_ass_pos<-assortativity(pos_net, V(pos_net)$VT_s2, directed = F)
biol_ass_pos_1<-assortativity(pos_net, V(pos_net)$VT_s1, directed = F)

biol_ass_neg<-assortativity(neg_net, V(neg_net)$VT_s2, directed = F)
biol_ass_neg_1<-assortativity(neg_net, V(neg_net)$VT_s1, directed = F)

data_fr<-data.frame(matrix(ncol = 2, nrow = 2))
colnames(data_fr)<-c("Association", "Assortativity")
data_fr$Association<-c("Positive", "Negative")
data_fr$Assortativity<-c(biol_ass_pos, biol_ass_neg)
data_fr

df_as<-subset(df_as, Method == "Random_abund" & Step == 1)
ggplot(data = df_as)+
  facet_wrap(~Association)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
  geom_violin(aes(y = Assortativity, x = Method, fill = Step))+
  geom_hline(data = data_fr, aes(yintercept = Assortativity), col = "#b348c7", size = 1.5)+
  theme_classic()+
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
  )+
  scale_fill_manual(values = "#b49242")+
  guides(fill = "none")  # 4 - 5

df_as %>%
  group_by(Association) %>%
  mutate(meen = mean(Assortativity),
         siid = sd(Assortativity)) %>%
  select(Association, meen, siid) %>%
  unique()
  

df_as_rabund_pos<-subset(df_as, Method == "Random_abund" & Association == "Positive")
sum(df_as_rabund_pos$Assortativity>biol_ass_pos)/length(df_as_rabund_pos$Assortativity) # p value = 0
df_as_rabund_neg<-subset(df_as, Method == "Random_abund" & Association == "Negative")
sum(df_as_rabund_neg$Assortativity<biol_ass_neg)/length(df_as_rabund_neg$Assortativity) # p value = 0.512




##################################################
######### Bootstrapping our own distribution #####
##################################################

reps_b<-10
vs<-length(V(net_ig))
b1<-round(vs/reps_b)

df_pos<-data.frame(matrix(ncol = 2, nrow =vs))
colnames(df_pos)<-c("Boot_iteration", "Vertex")
df_pos$Boot_iteration<-c(rep(1:10, b1), 11, 11)

df_neg<-data.frame(matrix(ncol = 2, nrow =vs))
colnames(df_neg)<-c("Boot_iteration", "Vertex")
df_neg$Boot_iteration<-c(rep(1:10, b1), 11, 11)

df_as<-subset(df_as, Method == "Biological" | Replicate < 11)

for (i in 1:reps_b){
  vers<-V(pos_net)$name
  vers<-subset(vers, !(vers %in% df_pos$Vertex))
  df_pos[df_pos$Boot_iteration==i, "Vertex"]<-sample(vers, b1, F)
  V(pos_net)$boots<-V(pos_net)$name %in% df_pos[df_pos$Boot_iteration==i, "Vertex"]
  sub_net_pos<-induced_subgraph(pos_net, V(pos_net)$boots==TRUE)#, delete.vertices = F)
  
  df_as[df_as$Method == "Biological" & df_as$Replicate == i & df_as$Association == "Positive", "Assortativity"]<-
    assortativity(sub_net_pos, V(sub_net_pos)$VT, directed = F)
  
  vers<-V(neg_net)$name
  vers<-subset(vers, !(vers %in% df_neg$Vertex))
  df_neg[df_neg$Boot_iteration==i, "Vertex"]<-sample(vers, b1, F)
  V(neg_net)$boots<-V(neg_net)$name %in% df_neg[df_neg$Boot_iteration==i, "Vertex"]
  sub_net_neg<-induced_subgraph(neg_net, V(neg_net)$boots==TRUE)#, delete.vertices = F)
  
  df_as[df_as$Method == "Biological" & df_as$Replicate == i & df_as$Association == "Negative", "Assortativity"]<-
    assortativity(sub_net_neg, V(sub_net_neg)$VT, directed = F)
  
  print(i)
}

ggplot(df_as, aes(y = Assortativity, x = Method, fill = Step))+
  facet_wrap(~Association)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # geom_jitter(alpha = 0.1)+
  geom_violin()+
  # geom_boxplot(outlier.colour = NA)+
  geom_hline(yintercept = biol_ass_pos, col = "red")+
  geom_hline(yintercept = biol_ass_neg, col = "blue")


strains<-list()
replicates<-10e2
asss<-data.frame(matrix(ncol = 2, nrow = replicates))
colnames(asss)<-c("Rep", "Assortativity")
asss$Rep<-1:replicates

for (i in 1:replicates){
  strains[[i]]<-sample(x = c(0, 1), replace = T, size = ntaxa(ps_D))
  V(net_ig)$random_VT<-strains[[i]]
  asss[asss$Rep == i, "Assortativity"]<-assortativity(neg_net, V(neg_net)$random_VT, directed=F)
  if (i %% 100 == 0){
    print(paste0(round(i*100/replicates, 4), " %"))
  }
}
asss
hist(asss$Assortativity)

biol_ass<-assortativity(pos_net, V(pos_net)$VT_s2, directed=F)

sum(asss$Assortativity<biol_ass)/length(asss$Assortativity) # p value = 0

## Positive net == 0
## Negative net == 0

### Missing compare to random strains from random model

assortativity(net_ig, V(net_ig)$VT_prob_random, directed=F)


# How many +ve interactions?
sum(netw %e% "weight" >0)/length(netw %e% "weight") #


##############################
#### Calculate abundances ####
##############################

df2<-data.frame(Genus = V(net_ig)$Genus, VT = V(net_ig)$VT_s2, random_VT = V(net_ig)$random_VT)
df2_pos<-data.frame(Genus = V(pos_net)$Genus, VT = V(pos_net)$VT_s2, random_VT = V(pos_net)$random_VT)
df2_neg<-data.frame(Genus = V(neg_net)$Genus, VT = V(neg_net)$VT_s2, random_VT = V(neg_net)$random_VT)

df1<-cbind(df, df2)
df1_pos<-cbind(df_pos00, df2_pos)
df1_neg<-cbind(df_neg00, df2_neg)

ggpairs(df1[,-which(colnames(df1)=="Genus")])

ggplot(df1, aes(y = Abundance, x = VT, col = Genus))+
  geom_jitter(height = 0)+
  scale_y_log10()

ggplot(df1, aes(y = Abundance, x = VT))+
  geom_boxplot()+
  scale_y_log10()

ggplot(df1, aes(y = Abundance, x = hubs, col = VT))+
  geom_point()+
  geom_smooth(method = "loess")+
  scale_y_log10()+scale_x_log10()

ggplot(df1, aes(y = hubs, x = VT, fill = VT))+
  geom_violin()+
  scale_y_log10()

lm_degree<-lm(data = df1, log(Degree) ~ log(Abundance) * VT )
summary(aov(lm_degree))
summary(lm_degree)

lm_degree_pos<-lm(data = df1_pos, log(Degree) ~ log(Abundance) * VT )
summary(aov(lm_degree_pos))
summary(lm_degree_pos)

lm_degree_neg<-lm(data = df1_neg, (Degree) ~ (Abundance) * VT )
summary(aov(lm_degree_neg))
summary(lm_degree_neg)

###########################################
#### Plotting nets +ve / -ve edges ########
###########################################
V(pos_net)$"Genus"

p <- ggnet2(net_ig, node.color = "VT", node.alpha = 0.4, size.min = 6,
            label = F, node.size = "degree", edge.size = 0.3,#((log(abs(netw %e% "weight"))+10)/3),
            label.size = 2, edge.color = "color") + guides(color=guide_legend(title="Vertically transmitted"), size = FALSE) +
  scale_color_manual(values = c25)
set.seed(122)
p2 <- ggnet2(pos_net, node.color = "VT_s2", node.alpha = 0.6, size.min = 6,
             label = F, node.size = "degree", edge.size = 0.3,#((log(abs(netw %e% "weight"))+10)/3),
             label.size = 4, edge.color = "color") + guides(color=guide_legend(title="Vertically transmitted"), size = FALSE) +
  scale_color_manual(values = c("darkgrey", vt_cols[2]))
theme(legend.position = "none")
p2

set.seed(122)
p2b <- ggnet2(pos_net, node.color = "Genus", node.alpha = 0.6, size.min = 6,
              label = F, node.size = "degree", edge.size = 0.3,#((log(abs(netw %e% "weight"))+10)/3),
              label.size = 4, edge.color = "color") + guides(color=guide_legend(title="Genus"), size = FALSE) +
  scale_color_manual(values = c25)
theme(legend.position = "none")
p2b

library(gridExtra)
grid.arrange(p2, p2b, nrow = 1)

p2
p2b




