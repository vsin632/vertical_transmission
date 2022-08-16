#### Step 3 receiver automated ####

Lins<-c("Yellow", "Red", "Green", "Blue")
Letters<-c("Y", "R", "G", "B")
Colours<-c("yellow3", "brown1", "green3", "deepskyblue")
list_df_Q<-list()
list_df_W<-list()
list_df_L<-list()
plot_Q<-list()
plot_W<-list()
plot_L<-list()
plot_Q_clr<-list()
plot_W_clr<-list()
plot_L_clr<-list()
Normalized_Q<-list()
plot_Q_clr_2<-list()
plot_Q_clr_3<-list()
plot_Q_clr_4<-list()
ps_3_rel<-transform_sample_counts(ps_3_clean, function(x) x/sum(x))
ps_3_g<-tax_glom(ps_3_rel, "Genus")

symlog = function(x, C){    
  #Where x is your data with positive and negative values
  #where the scaling constant C determines the resolution of the data around zero.  
  t = sign(x)*log(1+abs(x)/10^C)
  return(t)
}


theme_set(theme_classic()+theme(axis.text.x = element_text(size=13, angle = 90),
                                axis.text.y = element_blank(),
                                axis.title = element_blank(),
                                title = element_text(size = 15),
                                legend.text = element_text(size = 13),
                                legend.title = element_text(size = 14)))


df<-data.frame(matrix(ncol=length(taxa_names(ps_3_clean))+1, nrow = (nsamples(ps_3_clean))))
df$X1<-sample_names(ps_3_clean)
colnames(df)<-c("Sample", taxa_names(ps_3_clean))

for (Lin in 1:4){
  Letter<-Letters[Lin]
  for (i in sample_names(get(paste0("ps_", Letter, "_M1")))){
    ps<-prune_samples(sample_names(ps_3_rel) %in% i, ps_3_rel)
    df[df$Sample==i, 2:ncol(df)]<-taxa_sums(ps)
    print(i)
  }
}
df2<-data.frame(t(df))

df2$Yellow_transmission<-rownames(df2) %in% tax_Y_s1
df2$Blue_transmission<-rownames(df2) %in% tax_B_s1
df2$Red_transmission<-rownames(df2) %in% tax_R_s1
df2$Green_transmission<-rownames(df2) %in% tax_G_s1
colnames(df2)[1:63]<-paste0("Sample_", sample_names(ps_3_clean))
df3<-df2[2:nrow(df2),]

plot_values<-seq(-1e-2, 1e-2, 1e-3)
plot_values<-c(-0.01, -0.001, -1e-4, 0, 1e-4, 0.001, 0.01)
plot_values2<-c(plot_values[1:3], "Genus mean", plot_values[5:7])
plot_ticks<-symlog(plot_values, -5)
plot_values3<-c(plot_values2[1:3], "0", plot_values2[5:7])
plot_values3<-c(-0.01, -0.001, -1e-4, -1e-5, 0, 1e-5, 1e-4, 0.001, 0.01)
plot_ticks3<-symlog(plot_values3, -5)
# df_yellow is list_df[[1]]
list_df<-list()
for(Lin in 1:4){
  Lineage<-Lins[Lin]
  Letter<-Letters[Lin]
  colour<-Colours[Lin]
  samples<-sample_names(get(paste0("ps_", Letter, "_M1")))
  castes<-sample_data(get(paste0("ps_", Letter, "_M1")))$Caste_composite
  meta_df<-data.frame(samples = samples, castes = castes)
  list_df[[Lin]]<-select(df3, c(paste0("Sample_", sample_names(get(paste0("ps_", Letter, "_M1")))), paste0(Lineage, "_transmission")))
  cnames<-colnames(list_df[[Lin]])
  colnames(list_df[[Lin]])<-c(cnames[1:2], "Transmission")
  for (i in 1:nrow(list_df[[Lin]])){ # nrow(list_df[[Lin]]) ###Get means of abundances per caste
    Queen_samples<-list_df[[Lin]][i, paste0("Sample_", meta_df[meta_df$castes=="Mother colony Worker", "samples"])]
    list_df[[Lin]][i, "Mean_Q"]<-mean(c(as.numeric(as.character(Queen_samples[,1])),
                                        as.numeric(as.character(Queen_samples[,2]))))
    
    r00<-rownames(list_df[[Lin]])[i] # Get taxonomic rank 
    list_df[[Lin]][i, "Taxon"]<-paste0(tax_table(ps_3_clean)[r00, 1:5], collapse = ":_")
    list_df[[Lin]][i, "Family"]<-paste0(tax_table(ps_3_clean)[r00, 4], collapse = "_")
    list_df[[Lin]][i, "Genus"]<-paste0(tax_table(ps_3_clean)[r00, 5], collapse = "_")
    if (i %%1000 ==0){
      print(i)
    }
  }
  
  list_df_Q[[Lin]]<-subset(list_df[[Lin]], list_df[[Lin]]$Mean_Q>0)
  
  list_df_Q[[Lin]] %>% group_by(Taxon) %>% dplyr::summarise(mean_tax = mean(Mean_Q)) -> Q_means
  
  for (y in unique(Q_means$Taxon)){
    mn<-Q_means[Q_means$Taxon==y,"mean_tax"]
    list_df_Q[[Lin]][list_df_Q[[Lin]][,"Taxon"]==y, "Mean_Q_clr"]<-
      (list_df_Q[[Lin]][list_df_Q[[Lin]][,"Taxon"]==y,"Mean_Q"]-c(mn$mean_tax))#/sd_sd$sd) # Problem, not calculable for taxa with a single strain.
  }
  
  
  plot_Q[[Lin]]<-ggplot(data=list_df_Q[[Lin]], aes(x=Mean_Q, fill=Transmission, alpha=0.6))+
    geom_density(outline.type = "lower")+scale_x_log10()+ guides(alpha=F, fill=F)+
    scale_fill_manual(values = c("black", colour))+
    ylab(NULL)+xlab("Strain abundance")
  
  plot_Q_clr[[Lin]]<-ggplot(data=list_df_Q[[Lin]], aes(x=Mean_Q_clr, fill=Transmission, alpha=0.6))+
    geom_density(outline.type = "lower")+scale_x_log10()+ guides(alpha=F, fill=F)+
    scale_fill_manual(values = c("black", colour))+
    ylab(NULL)+xlab("Strain abundance")
  
  Normalized_Q[[Lin]] <- symlog(list_df_Q[[Lin]]$Mean_Q_clr, -4)
  list_df_Q[[Lin]] %>% mutate(Mean_Q_clr_2 = Normalized_Q[[Lin]]) -> list_df_Q[[Lin]]
  list_df_Q[[Lin]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_2))  -> M_Q
  plot_Q_clr_2[[Lin]]<-ggplot(list_df_Q[[Lin]], aes(x=Mean_Q_clr_2, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", colour))+
    geom_vline(xintercept = M_Q$M[1], col="black", size = 2) + geom_vline(xintercept = M_Q$M[2], col=colour, size = 2)+
    xlab(NULL)+
    scale_x_continuous(breaks = plot_ticks, labels = plot_values2)
  
  Normalized_Q[[Lin]] <- symlog(list_df_Q[[Lin]]$Mean_Q_clr, -5)
  list_df_Q[[Lin]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[Lin]]) -> list_df_Q[[Lin]]
  list_df_Q[[Lin]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  plot_Q_clr_3[[Lin]]<-ggplot(list_df_Q[[Lin]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", colour))+
    geom_vline(xintercept = M_Q$M[1], col="black", size = 2) + geom_vline(xintercept = M_Q$M[2], col=colour, size = 2)+
    xlab(NULL)+
    scale_x_continuous(breaks = plot_ticks, labels = plot_values2)
  
  Normalized_Q[[Lin]] <- symlog(list_df_Q[[Lin]]$Mean_Q_clr, -6)
  list_df_Q[[Lin]] %>% mutate(Mean_Q_clr_4 = Normalized_Q[[Lin]]) -> list_df_Q[[Lin]]
  list_df_Q[[Lin]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_4))  -> M_Q
  plot_Q_clr_4[[Lin]]<-ggplot(list_df_Q[[Lin]], aes(x=Mean_Q_clr_4, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", colour))+
    geom_vline(xintercept = M_Q$M[1], col="black", size = 2) + geom_vline(xintercept = M_Q$M[2], col=colour, size = 2)+
    xlab(NULL)+
    scale_x_continuous(breaks = plot_ticks, labels = plot_values2)
}

library(gridExtra)
do.call(grid.arrange, plot_Q)
do.call(grid.arrange, plot_Q_clr)
do.call(grid.arrange, plot_Q_clr_2)
do.call(grid.arrange, plot_Q_clr_3)
do.call(grid.arrange, plot_Q_clr_4)7

grid.arrange(plot_Q_clr_3[[1]], plot_Q_clr_3[[2]], plot_Q_clr_3[[3]], plot_Q_clr_3[[4]], nrow = 1) # 10-3

pdf("~/Desktop/R/Vertical_transmission_2020/final_figs/3_4_main/step_1_emitter.pdf",
    width = 10, height = 3, onefile = T)
grid.arrange(plot_Q_clr_3[[1]], plot_Q_clr_3[[2]], plot_Q_clr_3[[3]], plot_Q_clr_3[[4]], nrow = 1) # 10-3
dev.off()

# grid.arrange(plot_W[[1]], plot_W[[2]], plot_W[[3]], plot_W[[4]], 
#              plot_Q[[1]], plot_Q[[2]], plot_Q[[3]], plot_Q[[4]], 
#              plot_L[[1]], plot_L[[2]], plot_L[[3]], plot_L[[4]], nrow = 3)
# 
# grid.arrange(plot_Q_clr[[1]], plot_Q_clr[[2]], plot_Q_clr[[3]], plot_Q_clr[[4]], 
#              plot_L_clr[[1]], plot_L_clr[[2]], plot_L_clr[[3]], plot_L_clr[[4]], 
#              plot_W_clr[[1]], plot_W_clr[[2]], plot_W_clr[[3]], plot_W_clr[[4]], nrow = 3)
# 
# grid.arrange(plot_Q_clr_2[[1]], plot_Q_clr_2[[2]], plot_Q_clr_2[[3]], plot_Q_clr_2[[4]], 
#              plot_L_clr_2[[1]], plot_L_clr_2[[2]], plot_L_clr_2[[3]], plot_L_clr_2[[4]], 
#              plot_W_clr_2[[1]], plot_W_clr_2[[2]], plot_W_clr_2[[3]], plot_W_clr_2[[4]], nrow = 3)
# 
# grid.arrange(plot_Q_clr_3[[1]], plot_Q_clr_3[[2]], plot_Q_clr_3[[3]], plot_Q_clr_3[[4]], 
#              plot_L_clr_3[[1]], plot_L_clr_3[[2]], plot_L_clr_3[[3]], plot_L_clr_3[[4]], 
#              plot_W_clr_3[[1]], plot_W_clr_3[[2]], plot_W_clr_3[[3]], plot_W_clr_3[[4]], nrow = 3)
# 
# grid.arrange(plot_Q_clr_4[[1]], plot_Q_clr_4[[2]], plot_Q_clr_4[[3]], plot_Q_clr_4[[4]], 
#              plot_L_clr_4[[1]], plot_L_clr_4[[2]], plot_L_clr_4[[3]], plot_L_clr_4[[4]], 
#              plot_W_clr_4[[1]], plot_W_clr_4[[2]], plot_W_clr_4[[3]], plot_W_clr_4[[4]], nrow = 3)


# p_W_Y <- list()
p_Q_Y <- list()
# p_L_Y <- list()
# p_W_R <- list()
p_Q_R <- list()
# p_L_R <- list()
# p_W_G <- list()
p_Q_G <- list()
# p_L_G <- list()
# p_W_B <- list()
p_Q_B <- list()
# p_L_B <- list()

# ps_W<-list()
ps_Q<-list()
# ps_L<-list()
# length(ps_W)<-4
length(ps_Q)<-4
# length(ps_L)<-4
# length(ps_Q)

for(Lin in 1:4){
  Lineage<-Lins[Lin]
  Letter<-Letters[Lin]
  colour<-Colours[Lin]
  samples<-sample_names(get(paste0("ps_", Letter, "_M1")))
  castes<-sample_data(get(paste0("ps_", Letter, "_M1")))$Caste_composite
  meta_df<-data.frame(samples = samples, castes = castes)
  
  cnames<-colnames(list_df[[Lin]])
  colnames(list_df[[Lin]])<-c(cnames[1:2], "Transmission", cnames[4:length(cnames)])
  
  # for (CAST in unique(castes)){
  ps_3_g_sub<-subset_samples(ps_3_g, sample_data(ps_3_g)[,"Colony"]==Lineage &
                               sample_data(ps_3_g)[,"Caste_composite"]=="Mother colony Worker")
  ps_10<-prune_taxa(names(sort(decreasing = T, taxa_sums(ps_3_g_sub)))[1:20], ps_3_g)
  # Caste<-str_split(CAST, "")[[1]][11]
  for (i in 1:20){
    df_ali<-subset(list_df[[Lin]], Taxon==paste0(tax_table(ps_10)[i,1:5], collapse=":_"))
    colnames(df_ali)<-c(colnames(df_ali)[1:2], "Transmission", colnames(df_ali)[4:7])
    
    Normalized_Q[[Lin]] <- symlog(df_ali$Mean_Q, -5)
    df_ali %>% mutate(Mean_Q_clr_2 = Normalized_Q[[Lin]]) -> df_ali
    df_ali %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_2))  -> M_Q
    
    ps_Q[[Lin]][[i]]<-ggplot(data=df_ali, aes(x=Mean_Q_clr_2, fill=Transmission, alpha=0.6))+
      geom_density(outline.type = "lower")+guides(alpha=F, fill=F)+
      scale_fill_manual(values = c("black", colour))+
      ggtitle(tax_table(ps_10)[i,5],
              tax_table(ps_10)[i,4])+ylab(NULL)+xlab(NULL)+
      scale_x_continuous(breaks = plot_ticks3, labels = plot_values3, limits = c(0, NA))+
      geom_vline(xintercept = M_Q$M[1], col="black", size =2) + geom_vline(xintercept = M_Q$M[2], col=colour, size = 2)
    
  }
}

pdf("~/Desktop/R/Vertical_transmission_2020/final_figs/3_4_supplementary/step_1_emitter_yellow.pdf",
    width = 14, height = 12, onefile = T)
do.call(grid.arrange, ps_Q[[1]]) # 14 - 12
dev.off()
pdf("~/Desktop/R/Vertical_transmission_2020/final_figs/3_4_supplementary/step_1_emitter_red.pdf",
    width = 14, height = 12, onefile = T)
do.call(grid.arrange, ps_Q[[2]])
dev.off()
pdf("~/Desktop/R/Vertical_transmission_2020/final_figs/3_4_supplementary/step_1_emitter_green.pdf",
    width = 14, height = 12, onefile = T)
do.call(grid.arrange, ps_Q[[3]])
dev.off()
pdf("~/Desktop/R/Vertical_transmission_2020/final_figs/3_4_supplementary/step_1_emitter_blue.pdf",
    width = 14, height = 12, onefile = T)
do.call(grid.arrange, ps_Q[[4]])
dev.off()


do.call(grid.arrange, ps_Q[[1]]) # 14 - 12
do.call(grid.arrange, ps_Q[[2]])
do.call(grid.arrange, ps_Q[[3]])
do.call(grid.arrange, ps_Q[[4]])

#### Statistics
library(nlme)
m0_Q<-data.frame(matrix(ncol=1,nrow=length(list_df_Q)))

for (Lin in 1:length(list_df_Q)){
  m1 <- lme(Transmission~Mean_Q_clr,random=~1|Taxon,data=list_df_Q[[Lin]])
  a<-anova(m1) ### Try summary as an alternative
  list_df_Q[[Lin]]$p_val_lme<-a[2, "p-value"]
  list_df_Q[[Lin]]$f_val_lme<-a[2, "F-value"]
  list_df_Q[[Lin]]$df_lme<-a[2, "denDF"]
  list_df_Q[[Lin]]$coef_lme<-c(m1$coefficients$fixed[2])
  
  k<-kruskal.test(list_df_Q[[Lin]]$Transmission~list_df_Q[[Lin]]$Mean_Q_clr)
  list_df_Q[[Lin]]$p_val_kw<-c(k$p.value)
  list_df_Q[[Lin]]$chi_sq_kw<-c(k$statistic)
  list_df_Q[[Lin]]$df_kw<-c(k$parameter)
  
  m0_Q[Lin, "kw_p"]<-unique(list_df_Q[[Lin]][, "p_val_kw"])
  m0_Q[Lin, "kw_chi_sq"]<-unique(list_df_Q[[Lin]][, "chi_sq_kw"])
  m0_Q[Lin, "kw_df"]<-unique(list_df_Q[[Lin]][, "df_kw"])
  m0_Q[Lin, "lme_p"]<-unique(list_df_Q[[Lin]][, "p_val_lme"])
  m0_Q[Lin, "lme_f"]<-unique(list_df_Q[[Lin]][, "f_val_lme"])
  m0_Q[Lin, "lme_df"]<-unique(list_df_Q[[Lin]][, "df_lme"])
  m0_Q[Lin, "lin"]<-  colour<-Colours[Lin]
  m0_Q[Lin, "lme_coef"]<-unique(list_df_Q[[Lin]][, "coef_lme"])
}
# for (Lin in 1:length(list_df_W)){
#   m1 <- lme(Transmission~Mean_W_clr,random=~1|Taxon,data=list_df_W[[Lin]])
#   a<-anova(m1) ### Try summary as an alternative
#   list_df_W[[Lin]]$p_val_lme<-a[2, "p-value"]
#   list_df_W[[Lin]]$f_val_lme<-a[2, "F-value"]
#   list_df_W[[Lin]]$df_lme<-a[2, "denDF"]
#   list_df_W[[Lin]]$coef_lme<-c(m1$coefficients$fixed[2])
#   
#   k<-kruskal.test(list_df_W[[Lin]]$Transmission~list_df_W[[Lin]]$Mean_W_clr)
#   list_df_W[[Lin]]$p_val_kw<-c(k$p.value)
#   list_df_W[[Lin]]$chi_sq_kw<-c(k$statistic)
#   list_df_W[[Lin]]$df_kw<-c(k$parameter)
#   
#   m0_W[Lin, "kw_p"]<-unique(list_df_W[[Lin]][, "p_val_kw"])
#   m0_W[Lin, "kw_chi_sq"]<-unique(list_df_W[[Lin]][, "chi_sq_kw"])
#   m0_W[Lin, "kw_df"]<-unique(list_df_W[[Lin]][, "df_kw"])
#   m0_W[Lin, "lme_p"]<-unique(list_df_W[[Lin]][, "p_val_lme"])
#   m0_W[Lin, "lme_f"]<-unique(list_df_W[[Lin]][, "f_val_lme"])
#   m0_W[Lin, "lme_df"]<-unique(list_df_W[[Lin]][, "df_lme"])
#   m0_W[Lin, "lin"]<-  colour<-Colours[Lin]
#   m0_W[Lin, "lme_coef"]<-unique(list_df_W[[Lin]][, "coef_lme"])
# }
# for (Lin in 1:length(list_df_L)){
#   m1 <- lme(Transmission~Mean_L_clr,random=~1|Taxon,data=list_df_L[[Lin]])
#   a<-anova(m1) ### Try summary as an alternative
#   list_df_L[[Lin]]$p_val_lme<-a[2, "p-value"]
#   list_df_L[[Lin]]$f_val_lme<-a[2, "F-value"]
#   list_df_L[[Lin]]$df_lme<-a[2, "denDF"]
#   list_df_L[[Lin]]$coef_lme<-c(m1$coefficients$fixed[2])
#   
#   k<-kruskal.test(list_df_L[[Lin]]$Transmission~list_df_L[[Lin]]$Mean_L_clr)
#   list_df_L[[Lin]]$p_val_kw<-c(k$p.value)
#   list_df_L[[Lin]]$chi_sq_kw<-c(k$statistic)
#   list_df_L[[Lin]]$df_kw<-c(k$parameter)
#   
#   m0_L[Lin, "kw_p"]<-unique(list_df_L[[Lin]][, "p_val_kw"])
#   m0_L[Lin, "kw_chi_sq"]<-unique(list_df_L[[Lin]][, "chi_sq_kw"])
#   m0_L[Lin, "kw_df"]<-unique(list_df_L[[Lin]][, "df_kw"])
#   m0_L[Lin, "lme_p"]<-unique(list_df_L[[Lin]][, "p_val_lme"])
#   m0_L[Lin, "lme_f"]<-unique(list_df_L[[Lin]][, "f_val_lme"])
#   m0_L[Lin, "lme_df"]<-unique(list_df_L[[Lin]][, "df_lme"])
#   m0_L[Lin, "lin"]<-  colour<-Colours[Lin]
#   m0_L[Lin, "lme_coef"]<-unique(list_df_L[[Lin]][, "coef_lme"])
# }


m0_Q$kw_adj<-p.adjust(m0_Q$kw_p)
m0_Q$lme_adj<-p.adjust(m0_Q$lme_p)
m0_Q$kw_sig<-m0_Q$kw_adj<0.05
m0_Q$lme_sig<-m0_Q$lme_adj<0.05

# m0_W$kw_adj<-p.adjust(m0_W$kw_p)
# m0_W$lme_adj<-p.adjust(m0_W$lme_p)
# m0_W$kw_sig<-m0_W$kw_adj<0.05
# m0_W$lme_sig<-m0_W$lme_adj<0.05
# 
# m0_L$kw_adj<-p.adjust(m0_L$kw_p)
# m0_L$lme_adj<-p.adjust(m0_L$lme_p)
# m0_L$kw_sig<-m0_L$kw_adj<0.05
# m0_L$lme_sig<-m0_L$lme_adj<0.05

library(ggrepel)
ggplot(m0_Q, aes(y=lme_adj, x=lme_coef, col=lme_sig, label = round(lme_adj, 4)))+
  geom_point(size= 2.5)+scale_y_log10()+
  facet_grid(~lin)+
  geom_text_repel(size = 4)+ theme_bw()+
  theme(axis.text.y = element_text(size=13))+
  ggtitle("Queen")

ggplot(m0_W, aes(y=lme_adj, x=lme_coef, col=lme_sig, label = round(lme_adj, 4)))+
  geom_point(size= 2.5)+scale_y_log10()+
  facet_grid(~lin)+
  geom_text_repel(size = 4)+ theme_bw()+
  theme(axis.text.y = element_text(size=13))+
  ggtitle("Worker")

ggplot(m0_L, aes(y=lme_adj, x=lme_coef, col=lme_sig, label = round(lme_adj, 4)))+
  geom_point(size= 2.5)+scale_y_log10()+
  facet_grid(~lin)+
  geom_text_repel(size = 4)+ theme_bw()+
  theme(axis.text.y = element_text(size=13))+
  ggtitle("Larva")



#### Stats together ####

list_df_Q[[1]]$Lineage <- "Yellow"
list_df_Q[[2]]$Lineage <- "Red"
list_df_Q[[3]]$Lineage <- "Green"
list_df_Q[[4]]$Lineage <- "Blue"

df_Q_new<-rbind(list_df_Q[[1]][,3:16],
                list_df_Q[[2]][,3:16],
                list_df_Q[[3]][,3:16],
                list_df_Q[[4]][,3:16])

dim(df_Q_new)
m1 <- lme(Transmission~Mean_Q_clr,random=list(~1|Taxon, ~1|Lineage), data=df_Q_new)
a<-anova(m1)
a

