library(tidyverse)
library(grid)
library(gridExtra)
library(ggpattern)
library(nlme)
library(ggrepel)

theme_set(theme_classic()+theme(axis.text.x = element_text(size=13),
                                axis.text.y = element_blank(),
                                axis.title = element_blank(),
                                title = element_text(size = 15),
                                legend.text = element_text(size = 13),
                                legend.title = element_text(size = 14)))

# Mother step 1
df_b_s1<-readRDS("~/Desktop/Bin_stats/Both_together/df_Blue_s1_25_notip.rds")
df_r_s1<-readRDS("~/Desktop/Bin_stats/Both_together/df_Red_s1_25_notip.rds")
df_g_s1<-readRDS("~/Desktop/Bin_stats/Both_together/df_Green_s1_25_notip.rds")
df_y_s1<-readRDS("~/Desktop/Bin_stats/Both_together/df_Yellow_s1_25_notip.rds")

# Mother step 2
df_b_s2<-readRDS("~/Desktop/Bin_stats/Both_together/df_Blue_Step2_25_notip.rds")
df_r_s2<-readRDS("~/Desktop/Bin_stats/Both_together/df_Red_step2_25_notip.rds")
df_g_s2<-readRDS("~/Desktop/Bin_stats/Both_together/df_Green_Step2_25_notip.rds")
df_y_s2<-readRDS("~/Desktop/Bin_stats/Both_together/df_yellow_s2_25_notip.rds")

# Father step 1
df_f_s1_f<-readRDS("~/Desktop/Bin_stats/Both_together/df_Father_s1_25_notip.rds")

# Father step 2
df_b_s2f<-readRDS("~/Desktop/Bin_stats/Both_together/df_Blue_F_s2_25_notip.rds")
df_r_s2f<-readRDS("~/Desktop/Bin_stats/Both_together/df_Red_F_s2_25_notip.rds")
df_g_s2f<-readRDS("~/Desktop/Bin_stats/Both_together/df_Green_F_s2_25_notip.rds")
df_y_s2f<-readRDS("~/Desktop/Bin_stats/Both_together/df_yellow_F_s2_25_notip.rds")

DF_b_s1<-dplyr::select(df_b_s1, c(Blue_transmission, Mean_worker, Mean_overall, Mean_Q, Taxon))
DF_r_s1<-dplyr::select(df_r_s1, c(Red_transmission, Mean_worker, Mean_overall, Mean_Q, Taxon))
DF_g_s1<-dplyr::select(df_g_s1, c(Green_transmission, Mean_worker, Mean_overall, Mean_Q, Taxon))
DF_y_s1<-dplyr::select(df_y_s1, c(Yellow_transmission, Mean_worker, Mean_overall, Mean_Q, Taxon))
colnames(DF_b_s1)<-c("Transmission", "Mean", "Mean_overall", "Mean_Q", "Taxon")
colnames(DF_r_s1)<-c("Transmission", "Mean", "Mean_overall", "Mean_Q", "Taxon")
colnames(DF_g_s1)<-c("Transmission", "Mean", "Mean_overall", "Mean_Q", "Taxon")
colnames(DF_y_s1)<-c("Transmission", "Mean", "Mean_overall", "Mean_Q", "Taxon")
DF_b_s1$color<-"deepskyblue"
DF_r_s1$color<-"brown1"
DF_g_s1$color<-"green3"
DF_y_s1$color<-"yellow3"

DF_f_s1<-dplyr::select(df_f_s1_f, c(F_transmission, Mean_worker, Taxon, Family, Genus))
colnames(DF_f_s1)<-c("Transmission", "Mean","Taxon", "Family", "Genus")
DF_f_s1$color<-"black"

DF_b_s2<-dplyr::select(df_b_s2, c(Blue_transmission, Mean_alate, Taxon, Family, Genus))
DF_r_s2<-dplyr::select(df_r_s2, c(Red_transmission, Mean_alate, Taxon, Family, Genus))
DF_g_s2<-dplyr::select(df_g_s2, c(Green_transmission, Mean_worker, Taxon, Family, Genus))
DF_y_s2<-dplyr::select(df_y_s2, c(Yellow_transmission, Mean_alate, Taxon, Family, Genus))
colnames(DF_b_s2)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
colnames(DF_r_s2)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
colnames(DF_g_s2)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
colnames(DF_y_s2)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
DF_b_s2$color<-"deepskyblue"
DF_r_s2$color<-"brown1"
DF_g_s2$color<-"green3"
DF_y_s2$color<-"yellow3"

DF_b_s2f<-dplyr::select(df_b_s2f, c(Blue_F_transmission, Mean_alate, Taxon, Family, Genus))
DF_r_s2f<-dplyr::select(df_r_s2f, c(Red_F_transmission, Mean_alate, Taxon, Family, Genus))
DF_g_s2f<-dplyr::select(df_g_s2f, c(Green_F_transmission, Mean_alate, Taxon, Family, Genus))
DF_y_s2f<-dplyr::select(df_y_s2f, c(Yellow_F_transmission, Mean_alate, Taxon, Family, Genus))
colnames(DF_b_s2f)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
colnames(DF_r_s2f)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
colnames(DF_g_s2f)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
colnames(DF_y_s2f)<-c("Transmission", "Mean", "Taxon", "Family", "Genus")
DF_b_s2f$color<-"deepskyblue"
DF_r_s2f$color<-"brown1"
DF_g_s2f$color<-"green3"
DF_y_s2f$color<-"yellow3"


all<-list(DF_b_s1, DF_r_s1, DF_y_s1, DF_g_s1,
          DF_b_s2, DF_r_s2, DF_y_s2, DF_g_s2,
          # DF_b_s3, DF_r_s3, DF_y_s3, DF_g_s3,
          DF_b_s2f, DF_r_s2f, DF_y_s2f, DF_g_s2f,
          # DF_b_s3f, DF_r_s3f, DF_y_s3f, DF_g_s3f,
          DF_f_s1
)

names(all)<-c("b_s1", "r_s1", "y_s1", "g_s1",
              "b_s2", "r_s2", "y_s2", "g_s2",
              # "b_s3", "r_s3", "y_s3", "g_s3",
              "b_s2f", "r_s2f", "y_s2f", "g_s2f",
              # "b_s3f", "r_s3f", "y_s3f", "g_s3f",
              "f_s1"
)

###### log ratio transformation & test ####

df_genus_norm<-list()

for (i in 1:length(all)){
  DF<-filter(all[[i]], Mean>0)
  DF %>%
    group_by(Taxon) %>%
    dplyr::summarise(mean_tax=mean(Mean),
                     sd=sd(Mean)) -> Tx_means
  df_genus_norm[[i]]<-data.frame(matrix(ncol = 2, nrow = length(unique(DF$Taxon))))
  colnames(df_genus_norm[[i]])<-c("Taxon", "Ab_norm")
  df_genus_norm[[i]][,"Taxon"]<-unique(DF$Taxon)
  for (y in unique(DF$Taxon)){
    mn<-Tx_means[Tx_means$Taxon==y,"mean_tax"]
    DF[DF$Taxon==y, "Mean_clr"]<-((DF[DF$Taxon==y,"Mean"]-c(mn$mean_tax)))+1e-6#/sd_sd$sd) # Problem, not calculable for taxa with a single strain.
    df_genus_norm[[i]][df_genus_norm[[i]]$Taxon == y, "Ab_norm"]<-1e-6-c(mn$mean_tax)
  }
  for (y in unique(DF$Taxon)){
    tryCatch({
      mn<-Tx_means[Tx_means$Taxon==y,"mean_tax"]
      DF[DF$Taxon==y, "Mean_Q_clr"]<-((DF[DF$Taxon==y,"Mean_Q"]-c(mn$mean_tax)))+1e-6#/sd_sd$sd) # Problem, not calculable for taxa with a single strain.
    }, error = function(e){print(e)})
  }
  m1 <- lme(Transmission~Mean_clr,random=~1|Taxon,data=DF)
  a<-anova(m1) ### Try summary as an alternative
  DF$p_val_lme<-a[2, "p-value"]
  DF$f_val_lme<-a[2, "F-value"]
  DF$df_lme<-a[2, "denDF"]
  DF$coef_lme<-c(m1$coefficients$fixed[2])
  
  k<-kruskal.test(DF$Transmission~DF$Mean_clr)
  DF$p_val_kw<-c(k$p.value)
  DF$chi_sq_kw<-c(k$statistic)
  DF$df_kw<-c(k$parameter)
  
  all[[i]]<-DF
}

summary(df_genus_norm)
# saveRDS(df_genus_norm, "~/Desktop/Bin_stats/df_genus_norm.rds")

# Do pvals correlate??
m0<-data.frame(matrix(ncol=1,nrow=length(all)))
colnames(m0)<-c("kw_p")
for (i in 1:length(all)){
  m0[i, "kw_p"]<-unique(all[[i]][, "p_val_kw"])
  m0[i, "kw_chi_sq"]<-unique(all[[i]][, "chi_sq_kw"])
  m0[i, "kw_df"]<-unique(all[[i]][, "df_kw"])
  m0[i, "lme_p"]<-unique(all[[i]][, "p_val_lme"])
  m0[i, "lme_f"]<-unique(all[[i]][, "f_val_lme"])
  m0[i, "lme_df"]<-unique(all[[i]][, "df_lme"])
  m0[i, "lin"]<-unique(all[[i]][, "color"])
  m0[i, "lme_coef"]<-unique(all[[i]][, "coef_lme"])
}
m0$group<-names(all)
m0$kw_adj<-p.adjust(m0$kw_p)
m0$lme_adj<-p.adjust(m0$lme_p)
m0$kw_sig<-m0$kw_adj<0.05
m0$lme_sig<-m0$lme_adj<0.05

ggplot(m0, aes(y=lme_p, x=lme_coef, col=lme_sig, label=group))+
  geom_point(size= 2.5)+scale_y_log10()+
  facet_grid(~lin)+
  geom_text_repel(size = 4)+ theme_bw()

m0
# write.csv(m0, "~/Desktop/Bin_stats/VTT_results/c_25_kw_lme.csv")


##################################
##### Stats overall ##############
##################################
summary(all[[1]])
names(all)

all[[1]]$Lineage<-"Blue"
all[[2]]$Lineage<-"Red"
all[[3]]$Lineage<-"Yellow"
all[[4]]$Lineage<-"Green"
all[[13]]$Lineage<-"Father"

s1_fem<-rbind(all[[1]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")],
              all[[2]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")],
              all[[3]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")],
              all[[4]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")])

m1 <- lme(Transmission~Mean_clr,random=list(~1|Lineage, ~1|Taxon), data=s1_fem)
a<-anova(m1)
a

s1_fem
colnames(s1_fem)
s1_fem %>%
  group_by(Transmission, Lineage) %>%
  summarise(meen = mean(Mean_clr),
            siid = sd(Mean_clr))

ggplot(s1_fem, aes(x=Mean_clr, fill = Transmission, alpha = 0.4))+
  geom_density()+
  scale_x_log10()

##
all[[5]]$Lineage<-"Blue"
all[[6]]$Lineage<-"Red"
all[[7]]$Lineage<-"Yellow"
all[[8]]$Lineage<-"Green"

s2_fem<-rbind(all[[5]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")],
              all[[6]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")],
              all[[7]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")],
              all[[8]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")])
# all[[13]][,c("Transmission", "Mean_clr", "Taxon", "Lineage")])

m1 <- lme(Transmission~Mean_clr,random=list(~1|Lineage, ~1|Taxon), data=s2_fem)
summary(m1)
a<-anova(m1)
a

s2_fem %>%
  group_by(Transmission, Lineage) %>%
  summarise(meen = mean(Mean_clr),
            siid = sd(Mean_clr))

ggplot(s2_fem, aes(x=Mean_clr, fill = Transmission, alpha = 0.4))+
  geom_density()+
  scale_x_log10()





##################################
##### Plots overall ##############
##################################

sst<-subset(all[[13]], Transmission==T)
plot_DF_f_s1<-ggplot(all[[13]], aes(x=Mean_clr, fill=Transmission, alpha = 0.6))+
  geom_density(outline.type = "lower")+
  scale_x_log10(breaks=c(1e-4, 1e-3))+theme_classic()+guides(alpha=F, fill=F)+
  scale_fill_manual(values = c("black", "black"))+
  ylab(NULL)+xlab("Strain abundance")+
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        title = element_text(size = 15, vjust=-10),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        plot.caption = element_text(vjust=88, hjust=0.9, size=18),
        plot.caption.position = "plot")+
  geom_density_pattern(data=sst,
                       aes(x=Mean, fill=Transmission, pattern_fill=Transmission, alpha=0.6),
                       pattern="crosshatch",
                       pattern_fill="grey",
                       colour=NA,
                       pattern_alpha=0.4)
# labs(caption= paste0("p-value = ", round(unique(all[[13]][, "p_val"]),8),
#                      "\n",
#                ". F-value = ", round(unique(all[[13]][, "f_val"]), 2)
#                ))
plot_DF_f_s1
#### d ####
p<-list()
p_2<-list()
p_3<-list()
p_4<-list()
Normalized_Q<-list()
symlog = function(x, C){    
  #Where x is your data with positive and negative values
  #where the scaling constant C determines the resolution of the data around zero.  
  t = sign(x)*log(1+abs(x)/10^C)
  return(t)
}

for (i in 1:4){
  p[[i]]<-ggplot(all[[i]], aes(x=Mean_clr, fill=Transmission, alpha = 0.6))+ # step 1 worker mother
    geom_density(outline.type = "lower")+
    scale_x_log10(breaks=c(1e-4, 1e-3))+
    theme_classic()+guides(alpha=F, fill=F)+
    scale_fill_manual(values = c("black", unique(all[[i]][,"color"])))+
    ylab(NULL)+xlab("Strain abundance")+
    theme(axis.text.x = element_text(size=13),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14))
  
  ##
  Normalized_Q[[i]] <- symlog(all[[i]]$Mean_clr, -4)
  all[[i]] %>% mutate(Mean_Q_clr_2 = Normalized_Q[[i]]) -> all[[i]]
  all[[i]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_2))  -> M_Q
  p_2[[i]]<-ggplot(all[[i]], aes(x=Mean_Q_clr_2, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i]]$Mean_clr, -5)
  all[[i]] %>% mutate(Mean_Q_clr_2 = Normalized_Q[[i]]) -> all[[i]]
  all[[i]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_2))  -> M_Q
  p_3[[i]]<-ggplot(all[[i]], aes(x=Mean_Q_clr_2, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i]]$Mean_clr, -6)
  all[[i]] %>% mutate(Mean_Q_clr_2 = Normalized_Q[[i]]) -> all[[i]]
  all[[i]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_2))  -> M_Q
  p_4[[i]]<-ggplot(all[[i]], aes(x=Mean_Q_clr_2, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i]]$color))+
    xlab(NULL)
  ##
  
  p[[(i+4)]]<-ggplot(all[[i]], aes(x=Mean_Q_clr, fill=Transmission, alpha = 0.6))+ # step 1 queen mother
    geom_density(outline.type = "lower")+
    scale_x_log10(breaks=c(1e-4, 1e-3))+theme_classic()+guides(alpha=F, fill=F)+
    scale_fill_manual(values = c("black", unique(all[[i]][,"color"])))+
    ylab(NULL)+xlab("Strain abundance")+
    theme(axis.text.x = element_text(size=13),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14))
  
  ##
  Normalized_Q[[i]] <- symlog(all[[i]]$Mean_Q_clr, -4)
  all[[i]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i]]
  all[[i]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_2[[i+4]]<-ggplot(all[[i]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i]]$Mean_Q_clr, -5)
  all[[i]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i]]
  all[[i]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_3[[i+4]]<-ggplot(all[[i]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i]]$Mean_Q_clr, -6)
  all[[i]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i]]
  all[[i]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_4[[i+4]]<-ggplot(all[[i]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  ##
  
  p[[(i+8)]]<-ggplot(all[[i+4]], aes(x=Mean_clr, fill=Transmission, alpha = 0.6))+ # step 2 mother
    geom_density(outline.type = "lower")+
    scale_x_log10(breaks=c(1e-4, 1e-3))+theme_classic()+guides(alpha=F, fill=F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]][,"color"])))+
    ylab(NULL)+xlab("Strain abundance")+
    theme(axis.text.x = element_text(size=13),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14))
  
  ##
  Normalized_Q[[i]] <- symlog(all[[i+4]]$Mean_clr, -4)
  all[[i+4]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i+4]]
  all[[i+4]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_2[[i+8]]<-ggplot(all[[i+4]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i+4]]$Mean_clr, -5)
  all[[i+4]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i+4]]
  all[[i+4]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_3[[i+8]]<-ggplot(all[[i+4]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i+4]]$Mean_clr, -6)
  all[[i+4]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i+4]]
  all[[i+4]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_4[[i+8]]<-ggplot(all[[i+4]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  ##
  
  p[[(i+12)]]<-ggplot(all[[i+8]], aes(x=Mean_clr, fill=Transmission, alpha = 0.6))+ # step 2 father
    geom_density(outline.type = "lower")+
    scale_x_log10(breaks=c(1e-4, 1e-3))+theme_classic()+guides(alpha=F, fill=F)+
    scale_fill_manual(values = c("black", unique(all[[i+8]][,"color"])))+
    ylab(NULL)+xlab("Strain abundance")+
    theme(axis.text.x = element_text(size=13),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          title = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14))
  
  ##
  Normalized_Q[[i]] <- symlog(all[[i+8]]$Mean_clr, -4)
  all[[i+8]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i+8]]
  all[[i+8]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_2[[i+12]]<-ggplot(all[[i+8]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i+8]]$Mean_clr, -5)
  all[[i+8]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i+8]]
  all[[i+8]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_3[[i+12]]<-ggplot(all[[i+8]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  Normalized_Q[[i]] <- symlog(all[[i+8]]$Mean_clr, -6)
  all[[i+8]] %>% mutate(Mean_Q_clr_3 = Normalized_Q[[i]]) -> all[[i+8]]
  all[[i+8]] %>% group_by(Transmission) %>% summarise(M = mean(Mean_Q_clr_3))  -> M_Q
  p_4[[i+12]]<-ggplot(all[[i+8]], aes(x=Mean_Q_clr_3, fill = Transmission, alpha = 0.6))+
    geom_density()+
    guides(fill = F, alpha = F)+
    scale_fill_manual(values = c("black", unique(all[[i+4]]$color)))+
    geom_vline(xintercept = M_Q$M[1], col="black") + geom_vline(xintercept = M_Q$M[2], col=unique(all[[i+4]]$color))+
    xlab(NULL)
  
  ##
  
}

blank <- grid.rect(gp=gpar(col="white"))

grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]],  # step 1 worker mother
             blank, blank, blank, plot_DF_f_s1,
             
             p[[5]], p[[6]], p[[7]], p[[8]], # step 1 queen mother
             blank, blank, blank, blank,
             
             p[[9]], p[[10]], p[[11]], p[[12]], # step 2 mother
             p[[13]], p[[14]], p[[15]], p[[16]], # step 2 father
             
             # p[[17]], p[[18]], p[[19]], p[[20]], # step 3 mother
             # p[[21]], p[[22]], p[,[23]], p[[24]], # step 3 father
             nrow = 3, top = "Both-batches: clr corrected")

grid.arrange(nrow=1, p[[9]], p[[10]], p[[11]], p[[12]])
grid.arrange(nrow=1, p_2[[9]], p_2[[10]], p_2[[11]], p_2[[12]])
grid.arrange(nrow=1, p_3[[9]], p_3[[10]], p_3[[11]], p_3[[12]])
grid.arrange(nrow=1, p_4[[9]], p_4[[10]], p_4[[11]], p_4[[12]])
grid.arrange(nrow=1, p[[1]], p[[2]], p[[3]], p[[4]])
grid.arrange(nrow=1, p_2[[1]], p_2[[2]], p_2[[3]], p_2[[4]])
grid.arrange(nrow=1, p_3[[1]], p_3[[2]], p_3[[3]], p_3[[4]])
grid.arrange(nrow=1, p_4[[1]], p_4[[2]], p_4[[3]], p_4[[4]])

library(cowplot)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]],  # step 1 worker mother
          blank, blank, blank, plot_DF_f_s1,
          
          p[[5]], p[[6]], p[[7]], p[[8]], # step 1 queen mother
          blank, blank, blank, blank,
          
          p[[9]], p[[10]], p[[11]], p[[12]], # step 2 mother
          p[[13]], p[[14]], p[[15]], p[[16]], # step 2 father
          
          p[[17]], p[[18]], p[[19]], p[[20]], # step 3 mother
          p[[21]], p[[22]], p[[23]], p[[24]], # step 3 father
          align='vh', nrow = 4)

new_list<-list(p[[1]], p[[2]], p[[3]], p[[4]],  # step 1 worker mother
               blank, blank, blank, plot_DF_f_s1,
               
               p[[5]], p[[6]], p[[7]], p[[8]], # step 1 queen mother
               blank, blank, blank, blank,
               
               p[[9]], p[[10]], p[[11]], p[[12]], # step 2 mother
               p[[13]], p[[14]], p[[15]], p[[16]], # step 2 father
               
               p[[17]], p[[18]], p[[19]], p[[20]], # step 3 mother
               p[[21]], p[[22]], p[[23]], p[[24]])

# pdf(paste0("~/Desktop/Bin_stats/Both_together/TESTING",".pdf", collapse=""),
#      onefile = TRUE,width = 12, height = 12)
do.call("grid.arrange", new_list)  
# dev.off()


