
############################
########## Fig 1 ###########
############################
library(reshape2)

taxa_s2<-list(tax_Y_s2b, tax_R_s2b, tax_G_s2b, tax_B_s2b,
              tax_Y_s2Fb, tax_R_s2Fb, tax_G_s2Fb, tax_B_s2Fb)
lins<-c("Yellow", "Red", "Green", "Blue",
        "Yellow", "Red", "Green", "Blue")

ps_norm<-transform_sample_counts(ps_3_clean, function(x) x/sum(x))
ps_norm_g<-tax_glom(ps_norm, "Genus")
sdf<-sample_data(ps_norm)

Ubitaxa<-subset(taxa_s2[[1]], taxa_s2[[1]] %in% taxa_s2[[2]] &
                  taxa_s2[[1]] %in% taxa_s2[[3]] &
                  taxa_s2[[1]] %in% taxa_s2[[4]] &
                  taxa_s2[[1]] %in% taxa_s2[[5]] &
                  taxa_s2[[1]] %in% taxa_s2[[6]] &
                  taxa_s2[[1]] %in% taxa_s2[[7]] &
                  taxa_s2[[1]] %in% taxa_s2[[8]]
)

get_abu<-function(x, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony==lin,]), ps_norm)
  ps<-prune_taxa(x, ps_norm_lin)
  ss<-sample_sums(ps)*100
  return(ss)
}
get_count<-function(x, lin){
  c<-c()
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony==lin,]), ps_norm) #Of that lineage
  ps<-prune_taxa(x, ps_norm_lin) # Select those taxa
  for(i in sample_names(ps)){
    a<-prune_samples(i, ps) # make ps of only those
    if(sum(taxa_sums(a))>0){ # Avoid error if zero
      b<-prune_taxa(taxa_sums(a)>0, a)
      b<-ntaxa(b) # How many of those taxa are in that sample?
    }else{b<-0} # Zero
    c<-c(c, b) # Add each sample
  }
  return(c)
}
get_count_pct<-function(x, lin){
  c<-c()
  c2<-c()
  
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony==lin,]), ps_norm) #Of that lineage
  ps<-prune_taxa(x, ps_norm_lin) # Select ps with those taxa
  
  x2<-subset(taxa_names(ps_norm_lin), !(taxa_names(ps_norm_lin) %in% x))
  ps2<-prune_taxa(x2, ps_norm_lin)
  
  for(i in sample_names(ps)){
    a<-prune_samples(i, ps) # make ps of only those
    if(sum(taxa_sums(a))>0){ # Avoid error if zero
      b<-prune_taxa(taxa_sums(a)>0, a)
      b<-ntaxa(b) # How many of those taxa are in that sample?
    }else{b<-0} # Zero
    c<-c(c, b)
    
    a<-prune_samples(i, ps2)
    if(sum(taxa_sums(a))>0){ # Avoid error if zero
      b2<-prune_taxa(taxa_sums(a)>0, a)
      b2<-ntaxa(b2) # How many of those taxa are in that sample?
    }else{b2<-0}
    c2<-c(c2, b2)
  }
  pct<-(c/(c+c2)*100)
  return(pct)
}

e<-data.frame()
for (i in 1:length(taxa_s2)){
  item<-taxa_s2[[i]]
  lina<-lins[i]
  d<-data.frame(rbind(cbind(get_abu(item, lina),
                            get_count_pct(item, lina),
                            get_count(item, lina)
  ),
  cbind(
    get_abu(item, "Father"),
    get_count_pct(item, "Father"),
    get_count(item, "Father")
  )))
  colnames(d)<-c("ASV_abundance", "ASV_pct", "ASV_count")#, "Genus_abundance")
  if (i %in% 1:4){d$Lineage<-lina}
  if (i %in% 5:8){d$Lineage<-paste0(lina, "_father")}
  d$SampleID<-rownames(d)
  e<-rbind(e, d)
}
sdf<-data.frame(sdf)
e
for (i in sample_names(ps_3_clean)){
  e[e$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  e[e$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
}

## Add male lineages 
f<-melt(e, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))
f

# f<-subset(f, !(Caste %in% c("T2 colony Queen", "Mother colony Queen", "Mother colony King", "T1 colony Larva")))
f<-subset(f, (Lineage %in% lins[5:8]) |
            !(Caste %in% c( "Mother colony Worker", "Female Alate")))
f<-subset(f, !(Lineage %in% lins[1:4]) |
            !(Caste %in% c("Father colony Worker", "Male Alate")))

f$Lineage<-factor(f$Lineage, levels = c("Yellow", "Yellow_father", "Red", "Red_father", "Green", "Green_father", "Blue", "Blue_father"))
f$SampleID<-factor(f$SampleID, levels = unique(f[order(f$Lineage, decreasing = T), "SampleID"]))
f$variable<-factor(f$variable, levels = c("ASV_abundance", "ASV_pct", "ASV_count"))


### Average male bars
f$Lineage<-as.character(f$Lineage)
for (i in unique(f[f$Caste %in% c("Male Alate", "Father colony Worker"), "SampleID"])){
  print(i)
  Caste<-unique(f[f$SampleID==i, "Caste"])[1]
  f[f$SampleID == i & f$variable=="ASV_abundance", "value"]<-mean(f[f$SampleID==i & f$variable=="ASV_abundance", "value"])
  f[f$SampleID == i & f$variable=="ASV_pct", "value"]<-mean(f[f$SampleID==i & f$variable=="ASV_pct", "value"])
  f[f$SampleID == i & f$variable=="ASV_count", "value"]<-mean(f[f$SampleID==i & f$variable=="ASV_count", "value"])
  f[f$SampleID==i, "Lineage"]<-"Father"
}

f$Lineage<-factor(f$Lineage, levels = c("Yellow", "Yellow_father", "Red", "Red_father", "Green", "Green_father", "Blue", "Blue_father", "Father"))
f<-subset(f, !(Caste %in% c("Mother colony Queen", "T1 colony Queen", "Mother colony King", "T2 colony Queen")) |
            !(str_detect(f$Lineage, "_father")))


### Daughter barplots
f[f$Caste == "Male Alate" , "value"]<- #& f$variable == "Abundance"
  f[f$Caste == "Male Alate", "value"]/4

f[f$Caste == "Father colony Worker" , "value"]<-
  f[f$Caste == "Father colony Worker" , "value"]/6


ggplot(f)+
  geom_col(aes(y=value, fill = Lineage, x= SampleID), position = "stack")+
  geom_point(data = subset(f, f$variable == "Count"),aes(y=value, color = Lineage, x= SampleID))+
  coord_flip()+xlab(NULL)+
  facet_grid(Pedigree+Caste~variable, scale = "free", space = "free_y")+
  scale_fill_manual(values = c("yellow3", "dark orange",
                               "red", "dark red",
                               "green", "dark green",
                               "blue", "dark blue", "darkgrey"))+
  scale_color_manual(values = c("yellow3", "dark orange",
                                "red", "dark red",
                                "green", "dark green",
                                "blue", "dark blue", "darkgrey"))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 0)
  )


####### Create overlap in daughter colonies
get_overlap<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite == "T1 colony Worker", ]), ps_norm)
  x3<-subset(x, x %in% x2 )
  x1<-subset(x, !(x %in% x2))
  x2<-subset(x2, !(x2 %in% x))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_overlap<-prune_taxa(x3, ps_norm_lin)
  ss_m<-sample_sums(ps_m)*100
  ss_f<-sample_sums(ps_f)*100
  ss_overlap<-sample_sums(ps_overlap)*100
  return(t(rbind(ss_m, ss_f, ss_overlap)))
}
get_overlap_ubi<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite == "T1 colony Worker", ]), ps_norm)
  x3<-subset(x, x %in% x2 & !(x %in% Ubitaxa))
  x1<-subset(x, !(x %in% x2) & !(x %in% Ubitaxa))
  x2<-subset(x2, !(x2 %in% x) & !(x2 %in% Ubitaxa))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_overlap<-prune_taxa(x3, ps_norm_lin)
  ps_ubi<-prune_taxa(Ubitaxa, ps_norm_lin)
  ss_m<-sample_sums(ps_m)*100
  ss_f<-sample_sums(ps_f)*100
  ss_overlap<-sample_sums(ps_overlap)*100
  ss_ubi<-sample_sums(ps_ubi)*100
  return(t(rbind(ss_m, ss_f, ss_overlap, ss_ubi)))
}
get_count_pct_overlap<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite == "T1 colony Worker", ]), ps_norm) 
  x3<-subset(x, x %in% x2)
  x1<-subset(x, !(x %in% x2))
  x2<-subset(x2, !(x2 %in% x))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_o<-prune_taxa(x3, ps_norm_lin)
  c_m<-c()
  c_f<-c()
  c_o<-c()
  c_tot<-c()
  namez<-c()
  for(i in sample_names(ps_norm_lin)){
    a_m<-prune_samples(i, ps_m) # make ps of only those
    if(sum(taxa_sums(a_m))>0){ # Avoid error if zero
      b_m<-prune_taxa(taxa_sums(a_m)>0, a_m)
      b_m<-ntaxa(b_m) # How many of those taxa are in that sample?
    }else{b_m<-0} # Zero
    c_m<-c(c_m, b_m)  
    
    a_f<-prune_samples(i, ps_f) # make ps of only those
    if(sum(taxa_sums(a_f))>0){ # Avoid error if zero
      b_f<-prune_taxa(taxa_sums(a_f)>0, a_f)
      b_f<-ntaxa(b_f) # How many of those taxa are in that sample?
    }else{b_f<-0} # Zero
    c_f<-c(c_f, b_f)  
    
    a_o<-prune_samples(i, ps_o) # make ps of only those
    if(sum(taxa_sums(a_o))>0){ # Avoid error if zero
      b_o<-prune_taxa(taxa_sums(a_o)>0, a_o)
      b_o<-ntaxa(b_o) # How many of those taxa are in that sample?
    }else{b_o<-0} # Zero
    c_o<-c(c_o, b_o)  
    
    ps_norm_lin_tax<-prune_taxa(taxa_sums(ps_norm_lin)>0, ps_norm_lin)
    a_tot<-prune_samples(i, ps_norm_lin_tax)
    b_tot<-ntaxa(a_tot)
    c_tot<-c(c_tot, b_tot)
    
    namez<-c(namez, sample_names(a_tot))
  }
  
  c_m_pct<-(c_m/c_tot)*100
  c_f_pct<-(c_f/c_tot)*100
  c_o_pct<-(c_o/c_tot)*100
  dsi<-t(rbind(c_m_pct, c_f_pct, c_o_pct))
  rownames(dsi)<-namez
  return(dsi)
}
get_count_pct_overlap_ubi<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite == "T1 colony Worker", ]), ps_norm) 
  x3<-subset(x, x %in% x2 & !(x %in% Ubitaxa))
  x1<-subset(x, !(x %in% x2) & !(x %in% Ubitaxa))
  x2<-subset(x2, !(x2 %in% x) & !(x2 %in% Ubitaxa))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_o<-prune_taxa(x3, ps_norm_lin)
  ps_ubi<-prune_taxa(Ubitaxa, ps_norm_lin)
  c_m<-c()
  c_f<-c()
  c_o<-c()
  c_tot<-c()
  c_ubi<-c()
  namez<-c()
  for(i in sample_names(ps_norm_lin)){
    a_m<-prune_samples(i, ps_m) # make ps of only those
    if(sum(taxa_sums(a_m))>0){ # Avoid error if zero
      b_m<-prune_taxa(taxa_sums(a_m)>0, a_m)
      b_m<-ntaxa(b_m) # How many of those taxa are in that sample?
    }else{b_m<-0} # Zero
    c_m<-c(c_m, b_m)  
    
    a_f<-prune_samples(i, ps_f) # make ps of only those
    if(sum(taxa_sums(a_f))>0){ # Avoid error if zero
      b_f<-prune_taxa(taxa_sums(a_f)>0, a_f)
      b_f<-ntaxa(b_f) # How many of those taxa are in that sample?
    }else{b_f<-0} # Zero
    c_f<-c(c_f, b_f)  
    
    a_o<-prune_samples(i, ps_o) # make ps of only those
    if(sum(taxa_sums(a_o))>0){ # Avoid error if zero
      b_o<-prune_taxa(taxa_sums(a_o)>0, a_o)
      b_o<-ntaxa(b_o) # How many of those taxa are in that sample?
    }else{b_o<-0} # Zero
    c_o<-c(c_o, b_o)  
    
    ps_norm_lin_tax<-prune_taxa(taxa_sums(ps_norm_lin)>0, ps_norm_lin)
    a_tot<-prune_samples(i, ps_norm_lin_tax)
    b_tot<-ntaxa(a_tot)
    c_tot<-c(c_tot, b_tot)
    
    a_ubi<-prune_samples(i, ps_ubi)
    b_ubi<-ntaxa(a_ubi)
    c_ubi<-c(c_ubi, b_ubi)
    
    namez<-c(namez, sample_names(a_tot))
  }
  
  c_m_pct<-(c_m/c_tot)*100
  c_f_pct<-(c_f/c_tot)*100
  c_o_pct<-(c_o/c_tot)*100
  c_ubi_pct<-(c_ubi/c_tot)*100
  dsi<-t(rbind(c_m_pct, c_f_pct, c_o_pct, c_ubi_pct))
  rownames(dsi)<-namez
  return(dsi)
}
get_count_overlap<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite == "T1 colony Worker", ]), ps_norm) 
  x3<-subset(x, x %in% x2)
  x1<-subset(x, !(x %in% x2))
  x2<-subset(x2, !(x2 %in% x))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_o<-prune_taxa(x3, ps_norm_lin)
  c_m<-c()
  c_f<-c()
  c_o<-c()
  c_tot<-c()
  namez<-c()
  for(i in sample_names(ps_norm_lin)){
    a_m<-prune_samples(i, ps_m) # make ps of only those
    if(sum(taxa_sums(a_m))>0){ # Avoid error if zero
      b_m<-prune_taxa(taxa_sums(a_m)>0, a_m)
      b_m<-ntaxa(b_m) # How many of those taxa are in that sample?
    }else{b_m<-0} # Zero
    c_m<-c(c_m, b_m)  
    
    a_f<-prune_samples(i, ps_f) # make ps of only those
    if(sum(taxa_sums(a_f))>0){ # Avoid error if zero
      b_f<-prune_taxa(taxa_sums(a_f)>0, a_f)
      b_f<-ntaxa(b_f) # How many of those taxa are in that sample?
    }else{b_f<-0} # Zero
    c_f<-c(c_f, b_f)  
    
    a_o<-prune_samples(i, ps_o) # make ps of only those
    if(sum(taxa_sums(a_o))>0){ # Avoid error if zero
      b_o<-prune_taxa(taxa_sums(a_o)>0, a_o)
      b_o<-ntaxa(b_o) # How many of those taxa are in that sample?
    }else{b_o<-0} # Zero
    c_o<-c(c_o, b_o)  
    
    ps_norm_lin_tax<-prune_taxa(taxa_sums(ps_norm_lin)>0, ps_norm_lin)
    a_tot<-prune_samples(i, ps_norm_lin_tax)
    b_tot<-ntaxa(a_tot)
    c_tot<-c(c_tot, b_tot)
    
    namez<-c(namez, sample_names(a_tot))
  }
  # c_m_pct<-(c_m/c_tot)*100
  # c_f_pct<-(c_f/c_tot)*100
  # c_o_pct<-(c_o/c_tot)*100
  dsi<-t(rbind(c_m, c_f, c_o))
  rownames(dsi)<-namez
  return(dsi)
}

get_count_overlap(taxa_s2[[1]], taxa_s2[[5]],"Yellow")

g<-data.frame()
for (i in 1:4){
  item<-taxa_s2[[i]]
  item2<-taxa_s2[[i+4]]
  lina<-lins[i]
  d<-data.frame(get_overlap(item, item2, lina))
  colnames(d)<-c("mother", "father", "Overlap")
  d$Lineage <- lina
  g<-rbind(d, g)
}
g$SampleID<-rownames(g)

g2<-data.frame()
for (i in 1:4){
  item<-taxa_s2[[i]]
  item2<-taxa_s2[[i+4]]
  lina<-lins[i]
  d<-data.frame(get_count_pct_overlap(item, item2, lina))
  colnames(d)<-c("mother", "father", "Overlap")
  d$Lineage <- lina
  g2<-rbind(d, g2)
}
g2$SampleID<-rownames(g2)

g3<-data.frame()
for (i in 1:4){
  item<-taxa_s2[[i]]
  item2<-taxa_s2[[i+4]]
  lina<-lins[i]
  d<-data.frame(get_count_overlap(item, item2, lina))
  colnames(d)<-c("mother", "father", "Overlap")
  d$Lineage <- lina
  g3<-rbind(d, g3)
}
g3$SampleID<-rownames(g3)


for (i in sample_names(ps_3_clean)){
  g[g$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  g[g$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
  g2[g2$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  g2[g2$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
  g3[g3$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  g3[g3$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
}

g2a<-melt(g, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))
g2b<-melt(g2, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))
g2c<-melt(g3, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))

g2a
g2b
g2c

g2c %>%
  group_by(variable) %>%
  summarize (miin = mean(value),
             seed = sd(value))
g2cb<-subset(g2c, variable != "Overlap")
lm11<-lm(data=g2cb, log(value)~variable+Lineage)
summary(aov(lm11))
aov(lm11)
lm11
shapiro.test(residuals(lm11))
lmtest::bptest(lm11)

f2<-subset(f, Caste == "Mother colony Worker" & variable == "ASV_pct")
f2 %>%
  group_by(Caste) %>%
  summarize (miin = mean(value),
             seed = sd(value))

f2<-subset(f, Caste %in% c("Male Alate", "Female Alate") & variable == "ASV_abundance")
f2 %>%
  group_by(Caste) %>%
  summarize (miin = mean(value),
             seed = sd(value))

f2<-subset(f, Caste %in% c("T1 colony Worker") & variable == "ASV_abundance")
f2 %>%
  group_by(Caste, Lineage) %>%
  summarize (miin = mean(value),
             seed = sd(value))

g2a$Lineage<-paste0(g2a$Lineage, "_", g2a$variable)
g2a$variable<-"ASV_abundance"
g2b$Lineage<-paste0(g2b$Lineage, "_", g2b$variable)
g2b$variable<-"ASV_pct"
g2c$Lineage<-paste0(g2c$Lineage, "_", g2c$variable)
g2c$variable<-"ASV_count"

f3<-subset(f, f$Caste!= "T1 colony Worker" )#| f$variable == "Count")
f2<-rbind(f3, g2a, g2b, g2c)

f2$Lineage<-as.character(f2$Lineage)
f2[f2$Lineage=="Yellow_mother", "Lineage"]<-"Yellow"
f2[f2$Lineage=="Red_mother", "Lineage"]<-"Red"
f2[f2$Lineage=="Green_mother", "Lineage"]<-"Green"
f2[f2$Lineage=="Blue_mother", "Lineage"]<-"Blue"

f2[f2$Lineage=="Blue_NA", "Lineage"]<-"Ubiquitous_taxa"
f2[f2$Lineage=="Yellow_NA", "Lineage"]<-"Ubiquitous_taxa"
f2[f2$Lineage=="Green_NA", "Lineage"]<-"Ubiquitous_taxa"
f2[f2$Lineage=="Red_NA", "Lineage"]<-"Ubiquitous_taxa"


f2$Lineage<-factor(f2$Lineage, levels = c("Yellow", "Yellow_father", "Yellow_Overlap",
                                          "Red", "Red_father", "Red_Overlap",
                                          "Green", "Green_father", "Green_Overlap",
                                          "Blue", "Blue_father", "Blue_Overlap",
                                          "Father"
))

f3<-subset(f2, f2$Caste != "T1 colony Queen" |
             !(f2$Lineage %in% paste0(lins, "_father")))

d2<-data.frame(variable = c("ASV abundance", "ASV count", "ASV_pct"),
               value = c(0, 0,0, 100, 3600, 100),
               SampleID = "52b")

f2b<-subset(f2, variable != "ASV_pct")
d2b<-subset(d2, variable != "ASV_pct")

f2b$variable<-as.character(f2b$variable)
f2b[f2b$variable=="ASV_count","variable"]<-"ASV count"
f2b[f2b$variable=="ASV_abundance","variable"]<-"ASV abundance"
f2b$variable<-factor(f2b$variable, levels = c("ASV count", "ASV abundance"))


# #### Maternally transmitted strains outnumber paternal strains in founding reproductives 
# 
# found<-subset(f3, variable == "ASV_count" & Caste %in% c("Male Alate", "Female Alate"))
# found[found$Caste=="Male Alate", "value"]<-found[found$Caste=="Male Alate", "value"]*4
# found

p<-ggplot(f2b)+
  geom_col(data = f2b, aes(y=value, fill = Lineage, x= SampleID), position = "stack", width = 0.85)+
  coord_flip()+xlab(NULL)+
  scale_fill_manual(values = c("yellow3", "darkorange3", "dark orange",
                               "brown1", "brown4", "brown3",
                               "green3", "dark green", "chartreuse4",
                               "deepskyblue", "dark blue", "blue",
                               "darkgrey", "black"))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 0)
  )+
  geom_blank(data = d2b, aes(y=value))+
  geom_hline(yintercept = 0)

p+facet_grid(Pedigree+Caste~variable, scale = "free", space = "free_y")+
  ylab(NULL)

colnames(f2)

f22<-f2
f2w<-subset(f22, f22$Caste=="T1 colony Worker")
f2w %>%
  group_by(SampleID, variable) %>%
  summarise(sam = sum(value)) -> f2w
as.data.frame(f2w[order(f2w$variable),]) 
as.data.frame(f2w[order(f2w$variable),]) %>%
  group_by(variable) %>%
  summarise(miin = mean (sam),
            saam = sd(sam))


f22%>%
  group_by(Caste, variable) %>%
  summarise(miin = mean(value),
            stan = sd(value)) -> f22
f22[f22$Caste == "Father colony Worker", c("miin", "stan")] <-
  f22[f22$Caste == "Father colony Worker", c("miin", "stan")]*6
f22[f22$Caste == "Male Alate", c("miin", "stan")] <-
  f22[f22$Caste == "Male Alate", c("miin", "stan")]*4

View(f22[order(f22$variable),])


######## Supplementary fig #########
get_overlap_ubi<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite %in% c("T1 colony Worker", "T1 colony Larva"), ]), ps_norm)
  x3<-subset(x, x %in% x2 & !(x %in% Ubitaxa))
  x1<-subset(x, !(x %in% x2) & !(x %in% Ubitaxa))
  x2<-subset(x2, !(x2 %in% x) & !(x2 %in% Ubitaxa))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_overlap<-prune_taxa(x3, ps_norm_lin)
  ps_ubi<-prune_taxa(Ubitaxa, ps_norm_lin)
  ss_m<-sample_sums(ps_m)*100
  ss_f<-sample_sums(ps_f)*100
  ss_overlap<-sample_sums(ps_overlap)*100
  ss_ubi<-sample_sums(ps_ubi)*100
  return(t(rbind(ss_m, ss_f, ss_overlap, ss_ubi)))
}
get_count_pct_overlap_ubi<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite %in% c("T1 colony Worker", "T1 colony Larva"), ]), ps_norm) 
  x3<-subset(x, x %in% x2 & !(x %in% Ubitaxa))
  x1<-subset(x, !(x %in% x2) & !(x %in% Ubitaxa))
  x2<-subset(x2, !(x2 %in% x) & !(x2 %in% Ubitaxa))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_o<-prune_taxa(x3, ps_norm_lin)
  ps_ubi<-prune_taxa(Ubitaxa, ps_norm_lin)
  c_m<-c()
  c_f<-c()
  c_o<-c()
  c_tot<-c()
  c_ubi<-c()
  namez<-c()
  for(i in sample_names(ps_norm_lin)){
    a_m<-prune_samples(i, ps_m) # make ps of only those
    if(sum(taxa_sums(a_m))>0){ # Avoid error if zero
      b_m<-prune_taxa(taxa_sums(a_m)>0, a_m)
      b_m<-ntaxa(b_m) # How many of those taxa are in that sample?
    }else{b_m<-0} # Zero
    c_m<-c(c_m, b_m)  
    
    a_f<-prune_samples(i, ps_f) # make ps of only those
    if(sum(taxa_sums(a_f))>0){ # Avoid error if zero
      b_f<-prune_taxa(taxa_sums(a_f)>0, a_f)
      b_f<-ntaxa(b_f) # How many of those taxa are in that sample?
    }else{b_f<-0} # Zero
    c_f<-c(c_f, b_f)  
    
    a_o<-prune_samples(i, ps_o) # make ps of only those
    if(sum(taxa_sums(a_o))>0){ # Avoid error if zero
      b_o<-prune_taxa(taxa_sums(a_o)>0, a_o)
      b_o<-ntaxa(b_o) # How many of those taxa are in that sample?
    }else{b_o<-0} # Zero
    c_o<-c(c_o, b_o)  
    
    ps_norm_lin_tax<-prune_taxa(taxa_sums(ps_norm_lin)>0, ps_norm_lin)
    a_tot<-prune_samples(i, ps_norm_lin_tax)
    b_tot<-ntaxa(a_tot)
    c_tot<-c(c_tot, b_tot)
    
    a_ubi<-prune_samples(i, ps_ubi)
    b_ubi<-ntaxa(a_ubi)
    c_ubi<-c(c_ubi, b_ubi)
    
    namez<-c(namez, sample_names(a_tot))
  }
  
  c_m_pct<-(c_m/c_tot)*100
  c_f_pct<-(c_f/c_tot)*100
  c_o_pct<-(c_o/c_tot)*100
  c_ubi_pct<-(c_ubi/c_tot)*100
  dsi<-t(rbind(c_m_pct, c_f_pct, c_o_pct, c_ubi_pct))
  rownames(dsi)<-namez
  return(dsi)
}
get_count_overlap_ubi<-function(x, x2, lin){
  ps_norm_lin<-prune_samples(rownames(sdf[sdf$Colony == lin & sdf$Caste_composite %in% c("T1 colony Worker", "T1 colony Larva"), ]), ps_norm) 
  x3<-subset(x, x %in% x2 & !(x %in% Ubitaxa))
  x1<-subset(x, !(x %in% x2) & !(x %in% Ubitaxa))
  x2<-subset(x2, !(x2 %in% x) & !(x2 %in% Ubitaxa))
  ps_m<-prune_taxa(x1, ps_norm_lin)
  ps_f<-prune_taxa(x2, ps_norm_lin)
  ps_o<-prune_taxa(x3, ps_norm_lin)
  ps_u<-prune_taxa(Ubitaxa, ps_norm_lin)
  c_m<-c()
  c_f<-c()
  c_o<-c()
  c_u<-c()
  c_tot<-c()
  namez<-c()
  for(i in sample_names(ps_norm_lin)){
    a_m<-prune_samples(i, ps_m) # make ps of only those
    if(sum(taxa_sums(a_m))>0){ # Avoid error if zero
      b_m<-prune_taxa(taxa_sums(a_m)>0, a_m)
      b_m<-ntaxa(b_m) # How many of those taxa are in that sample?
    }else{b_m<-0} # Zero
    c_m<-c(c_m, b_m)  
    
    a_f<-prune_samples(i, ps_f) # make ps of only those
    if(sum(taxa_sums(a_f))>0){ # Avoid error if zero
      b_f<-prune_taxa(taxa_sums(a_f)>0, a_f)
      b_f<-ntaxa(b_f) # How many of those taxa are in that sample?
    }else{b_f<-0} # Zero
    c_f<-c(c_f, b_f)  
    
    a_o<-prune_samples(i, ps_o) # make ps of only those
    if(sum(taxa_sums(a_o))>0){ # Avoid error if zero
      b_o<-prune_taxa(taxa_sums(a_o)>0, a_o)
      b_o<-ntaxa(b_o) # How many of those taxa are in that sample?
    }else{b_o<-0} # Zero
    c_o<-c(c_o, b_o)  
    
    a_u<-prune_samples(i, ps_u) # make ps of only those
    if(sum(taxa_sums(a_u))>0){ # Avoid error if zero
      b_u<-prune_taxa(taxa_sums(a_u)>0, a_u)
      b_u<-ntaxa(b_u) # How many of those taxa are in that sample?
    }else{b_u<-0} # Zero
    c_u<-c(c_u, b_u)  
    
    ps_norm_lin_tax<-prune_taxa(taxa_sums(ps_norm_lin)>0, ps_norm_lin)
    a_tot<-prune_samples(i, ps_norm_lin_tax)
    b_tot<-ntaxa(a_tot)
    c_tot<-c(c_tot, b_tot)
    
    namez<-c(namez, sample_names(a_tot))
  }
  # c_m_pct<-(c_m/c_tot)*100
  # c_f_pct<-(c_f/c_tot)*100
  # c_o_pct<-(c_o/c_tot)*100
  dsi<-t(rbind(c_m, c_f, c_o, c_u))
  rownames(dsi)<-namez
  return(dsi)
}

get_count_overlap_ubi(item, item2, lina)

g<-data.frame()
for (i in 1:4){
  item<-taxa_s2[[i]]
  item2<-taxa_s2[[i+4]]
  lina<-lins[i]
  d<-data.frame(get_overlap_ubi(item, item2, lina))
  colnames(d)<-c("mother", "father", "Overlap", "Ubiquitous")
  d$Lineage <- lina
  g<-rbind(d, g)
}
g$SampleID<-rownames(g)

g2<-data.frame()
for (i in 1:4){
  item<-taxa_s2[[i]]
  item2<-taxa_s2[[i+4]]
  lina<-lins[i]
  d<-data.frame(get_count_pct_overlap_ubi(item, item2, lina))
  colnames(d)<-c("mother", "father", "Overlap", "Ubiquitous")
  d$Lineage <- lina
  g2<-rbind(d, g2)
}
g2$SampleID<-rownames(g2)

g3<-data.frame()
for (i in 1:4){
  item<-taxa_s2[[i]]
  item2<-taxa_s2[[i+4]]
  lina<-lins[i]
  d<-data.frame(get_count_overlap_ubi(item, item2, lina))
  colnames(d)<-c("mother", "father", "Overlap", "Ubiquitous")
  d$Lineage <- lina
  g3<-rbind(d, g3)
}
g3$SampleID<-rownames(g3)


for (i in sample_names(ps_3_clean)){
  g[g$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  g[g$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
  g2[g2$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  g2[g2$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
  g3[g3$SampleID==i, "Caste"]<-
    sdf[sdf$SampleID==i, "Caste_composite"]
  g3[g3$SampleID==i, "Pedigree"]<-
    sdf[sdf$SampleID==i, "Timepoint"]
}

g2a<-melt(g, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))
g2b<-melt(g2, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))
g2c<-melt(g3, id.vars = c("Lineage", "SampleID", "Caste", "Pedigree"))

g2a$Lineage<-paste0(g2a$Lineage, "_", g2a$variable)
g2a$variable<-"ASV_abundance"
g2b$Lineage<-paste0(g2b$Lineage, "_", g2b$variable)
g2b$variable<-"ASV_pct"
g2c$Lineage<-paste0(g2c$Lineage, "_", g2c$variable)
g2c$variable<-"ASV_count"


f3<-subset(f, f$Caste!= "T1 colony Worker" )#| f$variable == "Count")
f3<-subset(f3, f3$Caste!= "T1 colony Larva" )#| f$variable == "Count")
f2<-rbind(f3, g2a, g2b, g2c)

f2$Lineage<-as.character(f2$Lineage)
f2[f2$Lineage=="Yellow_mother", "Lineage"]<-"Yellow"
f2[f2$Lineage=="Red_mother", "Lineage"]<-"Red"
f2[f2$Lineage=="Green_mother", "Lineage"]<-"Green"
f2[f2$Lineage=="Blue_mother", "Lineage"]<-"Blue"

f2[f2$Lineage=="Blue_Ubiquitous", "Lineage"]<-"Ubiquitous_taxa"
f2[f2$Lineage=="Yellow_Ubiquitous", "Lineage"]<-"Ubiquitous_taxa"
f2[f2$Lineage=="Green_Ubiquitous", "Lineage"]<-"Ubiquitous_taxa"
f2[f2$Lineage=="Red_Ubiquitous", "Lineage"]<-"Ubiquitous_taxa"


f2$Lineage<-factor(f2$Lineage, levels = c("Yellow", "Yellow_father", "Yellow_Overlap",
                                          "Red", "Red_father", "Red_Overlap",
                                          "Green", "Green_father", "Green_Overlap",
                                          "Blue", "Blue_father", "Blue_Overlap",
                                          "Father", "Ubiquitous_taxa"
))

f3<-subset(f2, f2$Caste != "T1 colony Queen" |
             !(f2$Lineage %in% paste0(lins, "_father")))
f3<-subset(f3, Caste != "T2 colony Queen")

d2<-data.frame(variable = c("ASV abundance", "ASV count", "ASV pct"),
               value = c(0, 0,0, 100, 3600, 100),
               SampleID = "52b")

f2b<-f3
f2b$variable<-as.character(f2b$variable)
f2b[f2b$variable=="ASV_count","variable"]<-"ASV count"
f2b[f2b$variable=="ASV_abundance","variable"]<-"ASV abundance"
f2b[f2b$variable=="ASV_pct","variable"]<-"ASV pct"

f2b$variable<-factor(f2b$variable, levels = c("ASV count", "ASV abundance", "ASV pct"))

f2b<-subset(f2b, f2b$Caste != "T2 colony Queen")


fq<-subset(f2b, Caste == "Mother colony Queen")
fq %>%
  group_by(Lineage, variable) %>%
  transmute(value = mean(value)) ->fq2

fq2$Caste<-"Mother colony Queen"
fq2$Pedigree<-"Mature colony"
fq2$SampleID<-fq2$Lineage

fq2
f2b2<-subset(f2b, Caste!="Mother colony Queen")
fq3<-as.data.frame(fq2[,colnames(f2b2)])
f2b2<-rbind(f2b2, fq3)

p<-ggplot(f2b2)+
  geom_col(data = f2b2, aes(y=value, fill = Lineage, x= SampleID), position = "stack", width = 0.85)+
  coord_flip()+xlab(NULL)+
  scale_fill_manual(values = c("yellow3", "darkorange3", "dark orange",
                               "brown1", "brown4", "brown3",
                               "green3", "dark green", "chartreuse4",
                               "deepskyblue", "dark blue", "blue",
                               "darkgrey", "black"))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 0)
  )+
  geom_blank(data = d2, aes(y=value))+
  geom_hline(yintercept = 0)

p+facet_grid(Pedigree+Caste~variable, scale = "free", space = "free_y")+
  ylab(NULL) # 10 - 7


## Save Table S1 ##

# write.csv2(f2b2, "~/Desktop/Bin_stats/Sup_data_paper/table_s1.csv")
# write.csv2(tax_table(ps_3_clean), "~/Desktop/Bin_stats/Sup_data_paper/table_s2.csv")
# write.csv2(otu_table(ps_3_clean), "~/Desktop/Bin_stats/Sup_data_paper/table_s3.csv")
# write.csv2(sample_data(ps_3_clean), "~/Desktop/Bin_stats/Sup_data_paper/table_s4.csv")



ftest<-subset(f2b2, Caste %in% c("Male Alate", "Female Alate") &
                variable == "ASV count")

ftest[ftest$Caste=="Male Alate", "value"]<-ftest[ftest$Caste=="Male Alate", "value"]*4

ftest %>% group_by(Caste) %>%
  summarise(miin = mean(value))

lmtest<-lm(data= ftest, log(value)~Caste)
summary(aov(lmtest))
shapiro.test(residuals(lmtest))
lmtest::bptest(lmtest)
summary(lmtest)




ftest<-subset(f2b2, Caste %in% c("T1 colony Worker") &
                variable == "ASV pct")
ftest$Lineage
ftest$Group<-c("M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M",
               "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F",
               "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared", "Shared",
               "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi", "Ubi"
)

ftest %>% group_by(Group) %>%
  summarise(miin = mean(value))

lmtest<-lm(data= ftest, log(value)~Caste)
summary(aov(lmtest))
shapiro.test(residuals(lmtest))
lmtest::bptest(lmtest)

