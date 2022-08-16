##############################
######## Rarefaction #########
##############################

ps_3_clean
library(plyr)
library(reshape2)
library(ggplot2)
set.seed(123)

calculate_rarefaction_curves <- function(psFIN, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  estimate_rarified_richness <- function(psFIN, measures, depth) {
    if(max(sample_sums(psFIN)) < depth) return()
    psFIN <- prune_samples(sample_sums(psFIN) >= depth, psFIN)
    rarified_psFIN <- rarefy_even_depth(psFIN, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psFIN, measures = measures)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psFIN = psFIN, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}
reps<-rep(c(1, 10, 100, seq(1000, 100000, 1000)), each = 10)

rarefaction_curve_data_Obs <- calculate_rarefaction_curves(ps_3_clean, 'Observed', reps)

rarefaction_curve_data_summary_Obs <- ddply(rarefaction_curve_data_Obs, c('Depth', 'Sample', "Measure"), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_Obs$Measure<-unique(rarefaction_curve_data_Obs$Measure)

rareplot<-rbind(rarefaction_curve_data_summary_Obs, rarefaction_curve_data_summary_Obs)

new_samples<-c(paste0("X", sample_names(ps_3_clean)[1:55]), "Blue.A", "Blue.B", "Green.A", "Green.B", "Red.A", "Red.B", "Yellow.A", "Yellow.B")

sample_names(ps_3_clean)<-new_samples
sample_data(ps_3_clean)[,"SampleID"]<-sample_names(ps_3_clean)

sdf<-data.frame(sample_data(ps_3_clean))
rareplot$Caste<-NA
rareplot$Caste<-as.character(rareplot$Caste)
for (i in unique(rareplot$Sample)){
  rareplot[rareplot$Sample == i, "Caste"]<-as.character(sdf[sdf$SampleID == i, "Caste_composite"])
}

ggplot(data = subset(rareplot, Measure=="Observed" & Caste != "T2 colony Queen"), mapping = aes(
  x = Depth, y = Alpha_diversity_mean, col = Caste,
  ymin = Alpha_diversity_mean - Alpha_diversity_sd, ymax = Alpha_diversity_mean + Alpha_diversity_sd,
  group = Sample)) + geom_line(alpha=1) + ylab("Observed ASV Richness")+
  ggtitle("Observed Richness Rarefaction Plot") +
  theme_bw()+
  scale_color_manual(values = c25) # 8 - 4







