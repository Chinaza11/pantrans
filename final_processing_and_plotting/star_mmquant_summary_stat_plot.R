rm(list=ls())
ls() 

library(ggpubr)
library(ggplot2)
library(dplyr)
library(stringr)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# STAR
# =============================================================================

star = read.csv("final_processing_and_plotting/star_alignment.csv")

png(file="final_processing_and_plotting/star_alignment.png", width=14, height=8, units="in", res=500)

star_plot = ggplot(star, aes(x = Cat, y = Percentage, fill=factor(Categories, levels=c('Unmapped: other', 'Unmapped: too short', 'Unmapped: too many mismatches', 'Mapped to too many loci', 'Mapped to multiple loci', 'Uniquely mapped')))) + 
              geom_bar(stat = 'identity', position = 'stack') + 
              theme_bw() + 
              facet_grid(~ Samples) +
              theme(legend.position="top", 
                    legend.title = element_blank(), 
                    axis.title.x=element_blank(),
                    strip.text.x = element_text(size=8)) +
              scale_fill_manual(values = c("red", 
                                           "orange", 
                                           "blue", 
                                           "black", 
                                           "purple",
                                           "green"))
star_plot

dev.off()

# -----------------------------------------------------------------------------
# Unmapped categories grouped together

star_2 = read.csv("final_processing_and_plotting/star_alignment_2.csv")

png(file="final_processing_and_plotting/star_alignment_2.png", width=14, height=8, units="in", res=500)

star_plot_2 = ggplot(star_2, aes(x = Cat, y = Percentage, fill=factor(Categories, levels=c('Unmapped', 'Mapped to too many loci', 'Mapped to multiple loci', 'Uniquely mapped')))) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme_bw() + 
  facet_grid(~ Samples) +
  theme(legend.position="top", 
        legend.title = element_blank(), 
        legend.text=element_text(size=7),
        legend.key.size=unit(0.35, "cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x = element_text(size=4.75)) +
  scale_fill_manual(values = c("orange", 
                               "black", 
                               "purple",
                               "green"))
star_plot_2

dev.off()

# =============================================================================
# MMQUANT
# =============================================================================

mmquant = read.csv("final_processing_and_plotting/mmquant_percent.csv")

mmquant <- mmquant %>%
  mutate(Categories = str_replace_all(Categories, 
                                c("uniquely_mapped_reads" = "Uniquely mapped hits", 
                                  "ambiguous_hits" = "Ambiguous hits", 
                                  "non_uniquely_mapped_hits" = "Non-uniquely mapped hits",
                                  "unassigned_hits" = "Unassigned hits")))

png(file="final_processing_and_plotting/mmquant_percent.png", width=14, height=8, units="in", res=500)

mmquant_plot1 = ggplot(mmquant, aes(x = Cat, y = Percentage, fill=factor(Categories, levels=c('Unassigned hits', 'Non-uniquely mapped hits', 'Ambiguous hits', 'Uniquely mapped hits')))) + 
                  geom_bar(stat = 'identity', position = 'stack') + 
                  theme_bw() + 
                  facet_grid(~ samples) +
                  theme(legend.position="top", 
                        legend.title = element_blank(), 
                        axis.title.x=element_blank(),
                        strip.text.x = element_text(size=8)) +
                  scale_fill_manual(values = c("orange", 
                                               "purple", 
                                               "blue", 
                                               "green"))
mmquant_plot1

dev.off()


mmquant_2 = read.csv("final_processing_and_plotting/mmquant_no.csv")

mmquant_2 <- mmquant_2 %>%
  mutate(Categories = str_replace_all(Categories, 
                                      c("uniquely_mapped_reads" = "Uniquely mapped hits", 
                                        "ambiguous_hits" = "Ambiguous hits", 
                                        "non_uniquely_mapped_hits" = "Non-uniquely mapped hits",
                                        "unassigned_hits" = "Unassigned hits")))

png(file="final_processing_and_plotting/mmquant_no.png", width=14, height=8, units="in", res=500)
mmquant_plot2 = ggplot(mmquant_2, aes(x = Cat, y = Million, fill=factor(Categories, levels=c('Unassigned hits', 'Non-uniquely mapped hits', 'Ambiguous hits', 'Uniquely mapped hits')))) + 
                  geom_bar(stat = 'identity', position = 'stack') + 
                  theme_bw() + 
                  facet_grid(~ samples) +
                  theme(legend.position="top", 
                        legend.title = element_blank(), 
                        axis.title.x=element_blank(),
                        strip.text.x = element_text(size=8)) +
                  scale_fill_manual(values = c("orange", 
                                               "purple", 
                                               "blue", 
                                               "green")) +
                  scale_y_continuous(breaks = seq(0,30,by=5), limits=c(0,30))
mmquant_plot2

dev.off()

# -----------------------------------------------------------------------------
# Ambiguous hits deleted (reason: not a lot of them and to help readers focus on the core result)

mmquant_2.1 = mmquant_2 %>% filter(Categories != 'Ambiguous hits')

png(file="final_processing_and_plotting/mmquant_no_2.png", width=14, height=8, units="in", res=500)
mmquant_plot2.1 = ggplot(mmquant_2.1, aes(x = Cat, y = Million, fill=factor(Categories, levels=c('Unassigned hits', 'Non-uniquely mapped hits', 'Uniquely mapped hits')))) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme_bw() + 
  facet_grid(~ samples) +
  theme(legend.position="top", 
        legend.title = element_blank(), 
        legend.text=element_text(size=7),
        legend.key.size=unit(0.35, "cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        strip.text.x = element_text(size=4.75)) +
  scale_fill_manual(values = c("orange", 
                               "purple", 
                               "green")) +
  scale_y_continuous(breaks = seq(0,30,by=5), limits=c(0,30))
mmquant_plot2.1

dev.off()

# =============================================================================
# Final plots: version 1.0
# =============================================================================

png(file="final_processing_and_plotting/star_n_mmquant_plot.png", width=12, height=12, units="in", res=300)

ggarrange(star_plot, mmquant_plot2, 
          nrow=2,
          labels = "AUTO")

dev.off()

# =============================================================================
# Final plots: version 2.0 
# =============================================================================

png(file="final_processing_and_plotting/star_n_mmquant_plot_2.jpg", width=20, height=20, units="cm", res=300)

ggarrange(star_plot_2, mmquant_plot2.1, 
          nrow=2,
          labels = "AUTO")

dev.off()

