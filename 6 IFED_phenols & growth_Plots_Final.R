# Programmer: Danielle Stevens
# Last updated on: 11/17/2026
# 
# Purpose of code: This code reads in the results from IFED analyses
#                  and then generates figues for the paper 
#                  and a presentation for ISEE-NAC


#############################
## Plots for main analyses ##
#############################

library(dplyr)
library(tidyverse)
library(gtable)
library(reprex)
library(egg)
library(extrafont)
library(ggplotify)
library(stringr)
library(ggplot2)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
library(haven)
library(grid)
library(gridExtra)
library(scales)
library(ggpubr)
library(egg)
library(patchwork)


#### Exposure list
exposure.labels <- c('BPA', 'Benzophenone-3', '2,4-dichlorophenol', '2,5-dichlorophenol', 
                     'Methyl paraben', 'Propyl paraben', 'Triclosan')
exposures <- c('bpa', 'bp_3', 'dcp_24', 'dcp_25', 'm_pb', 'p_pb', 'tcs')
exposures.obrien.v68 <- paste0(exposures,"_obrien_v68")
exposures.obrien.v12 <- paste0(exposures,"_obrien_v12")
exposures.obrien <- paste0(exposures,"_obrien_gm")
exposures.obrien.ln <- paste0(exposures,"_obrien_gm_ln")
exposures.obrien.v68.ln <- paste0(exposures,"_obrien_v68_ln")
exposures.obrien.v12.ln <- paste0(exposures,"_obrien_v12_ln")




###########################
## Exposure correlations ##
###########################
library(pheatmap)

#Correlations between exposures
cormat <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Correlation Matrix.csv")
rownames(cormat) <- cormat[,1]
cormat <- cormat[,2:8]
colnames(cormat) <- rownames(cormat)

#Plot
num_labels <- matrix("", nrow=nrow(cormat), ncol=ncol(cormat))
num_labels[upper.tri((cormat))] <- cormat[upper.tri(cormat)]
cormat[upper.tri(cormat)] <- NA

p <- pheatmap(cormat,
              color=colorRampPalette(c('navy', "white", "firebrick3"))(100),
              breaks = seq(-1,1, length.out = 101),
              cluster_rows=F, 
              cluster_cols=F,
              gaps_row=c(1,2,4,6),
              gaps_col=c(1,2,4,6),
              legend=T,
              legend_breaks=-1:1,
              legend_labels=c(-1,0,1),
              cellwidth=60, 
              cellheight=60,
              display_numbers=num_labels,
              number_color = "black",
              na_col=,
              fontsize=25)


png(height=1800, width=1800, file="J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Correlation Matrix.png", type = "cairo")
p
dev.off()


##########################
## Trajectory Groupings ##
##########################

## Dataframes
ifed1 <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Adjusted Multinomial Model Estimates with Weights.csv")
ifed2 <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Adjusted QGCOMP RR Model Estimates with Weights.csv")

colnames(ifed1)<-tolower(colnames(ifed1))
colnames(ifed2)<-tolower(colnames(ifed2))

#### Edit dataframes
ifed1$x = ifelse(grepl("Inc",ifed1$y.level), "Delayed Increase",
                 ifelse(grepl("Dec",ifed1$y.level), "Decrease", "Rapid Increase"))
ifed2$x = ifelse(grepl("inc",ifed2$outcome), "Delayed Increase",
                 ifelse(grepl("dec",ifed2$outcome), "Decrease", "Rapid Increase"))
ifed2$or = ifed2$est
ifed2$term <- "Mixture"

##Merge and output datasets
ifed3 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(ifed1,ifed2))

## Order variables: exposures
ifed3$variable <- ifelse(ifed3$term == "bpa_obrien_gm_ln", "BPA",
                         ifelse(ifed3$term == "bp_3_obrien_gm_ln", "Benzophenone-3",
                                ifelse(ifed3$term == "dcp_24_obrien_gm_ln", "2,4-dichlorophenol",
                                       ifelse(ifed3$term == "dcp_25_obrien_gm_ln", "2,5-dichlorophenol",
                                              ifelse(ifed3$term == "m_pb_obrien_gm_ln", "Methyl paraben",
                                                     ifelse(ifed3$term == "p_pb_obrien_gm_ln", "Propyl paraben",
                                                            ifelse(ifed3$term == "tcs_obrien_gm_ln", "Triclosan",
                                                                   ifelse(ifed3$term == "Mixture", "Mixture", "Mixture"))))))))

exposure.labels <- c(exposure.labels,"Mixture")
ifed3$variable <- factor(ifed3$variable, levels = c(exposure.labels))
ifed3 <- ifed3[order(ifed3$variable),]

## Order variables: outcomes
ifed3$x <- factor(ifed3$x, levels = c('Rapid Increase','Delayed Increase','Decrease'))
ifed3 <- ifed3[order(ifed3$x),]

## Plot
ifed.plot.1 <- ggplot(ifed3, aes(x=as.numeric(or), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=variable, col=factor(x))) + 
  geom_point(position=position_dodge(width=0.5), shape = 19, size = 2) +
  geom_errorbar(width=0.2, position=position_dodge(width=0.5), linewidth = 1.5) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  labs(y = "", x = "Relative Risk Ratio (95% CI)\nfor a 1-IQR Increase in Phenol Concentration") + geom_vline(xintercept=1.0,linetype=2) +
  scale_y_discrete(labels = rev(exposure.labels),
                   limits = rev(exposure.labels)) +
  scale_x_continuous(trans=log10_trans(),breaks=c(0,0.25,0.5,1,2,4),limits=c(0.2,4))+
  facet_wrap(~x, nrow=1) +
  scale_colour_manual(values=c("#7570B3","#D95F02","#1B9E77"), name='x') +
  theme(text = element_text(size=25, color="black"),
        legend.position = 'none',
        panel.grid.minor = element_blank(),               #removes minor grid lines
        panel.grid.major = element_blank()) 

ggsave(ifed.plot.1, filename="J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Adjusted RRs for Primary Analysis.jpg", width = 18, height = 7, units = "in", dpi=600)




######################
## Mean Differences ##
######################

ifed4 <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\BMI z-score estimates.csv")
ifed5 <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Longitudinal QGCOMP Estimates.csv")

colnames(ifed4)<-tolower(colnames(ifed4))
colnames(ifed5)<-tolower(colnames(ifed5))

#### Edit dataframes
ifed5$estimate <- ifed5$est
ifed5$exposure <- "Mixture"
ifed4$lcl <- ifed4$lclmean
ifed4$ucl <- ifed4$uclmean

##Merge and output datasets
ifed6 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(ifed4,ifed5))

## Order variables: exposures
ifed6$variable <- ifelse(ifed6$exposure == "bpa_obrien_gm_ln", "BPA",
                         ifelse(ifed6$exposure == "bp_3_obrien_gm_ln", "Benzophenone-3",
                                ifelse(ifed6$exposure == "dcp_24_obrien_gm_ln", "2,4-dichlorophenol",
                                       ifelse(ifed6$exposure == "dcp_25_obrien_gm_ln", "2,5-dichlorophenol",
                                              ifelse(ifed6$exposure == "m_pb_obrien_gm_ln", "Methyl paraben",
                                                     ifelse(ifed6$exposure == "p_pb_obrien_gm_ln", "Propyl paraben",
                                                            ifelse(ifed6$exposure == "tcs_obrien_gm_ln", "Triclosan",
                                                                   ifelse(ifed6$exposure == "Mixture", "Mixture", "Mixture"))))))))

ifed6$variable <- factor(ifed6$variable, levels = c(exposure.labels))
ifed6 <- ifed6[order(ifed6$variable),]

## Plot
ifed7 <- subset(ifed6, label == "bmiz" | label == "Primary Analysis: Adjusted")
ifed.plot.2 <- ggplot(ifed7, aes(x=as.numeric(estimate), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=variable, col="black")) + 
  geom_point(position=position_dodge(width=0.5), shape = 19, size = 2, col="black") +
  geom_errorbar(width=0.2, position=position_dodge(width=0.5), linewidth = 1.5, col="black") +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  labs(y = "", x = "Mean Difference in BMI Z-Score (95% CI)\nfor a 1-IQR Increase in Phenol Concentration") + geom_vline(xintercept=0,linetype=2,alpha=0.2) +
  scale_y_discrete(labels = rev(exposure.labels),
                   limits = rev(exposure.labels)) +
  scale_x_continuous(breaks=c(-0.5,0,0.5,1.0),limits=c(-0.5,1.0)) +
  theme(text = element_text(size=25, color="black"),
        legend.position = 'none',
        panel.grid.minor = element_blank(),               
        panel.grid.major = element_blank())

ggsave(ifed.plot.2, filename="J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Adjusted Betas for Primary Analysis.jpg", width = 10, height = 7, units = "in", dpi=600)


#################################
## Plot of exposure prevalence ##
#################################

install.packages('ggmice')
library(ggmice)
exposure.labels <- c('BPA', 'Benzophenone-3', '2,4-dichlorophenol', '2,5-dichlorophenol', 
                     'Methyl paraben', 'Propyl paraben', 'Triclosan')


imp.dat.wide<-read.csv("U:\\IFED\\data\\Imputed dataset with outcomes wide.csv")
a = imp.dat.wide %>% select(c("ifedid", ".imp", ".id", all_of(exposures.obrien)))
b = a %>% pivot_longer(cols = all_of(exposures.obrien),
                       names_to = "analyte")
b$variable <- ifelse(b$analyte == "bpa_obrien_gm", "BPA",
                     ifelse(b$analyte == "bp_3_obrien_gm", "Benzophenone-3",
                            ifelse(b$analyte == "dcp_24_obrien_gm", "2,4-dichlorophenol",
                                   ifelse(b$analyte == "dcp_25_obrien_gm", "2,5-dichlorophenol",
                                          ifelse(b$analyte == "m_pb_obrien_gm", "Methyl paraben",
                                                 ifelse(b$analyte == "p_pb_obrien_gm", "Propyl paraben",
                                                        ifelse(b$analyte == "tcs_obrien_gm", "Triclosan", NA)))))))
b$lod <- ifelse(b$analyte == "bpa_obrien_gm", "BPA\n96%>LOD",
                ifelse(b$analyte == "bp_3_obrien_gm", "Benzophenone-3\n98%>LOD",
                       ifelse(b$analyte == "dcp_24_obrien_gm", "2,4-dichlorophenol\n95%>LOD",
                              ifelse(b$analyte == "dcp_25_obrien_gm", "2,5-dichlorophenol\n91%>LOD",
                                     ifelse(b$analyte == "m_pb_obrien_gm", "Methyl paraben\n99%>LOD",
                                            ifelse(b$analyte == "p_pb_obrien_gm", "Propyl paraben\n98%>LOD",
                                                   ifelse(b$analyte == "tcs_obrien_gm", "Triclosan\n64%>LOD", NA)))))))

new.lab<-c("BPA\n96%>LOD", "Benzophenone-3\n98%>LOD", "2,4-dichlorophenol\n95%>LOD",
           "2,5-dichlorophenol\n91%>LOD", "Methyl paraben\n99%>LOD", "Propyl paraben\n98%>LOD",
           "Triclosan\n64%>LOD")

b$lod <- factor(b$lod, levels = c(new.lab))
b <- b[order(b$lod),]

explot = ggmice(data = b, aes(x = value, fill = lod)) + 
  geom_boxplot(aes(x = value, fill = lod), color='black') + 
  scale_x_log10(labels = scales::comma) + xlab("Biomarker Concentration (ng/mL)") + 
  # geom_vline(data = b, aes(xintercept=LOD), linetype = "dashed") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_y_continuous(position = "right") +
  facet_grid(vars(lod), scales = "free", switch = "y") + 
  theme_bw() + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 18),
        legend.position = "none",
        strip.text.y.left = element_text(angle=0))

explot

ggsave(explot, filename = "J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Box Plot of Exposure Dist.jpg", width = 7, height = 8, units = "in", dpi=600)




