# Programmer: Danielle Stevens
# Last updated on: 11/11/2026
# 
# Purpose of code: This code reads in the IFED data,
#                  examines and plots outcomes, 
#                  and examines and plots the IFED growth trajectories using gmm
#                            


library(gtable)
library(gridExtra)
library(reprex)
library(egg)
library(extrafont)
library(ggplotify)
#font_import()
#y
windowsFonts(Arial = windowsFont("Arial"))
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()                       #vector of font family names
library(ggplot2)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
library(haven)
library(grid)
library(gridExtra)
library(scales)
library(stringr)
library(lcmm)
library(nnet)
library(stringr)
library(plyr)
library(tidyr)
library(dplyr) 
library(tidyverse)
library(pheatmap)
library(dichromat)
library(mice)
library(miceadds)
library(splines)
library(qgcomp)
library(mde)
library(data.table)
library(tableone)
library(ggpubr)
library(egg)
library(patchwork)

packageVersion('qgcomp') #2.9.0
packageVersion('mice') #3.14.0

#Read in sas dataset
long <- read_sas("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\ifed_long.sas7bdat")
wide <- read_sas("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\ifed_wide.sas7bdat")

#set working directory
setwd("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results")

#Get names to lowercase
names(long) <- tolower(names(long))
names(wide) <- tolower(names(wide))


###################
## Plot outcomes ##
###################
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Accent")[c(1,2,4)]
pal3 <- brewer.pal(8,"Set1")[c(1)]

my_palette <- c(pal1,pal2,pal3)

summary(long$bmiz)

long$n_infant_weight <- with(long, reorder(n_infant_weight, visit))
infant_weight <- long[, c('visit','ageinweeks','infant_weight',"n_infant_weight")]
infant_weight <- infant_weight[complete.cases(infant_weight),]
weightplot <- ggplot(data = infant_weight, aes(x = ageinweeks, y = infant_weight, color=as.factor(n_infant_weight))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant weight (kg)") + 
  scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(labels=label_number(accuracy = 1)) + 
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
weightplot

long$n_height <- with(long, reorder(n_height, visit))
height <- long[, c('visit','ageinweeks','height',"n_height")]
height <- height[complete.cases(height),]
heightplot <- ggplot(data = height, aes(x = ageinweeks, y = height, color=as.factor(n_height))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant length (cm)") + 
  scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(labels=label_number(accuracy = 1)) + 
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
heightplot

long$n_waz <- with(long, reorder(n_waz, visit))
waz <- long[, c('visit','ageinweeks','waz',"n_waz")]
waz <- waz[complete.cases(waz),]
wazplot <- ggplot(data = waz, aes(x = ageinweeks, y = waz, color=as.factor(n_waz))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant weight z-score") + 
  scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(labels=label_number(accuracy = 1)) + 
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
wazplot

long$n_haz <- with(long, reorder(n_haz, visit))
haz <- long[, c('visit','ageinweeks','haz',"n_haz")]
haz <- haz[complete.cases(haz),]
hazplot <- ggplot(data = haz, aes(x = ageinweeks, y = haz, color=as.factor(n_haz))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant length z-score") + 
  scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(labels=label_number(accuracy = 1)) + 
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
hazplot

long$n_bmiz <- with(long, reorder(n_bmiz, visit))
bmiz <- long[, c('visit','ageinweeks','bmiz',"n_bmiz")]
bmiz <- bmiz[complete.cases(bmiz),]
bmizplot <- ggplot(data = bmiz, aes(x = ageinweeks, y = bmiz, color=as.factor(n_bmiz))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant BMI z-score") + 
  scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(labels=label_number(accuracy = 1)) + 
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
bmizplot


comboplot = grid.arrange(weightplot, heightplot, wazplot, hazplot, ncol=2, widths = c(4,4))
comboplot

ggsave(plot = comboplot, filename = "wthtplot.jpg", width = 16, height = 12)


## Run BMI zscore plot with modifications suggested by co-authors on 08/22/2023
long$n_bmiz <- with(long, reorder(n_bmiz, visit))
bmiz <- long[, c('visit','ageinweeks','bmiz',"n_bmiz")]
bmiz <- bmiz[complete.cases(bmiz),]
bmizplot <- ggplot(data = bmiz, aes(x = ageinweeks, y = bmiz, color=as.factor(n_bmiz))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant BMI z-score") + 
  scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -1.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 1.96, linetype = 'dashed', color = 'grey') +
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
bmizplot
ggsave(plot = bmizplot, filename = "bmizplot.jpg", width = 8, height = 6)
# One potential outlier at -4 for bmiz


## Run BMI zscore plot with modifications suggested by co-authors on 08/22/2023 & then revmove color on 11/15/2023
long$n_bmiz <- with(long, reorder(n_bmiz, visit))
bmiz <- long[, c('visit','ageinweeks','bmiz',"n_bmiz")]
bmiz <- bmiz[complete.cases(bmiz),]
bmizplot <- ggplot(data = bmiz, aes(x = ageinweeks, y = bmiz, color=as.factor(n_bmiz))) + 
  geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant BMI z-score") + 
  scale_color_manual(values=c(rep('black',12))) + labs(color = 'Visit week (n)') +
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -1.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 1.96, linetype = 'dashed', color = 'grey') +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
bmizplot
ggsave(plot = bmizplot, filename = "bmizplot_nocol.jpg", width = 8, height = 6)




## Re-run for each group (males, females, breastmilk, cows-milk) ##
male <- subset(long, long$childsgender==2)
 male$n_male_weight <- with(male, reorder(n_male_weight, visit))
 male$n_male_height <- with(male, reorder(n_male_height, visit))
 male$n_male_waz <- with(male, reorder(n_male_waz, visit))
 male$n_male_haz <- with(male, reorder(n_male_haz, visit))
 male$n_male_bmiz <- with(male, reorder(n_male_bmiz, visit))

female <- subset(long, long$childsgender==1)
 female$n_female_weight <- with(female, reorder(n_female_weight, visit))
 female$n_female_height <- with(female, reorder(n_female_height, visit))
 female$n_female_waz <- with(female, reorder(n_female_waz, visit))
 female$n_female_haz <- with(female, reorder(n_female_haz, visit))
 female$n_female_bmiz <- with(female, reorder(n_female_bmiz, visit))

breastmilk <- subset(long, long$feedinggroup==1)
 breastmilk$n_breastmilk_weight <- with(breastmilk, reorder(n_breastmilk_weight, visit))
 breastmilk$n_breastmilk_height <- with(breastmilk, reorder(n_breastmilk_height, visit))
 breastmilk$n_breastmilk_waz <- with(breastmilk, reorder(n_breastmilk_waz, visit))
 breastmilk$n_breastmilk_haz <- with(breastmilk, reorder(n_breastmilk_haz, visit))
 breastmilk$n_breastmilk_bmiz <- with(breastmilk, reorder(n_breastmilk_bmiz, visit))

cowsmilk <- subset(long, long$feedinggroup==2)
 cowsmilk$n_cowsmilk_weight <- with(cowsmilk, reorder(n_cowsmilk_weight, visit))
 cowsmilk$n_cowsmilk_height <- with(cowsmilk, reorder(n_cowsmilk_height, visit))
 cowsmilk$n_cowsmilk_waz <- with(cowsmilk, reorder(n_cowsmilk_waz, visit))
 cowsmilk$n_cowsmilk_haz <- with(cowsmilk, reorder(n_cowsmilk_haz, visit))
 cowsmilk$n_cowsmilk_bmiz <- with(cowsmilk, reorder(n_cowsmilk_bmiz, visit))

plotfunction <- function(input,weightvar,heightvar,wazvar,hazvar,bmizvar,output1,output2){

  data <- input
  
  weightplot <- ggplot(data = data, aes(x = ageinweeks, y = infant_weight, color=as.factor(weightvar))) + 
    geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant weight (kg)") + 
    scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
    scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
    scale_y_continuous(labels=label_number(accuracy = 1)) + 
    theme_bw() +  
    theme(text = element_text(size=19),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
  weightplot
  
  heightplot <- ggplot(data = data, aes(x = ageinweeks, y = height, color=as.factor(heightvar))) + 
    geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant length (cm)") + 
    scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
    scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
    scale_y_continuous(labels=label_number(accuracy = 1)) + 
    theme_bw() +  
    theme(text = element_text(size=19),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
  heightplot
  
   wazplot <- ggplot(data = data, aes(x = ageinweeks, y = waz, color=as.factor(wazvar))) + 
    geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant weight z-score") + 
    scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
    scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
    scale_y_continuous(labels=label_number(accuracy = 1)) + 
     theme_bw() +  
     theme(text = element_text(size=19),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
   wazplot
  
  hazplot <- ggplot(data = data, aes(x = ageinweeks, y = haz, color=as.factor(hazvar))) + 
    geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant length z-score") + 
    scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
    scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
    scale_y_continuous(labels=label_number(accuracy = 1)) + 
    theme_bw() +  
    theme(text = element_text(size=19),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
  hazplot
  
   bmizplot <- ggplot(data = data, aes(x = ageinweeks, y = bmiz, color=as.factor(bmizvar))) + 
    geom_point(size=2) + geom_jitter() + xlab("Age (weeks)\n ") + ylab(" \nInfant BMI z-score") + 
    scale_color_manual(values=my_palette) + labs(color = 'Visit week (n)') +
    scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
    scale_y_continuous(labels=label_number(accuracy = 1)) + 
     theme_bw() +  
     theme(text = element_text(size=19),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 
   
  comboplot = grid.arrange(weightplot, heightplot, wazplot, hazplot, ncol=2, widths = c(4,4))
  ggsave(plot = comboplot, filename = output1, width = 16, height = 12)
  ggsave(plot = bmizplot, filename = output2, width = 8, height = 6)
  
  bmizplot
  }

maleplot<-plotfunction(male, male$n_male_weight, male$n_male_height, male$n_male_waz, male$n_male_haz, male$n_male_bmiz, "male_wthtplot.jpg", "male_bmizplot.jpg")
femaleplot<-plotfunction(female, female$n_female_weight, female$n_female_height, female$n_female_waz, female$n_female_haz, female$n_female_bmiz, "female_wthtplot.jpg", "female_bmizplot.jpg")
breastmilkplot<-plotfunction(breastmilk, breastmilk$n_breastmilk_weight, breastmilk$n_breastmilk_height, breastmilk$n_breastmilk_waz, breastmilk$n_breastmilk_haz, breastmilk$n_breastmilk_bmiz, "breastmilk_wthtplot.jpg", "breastmilk_bmizplot.jpg")
cowsmilkplot<-plotfunction(cowsmilk, cowsmilk$n_cowsmilk_weight, cowsmilk$n_cowsmilk_height, cowsmilk$n_cowsmilk_waz, cowsmilk$n_cowsmilk_haz, cowsmilk$n_cowsmilk_bmiz, "cowsmilk_wthtplot.jpg", "cowsmilk_bmizplot.jpg")

## Merge all bmizplots together
comboplot = grid.arrange(maleplot + ggtitle("Males"),
                         femaleplot + ggtitle("Females"),
                         breastmilkplot + ggtitle("Breastmilk diet"),
                         cowsmilkplot + ggtitle("Cows-milk diet"),
                         ncol=2, widths=c(4,4))
  
ggsave(plot = comboplot, filename = "subgroups_bmizplot.jpg", width = 16, height = 12)



######################################
## Estimate trajectories using GMM  ##
######################################

## Convert z-scores to Centered z-scores
zscore_data <- wide[,c('ifedid','mean_bmiz','mean_waz','mean_haz')]
long2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(long,zscore_data))
long2$norm_bmiz <- long2$bmiz - long2$mean_bmiz
long2$norm_waz <- long2$waz - long2$mean_waz
long2$norm_haz <- long2$haz - long2$mean_haz

ggplot(data=long2,aes(x=ageinweeks,y=norm_bmiz)) + geom_line()
ggplot(data=long2,aes(x=ageinweeks,y=norm_waz)) + geom_line()
ggplot(data=long2,aes(x=ageinweeks,y=norm_haz)) + geom_line()

## Count number of BMI measures
dat.bmiz <- long2[,c('ifedid','bmiz')][complete.cases(long2[,c('ifedid','bmiz')]),]
data_table <- data.table(dat.bmiz)  
n_bmiz <- data_table[ , .(count = length(unique(bmiz))), by = ifedid]
mean_bmiz <- wide[,c('ifedid', 'mean_bmiz')]
dat.bmiz <- Reduce(function(x, y) merge(x, y, all=TRUE), list(mean_bmiz, n_bmiz))

CreateTableOne(data=dat.bmiz,
               vars='mean_bmiz',
               factorVars='count',
               strata="count")


## Subset datasets by sex
male <- subset(long2, long2$childsgender==2)
female <- subset(long2, long2$childsgender==1)
male.wide <- subset(wide, wide$childsgender==2)
length(male.wide$ifedid) #n=101
female.wide <- subset(wide, wide$childsgender==1)
length(female.wide$ifedid) #n=98

## Subset datasets by diet
breastmilk <- subset(long2, long2$feedinggroup==1)
cowsmilk <- subset(long2, long2$feedinggroup==2)
breastmilk.wide <- subset(wide, wide$feedinggroup==1)
length(breastmilk.wide$ifedid) #n=90
cowsmilk.wide <- subset(wide, wide$feedinggroup==2)
length(cowsmilk.wide$ifedid) #n=109



#########################
## Run for BMI z-score ##
## By sex              ##
#########################

## Run gaussian model for Centered bmi z-score with ageinweeks as the time variable.
# We allow for random intercepts and random effects and run stratified by infant sex
# We initialize using the gridsearch function based on a 1-group model

set.seed(333)


## Test splines models: b-splines ##
## 2 knots: 5.85714, 19 weeks split the age distribution at 33 and 66th percentiles ##
## 3 knots: 4, 11.4286, 23 weeks split the age distribution at 25, 50, 75th percentiles ##
## 4 knots: 3.42857, 7.85714, 15.8571, 24.2857 weeks split the age distribution at 20, 40, 60, and 80th percentiles ##

## Males first ##
# Cubic B-splines #
males.m1.cs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=male)
males.m1.cs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=3) , random =~1, subject = "ifedid", data=male)
males.m1.cs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=3) , random =~1, subject = "ifedid", data=male)

males.m1.cs.4$pred$resid_ss_sq <- males.m1.cs.4$pred$resid_ss^2
males.m1.cs.4.sse <- sum(males.m1.cs.4$pred$resid_ss_sq) 

males.m1.cs.3$pred$resid_ss_sq <- males.m1.cs.3$pred$resid_ss^2
males.m1.cs.3.sse <- sum(males.m1.cs.3$pred$resid_ss_sq) 

males.m1.cs.2$pred$resid_ss_sq <- males.m1.cs.2$pred$resid_ss^2
males.m1.cs.2.sse <- sum(males.m1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
males.m1.qs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=2) , random =~1, subject = "ifedid", data=male)
males.m1.qs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=2) , random =~1, subject = "ifedid", data=male)
males.m1.qs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=2) , random =~1, subject = "ifedid", data=male)

males.m1.qs.4$pred$resid_ss_sq <- males.m1.qs.4$pred$resid_ss^2
males.m1.qs.4.sse <- sum(males.m1.qs.4$pred$resid_ss_sq) 

males.m1.qs.3$pred$resid_ss_sq <- males.m1.qs.3$pred$resid_ss^2
males.m1.qs.3.sse <- sum(males.m1.qs.3$pred$resid_ss_sq) 

males.m1.qs.2$pred$resid_ss_sq <- males.m1.qs.2$pred$resid_ss^2
males.m1.qs.2.sse <- sum(males.m1.qs.2$pred$resid_ss_sq) 

# Compare SSEs #
males.m1.cs.4.sse
males.m1.cs.3.sse
males.m1.cs.2.sse

males.m1.qs.4.sse
males.m1.qs.3.sse
males.m1.qs.2.sse
## Lowest SSE is from: cubic model with 4 knots##


## Does this agree with results for females? ##
# Cubic B-splines #
females.m1.cs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=female)
females.m1.cs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=3) , random =~1, subject = "ifedid", data=female)
females.m1.cs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=3) , random =~1, subject = "ifedid", data=female)

females.m1.cs.4$pred$resid_ss_sq <- females.m1.cs.4$pred$resid_ss^2
females.m1.cs.4.sse <- sum(females.m1.cs.4$pred$resid_ss_sq) 

females.m1.cs.3$pred$resid_ss_sq <- females.m1.cs.3$pred$resid_ss^2
females.m1.cs.3.sse <- sum(females.m1.cs.3$pred$resid_ss_sq) 

females.m1.cs.2$pred$resid_ss_sq <- females.m1.cs.2$pred$resid_ss^2
females.m1.cs.2.sse <- sum(females.m1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
females.m1.qs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=2) , random =~1, subject = "ifedid", data=female)
females.m1.qs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=2) , random =~1, subject = "ifedid", data=female)
females.m1.qs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=2) , random =~1, subject = "ifedid", data=female)

females.m1.qs.4$pred$resid_ss_sq <- females.m1.qs.4$pred$resid_ss^2
females.m1.qs.4.sse <- sum(females.m1.qs.4$pred$resid_ss_sq) 

females.m1.qs.3$pred$resid_ss_sq <- females.m1.qs.3$pred$resid_ss^2
females.m1.qs.3.sse <- sum(females.m1.qs.3$pred$resid_ss_sq) 

females.m1.qs.2$pred$resid_ss_sq <- females.m1.qs.2$pred$resid_ss^2
females.m1.qs.2.sse <- sum(females.m1.qs.2$pred$resid_ss_sq) 

females.m1.cs.4.sse
females.m1.cs.3.sse
females.m1.cs.2.sse

females.m1.qs.4.sse
females.m1.qs.3.sse
females.m1.qs.2.sse
## Same model has the best fit so let's use that: Cubic b-spline with 4 knots ##




## Determine Number of Groups ##
## Males ##
## NOTE: We get almost exactly (some differences for male.m3.s) even if we use splines at points that split the male age distribution ##
male.m1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=male)
male.m2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=male, ng=2, B=male.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=male.m1.s)
male.m3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=male, ng=3, B=male.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=male.m1.s)
male.m4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=male, ng=4, B=male.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=male.m1.s)
male.m5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=male, ng=5, B=male.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=male.m1.s)
male.m6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=male, ng=6, B=male.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=male.m1.s)
male.tab.s<-as.data.frame(summarytable(male.m1.s, male.m2.s, male.m3.s, male.m4.s, male.m5.s, male.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(male.m1.s, male.m2.s, male.m3.s, male.m4.s, male.m5.s, male.m6.s, which = c("AIC", "BIC", "entropy"))  
plot(male.m4.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE) 

## Number of participants/group?
male.wide <- Reduce(function(x, y) merge(x, y, all=TRUE), list(male.wide, male.m4.s$pprob))
table(male.wide$class)

## Plots
# Output predictions by class
male.pred <- data.frame(ageinweeks=seq(0,28))
male.fit <- predictY(male.m4.s,draws=TRUE,newdata=male.pred,var.time="ageinweeks")
male.fit.pred <- cbind(male.fit$time, male.fit$pred)
male.fit <- male.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
male.fit <- male.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
male.fit$norm_bmiz = male.fit$Ypred
male.fit$lower.Ypred <- ifelse(male.fit$lower.Ypred > -2.5, male.fit$lower.Ypred, NA)
male.fit$label <- ifelse(male.fit$class == 'class1', "Stable, n=57",
                           ifelse(male.fit$class =='class2', "Rapid Increase, n=15",
                                  ifelse(male.fit$class =='class3', "Decrease, n=21", "Delayed Increase, n=8")))

# Overlay with original data
male <- Reduce(function(x, y) merge(x, y, all=TRUE), list(male, male.m4.s$pprob))
male$class = str_c("Class ", male$class)
male$label <- ifelse(male$class == 'Class 1', "Stable, n=57",
                         ifelse(male$class =='Class 2', "Rapid Increase, n=15",
                                ifelse(male$class =='Class 3', "Decrease, n=21", "Delayed Increase, n=8")))


male.plot <- ggplot(male, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=male.fit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=male.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 


## Females ##
female.m1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=female)
female.m2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female, ng=2, B=female.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female.m1.s)
female.m3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female, ng=3, B=female.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female.m1.s)
female.m4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female, ng=4, B=female.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female.m1.s)
female.m5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female, ng=5, B=female.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female.m1.s)
female.m6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female, ng=6, B=female.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female.m1.s)
female.tab.s<-as.data.frame(summarytable(female.m1.s, female.m2.s, female.m3.s, female.m4.s, female.m5.s, female.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(female.m1.s, female.m2.s, female.m3.s, female.m4.s, female.m5.s, female.m6.s, which = c("AIC", "BIC", "entropy"))  
plot(female.m4.s, which="fit", var.time="ageinweeks", marg=TRUE, shades = TRUE) 

## Number of participants/group?
female.wide <- Reduce(function(x, y) merge(x, y, all=TRUE), list(female.wide, female.m4.s$pprob))
table(female.wide$class)

## Plots
# Output predictions by class
female.pred <- data.frame(ageinweeks=seq(0,36))
female.fit <- predictY(female.m4.s,draws=TRUE,newdata=female.pred,var.time="ageinweeks")
female.fit.pred <- cbind(female.fit$time, female.fit$pred)
female.fit <- female.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
female.fit <- female.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
female.fit$norm_bmiz = female.fit$Ypred
female.fit$label <- ifelse(female.fit$class == 'class1', "Stable, n=30",
                           ifelse(female.fit$class =='class2', "Decrease, n=24",
                                  ifelse(female.fit$class =='class3', "Rapid Increase, n=30", "Delayed Increase, n=14")))
                                         
# Overlay with original data
female <- Reduce(function(x, y) merge(x, y, all=TRUE), list(female, female.m4.s$pprob))
female$class = str_c("Class ", female$class)
female$label <- ifelse(female$class == 'Class 1', "Stable, n=30",
                           ifelse(female$class =='Class 2', "Decrease, n=24",
                                  ifelse(female$class =='Class 3', 'Rapid Increase, n=30', "Delayed Increase, n=14")))

female.plot <- ggplot(female, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=female.fit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=female.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 


## Output table of model fit statistics
male.tab.s$data <- "Males"
female.tab.s$data <- "Females"

tab.s <- rbind(male.tab.s,female.tab.s)
write.csv(tab.s, "tabforGMM_sex.csv")
## Let's pick the 4 group b/c it maxs entropy while min BIC, and we have good differentiation of groups ##
## We have only 7% in one of the groups for males but the trajectories are essentially the same regardless of sex ##




##################################################################
## Sensitivity analysis for females limiting to before 28 weeks ##
##################################################################

female28 <- subset(long2, long2$childsgender==1 & long2$ageinweeks <= 28)
female28.wide <- subset(wide, wide$childsgender==1)

length(unique(female28$ifedid)) #n=98

## Females ##
female28.m1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=female28)
female28.m2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female28, ng=2, B=female28.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female28.m1.s)
female28.m3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female28, ng=3, B=female28.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female28.m1.s)
female28.m4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female28, ng=4, B=female28.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female28.m1.s)
female28.m5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female28, ng=5, B=female28.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female28.m1.s)
female28.m6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=female28, ng=6, B=female28.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=female28.m1.s)
female28.tab.s<-as.data.frame(summarytable(female28.m1.s, female28.m2.s, female28.m3.s, female28.m4.s, female28.m5.s, female28.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(female28.m1.s, female28.m2.s, female28.m3.s, female28.m4.s, female28.m5.s, female28.m6.s, which = c("AIC", "BIC", "entropy"))  
plot(female28.m4.s, which="fit", var.time="ageinweeks", marg=TRUE, shades = TRUE) 


## Number of participants/group?
female28.wide <- Reduce(function(x, y) merge(x, y, all=TRUE), list(female28.wide, female28.m4.s$pprob))
table(female28.wide$class)

## Plots
# Output predictions by class
female28.pred <- data.frame(ageinweeks=seq(0,32))
female28.fit <- predictY(female28.m4.s,draws=TRUE,newdata=female28.pred,var.time="ageinweeks")
female28.fit.pred <- cbind(female28.fit$time, female28.fit$pred)
female28.fit <- female28.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
female28.fit <- female28.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
female28.fit$norm_bmiz = female28.fit$Ypred
female28.fit$label <- ifelse(female28.fit$class == 'class1', "Delayed Increase, n=12",
                           ifelse(female28.fit$class =='class2', "Decrease, n=30",
                                  ifelse(female28.fit$class =='class3', "Rapid Increase, n=22", "Stable, n=34")))

# Overlay with original data
female28 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(female28, female28.m4.s$pprob))
female28$class = str_c("Class ", female28$class)
female28$label <- ifelse(female28$class == 'Class 1', "Delayed Increase, n=12",
                       ifelse(female28$class =='Class 2', "Decrease, n=30",
                              ifelse(female28$class =='Class 3', "Rapid Increase, n=22", "Stable, n=34")))

female28.plot <- ggplot(female28, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=female28.fit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=female28.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 


## Output table of model fit statistics
female28.tab.s$data <- "Females"
write.csv(female28.tab.s, "tabforGMM_femalesens28wks.csv")


#########################
## Run for BMI z-score ##
## By infant diet      ##
#########################

## Run gaussian model for Centered bmi z-score with ageinweeks as the time variable.
# We allow for random intercepts and random effects and run stratified by infant sex
# We initialize using the gridsearch function based on a 1-group model

set.seed(333)


## breastmilk first ##
# Cubic B-splines #
breastmilk.m1.cs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=breastmilk)
breastmilk.m1.cs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=3) , random =~1, subject = "ifedid", data=breastmilk)
breastmilk.m1.cs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=3) , random =~1, subject = "ifedid", data=breastmilk)

breastmilk.m1.cs.4$pred$resid_ss_sq <- breastmilk.m1.cs.4$pred$resid_ss^2
breastmilk.m1.cs.4.sse <- sum(breastmilk.m1.cs.4$pred$resid_ss_sq) 

breastmilk.m1.cs.3$pred$resid_ss_sq <- breastmilk.m1.cs.3$pred$resid_ss^2
breastmilk.m1.cs.3.sse <- sum(breastmilk.m1.cs.3$pred$resid_ss_sq) 

breastmilk.m1.cs.2$pred$resid_ss_sq <- breastmilk.m1.cs.2$pred$resid_ss^2
breastmilk.m1.cs.2.sse <- sum(breastmilk.m1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
breastmilk.m1.qs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=2) , random =~1, subject = "ifedid", data=breastmilk)
breastmilk.m1.qs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=2) , random =~1, subject = "ifedid", data=breastmilk)
breastmilk.m1.qs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=2) , random =~1, subject = "ifedid", data=breastmilk)

breastmilk.m1.qs.4$pred$resid_ss_sq <- breastmilk.m1.qs.4$pred$resid_ss^2
breastmilk.m1.qs.4.sse <- sum(breastmilk.m1.qs.4$pred$resid_ss_sq) 

breastmilk.m1.qs.3$pred$resid_ss_sq <- breastmilk.m1.qs.3$pred$resid_ss^2
breastmilk.m1.qs.3.sse <- sum(breastmilk.m1.qs.3$pred$resid_ss_sq) 

breastmilk.m1.qs.2$pred$resid_ss_sq <- breastmilk.m1.qs.2$pred$resid_ss^2
breastmilk.m1.qs.2.sse <- sum(breastmilk.m1.qs.2$pred$resid_ss_sq) 

# Compare SSEs #
breastmilk.m1.cs.4.sse
breastmilk.m1.cs.3.sse
breastmilk.m1.cs.2.sse

breastmilk.m1.qs.4.sse
breastmilk.m1.qs.3.sse
breastmilk.m1.qs.2.sse
## Lowest SSE is from: cubic model with 4 knots##


## Does this agree with results for cowsmilk? ##
# Cubic B-splines #
cowsmilk.m1.cs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=cowsmilk)
cowsmilk.m1.cs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=3) , random =~1, subject = "ifedid", data=cowsmilk)
cowsmilk.m1.cs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=3) , random =~1, subject = "ifedid", data=cowsmilk)

cowsmilk.m1.cs.4$pred$resid_ss_sq <- cowsmilk.m1.cs.4$pred$resid_ss^2
cowsmilk.m1.cs.4.sse <- sum(cowsmilk.m1.cs.4$pred$resid_ss_sq) 

cowsmilk.m1.cs.3$pred$resid_ss_sq <- cowsmilk.m1.cs.3$pred$resid_ss^2
cowsmilk.m1.cs.3.sse <- sum(cowsmilk.m1.cs.3$pred$resid_ss_sq) 

cowsmilk.m1.cs.2$pred$resid_ss_sq <- cowsmilk.m1.cs.2$pred$resid_ss^2
cowsmilk.m1.cs.2.sse <- sum(cowsmilk.m1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
cowsmilk.m1.qs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=2) , random =~1, subject = "ifedid", data=cowsmilk)
cowsmilk.m1.qs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=2) , random =~1, subject = "ifedid", data=cowsmilk)
cowsmilk.m1.qs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=2) , random =~1, subject = "ifedid", data=cowsmilk)

cowsmilk.m1.qs.4$pred$resid_ss_sq <- cowsmilk.m1.qs.4$pred$resid_ss^2
cowsmilk.m1.qs.4.sse <- sum(cowsmilk.m1.qs.4$pred$resid_ss_sq) 

cowsmilk.m1.qs.3$pred$resid_ss_sq <- cowsmilk.m1.qs.3$pred$resid_ss^2
cowsmilk.m1.qs.3.sse <- sum(cowsmilk.m1.qs.3$pred$resid_ss_sq) 

cowsmilk.m1.qs.2$pred$resid_ss_sq <- cowsmilk.m1.qs.2$pred$resid_ss^2
cowsmilk.m1.qs.2.sse <- sum(cowsmilk.m1.qs.2$pred$resid_ss_sq) 

cowsmilk.m1.cs.4.sse
cowsmilk.m1.cs.3.sse
cowsmilk.m1.cs.2.sse

cowsmilk.m1.qs.4.sse
cowsmilk.m1.qs.3.sse
cowsmilk.m1.qs.2.sse
## Same model has the best fit so let's use that: Cubic b-spline with 4 knots ##



## Determine Number of Groups ##
## breastmilk ##
count(unique(breastmilk$ifedid)) #90
breastmilk.m1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=breastmilk)
breastmilk.m2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=breastmilk, ng=2, B=breastmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=breastmilk.m1.s)
breastmilk.m3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=breastmilk, ng=3, B=breastmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=breastmilk.m1.s)
breastmilk.m4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=breastmilk, ng=4, B=breastmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=breastmilk.m1.s)
breastmilk.m5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=breastmilk, ng=5, B=breastmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=breastmilk.m1.s)
breastmilk.m6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=breastmilk, ng=6, B=breastmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=breastmilk.m1.s)
breastmilk.tab.s<-as.data.frame(summarytable(breastmilk.m1.s, breastmilk.m2.s, breastmilk.m3.s, breastmilk.m4.s, breastmilk.m5.s, breastmilk.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(breastmilk.m1.s, breastmilk.m2.s, breastmilk.m3.s, breastmilk.m4.s, breastmilk.m5.s, breastmilk.m6.s, which = c("AIC", "BIC", "entropy"))  
plot(breastmilk.m4.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE) 

## Number of participants/group?
breastmilk.wide <- Reduce(function(x, y) merge(x, y, all=TRUE), list(breastmilk.wide, breastmilk.m4.s$pprob))
table(breastmilk.wide$class)

## Plots
# Output predictions by class
breastmilk.pred <- data.frame(ageinweeks=seq(0,36))
breastmilk.fit <- predictY(breastmilk.m4.s,draws=TRUE,newdata=breastmilk.pred,var.time="ageinweeks")
breastmilk.fit.pred <- cbind(breastmilk.fit$time, breastmilk.fit$pred)
breastmilk.fit <- breastmilk.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
breastmilk.fit <- breastmilk.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
breastmilk.fit$norm_bmiz = breastmilk.fit$Ypred
breastmilk.fit$lower.Ypred <- ifelse(breastmilk.fit$lower.Ypred > -2.5, breastmilk.fit$lower.Ypred, NA)
breastmilk.fit$label <- ifelse(breastmilk.fit$class == 'class1', 'Decrease, n=22',
                         ifelse(breastmilk.fit$class =='class2', "Rapid Increase, n=10",
                                ifelse(breastmilk.fit$class =='class3', "Stable, n=49", "Delayed Increase, n=9")))

# Overlay with original data
breastmilk <- Reduce(function(x, y) merge(x, y, all=TRUE), list(breastmilk, breastmilk.m4.s$pprob))
breastmilk$class = str_c("Class ", breastmilk$class)
breastmilk$label <- ifelse(breastmilk$class == 'Class 1', "Decrease, n=22",
                               ifelse(breastmilk$class == 'Class 2', "Rapid Increase, n=10",
                                      ifelse(breastmilk$class == 'Class 3', "Stable, n=49", "Delayed Increase, n=9")))
# Plot
breastmilk.plot <- ggplot(breastmilk, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=breastmilk.fit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=breastmilk.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 


## cowsmilk ##
cowsmilk.m1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=cowsmilk)
cowsmilk.m2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=cowsmilk, ng=2, B=cowsmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=cowsmilk.m1.s)
cowsmilk.m3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=cowsmilk, ng=3, B=cowsmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=cowsmilk.m1.s)
cowsmilk.m4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=cowsmilk, ng=4, B=cowsmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=cowsmilk.m1.s)
cowsmilk.m5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=cowsmilk, ng=5, B=cowsmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=cowsmilk.m1.s)
cowsmilk.m6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=cowsmilk, ng=6, B=cowsmilk.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=cowsmilk.m1.s)
cowsmilk.tab.s<-as.data.frame(summarytable(cowsmilk.m1.s, cowsmilk.m2.s, cowsmilk.m3.s, cowsmilk.m4.s, cowsmilk.m5.s, cowsmilk.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(cowsmilk.m1.s, cowsmilk.m2.s, cowsmilk.m3.s, cowsmilk.m4.s, cowsmilk.m5.s, cowsmilk.m6.s, which = c("AIC", "BIC", "entropy"))  
plot(cowsmilk.m4.s, which="fit", var.time="ageinweeks", marg=TRUE, shades = TRUE) 
## NOTE:A 3-group solution may be slightly better but not enough for me to abandon the consistency from the 4-group solution

cowsmilk <- subset(long2, long2$feedinggroup==2)
cowsmilk.wide <- subset(wide, wide$feedinggroup==2)

## Number of participants/group?
cowsmilk.wide <- Reduce(function(x, y) merge(x, y, all=TRUE), list(cowsmilk.wide, cowsmilk.m4.s$pprob))
table(cowsmilk.wide$class)

## Plots
# Output predictions by class
cowsmilk.pred <- data.frame(ageinweeks=seq(0,36))
cowsmilk.fit <- predictY(cowsmilk.m4.s,draws=TRUE,newdata=cowsmilk.pred,var.time="ageinweeks")
cowsmilk.fit.pred <- cbind(cowsmilk.fit$time, cowsmilk.fit$pred)
cowsmilk.fit <- cowsmilk.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
cowsmilk.fit <- cowsmilk.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
cowsmilk.fit$norm_bmiz = cowsmilk.fit$Ypred
cowsmilk.fit$label <- ifelse(cowsmilk.fit$class == 'class1', "Decrease, n=21",
                             ifelse(cowsmilk.fit$class =='class2', "Rapid Increase, n=19",
                                    ifelse(cowsmilk.fit$class =='class3', "Stable, n=45", "Delayed Increase, n=24")))

# Overlay with original data
cowsmilk <- Reduce(function(x, y) merge(x, y, all=TRUE), list(cowsmilk, cowsmilk.m4.s$pprob))
cowsmilk$class = str_c("Class ", cowsmilk$class)
cowsmilk$label <- ifelse(cowsmilk$class == 'Class 1', "Decrease, n=21",
                             ifelse(cowsmilk$class =='Class 2', "Rapid Increase, n=19",
                                    ifelse(cowsmilk$class =='Class 3', "Stable, n=45", "Delayed Increase, n=24")))
# Plot
cowsmilk.plot <- ggplot(cowsmilk, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=cowsmilk.fit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=cowsmilk.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 



## Output table of model fit statistics
breastmilk.tab.s$data <- "breastmilk"
cowsmilk.tab.s$data <- "cowsmilk"

tab.s <- rbind(breastmilk.tab.s,cowsmilk.tab.s)
write.csv(tab.s, "tabforGMM_diet.csv")


## Merge all bmizplots together
comboplot = grid.arrange(female.plot + ggtitle("A) Females"),
                         male.plot + ggtitle("B) Males"),
                         breastmilk.plot + ggtitle("C) Breastmilk diet"),
                         cowsmilk.plot + ggtitle("D) Cows-milk diet"),
                         ncol=2, widths=c(4,4))

ggsave(plot = comboplot, filename = "subgroups_bmizplot_4groupings.jpg", width = 16, height = 12)


comboplot = grid.arrange(female.plot + ggtitle("A) Females"),
                         male.plot + ggtitle("B) Males"),
                         ncol=2, widths=c(4,4))

ggsave(plot = comboplot, filename = "bmizplot_sex.jpg", width = 16, height = 6)


comboplot = grid.arrange(female.plot + ggtitle("A) Females all ages"),
                         female28.plot + ggtitle("B) Females ≤28 weeks"),
                         male.plot + ggtitle("C) Males all ages"),
                         ncol=3, widths=c(4,4,4))
ggsave(plot = comboplot, filename = "bmizplot_sex_sens.jpg", width = 22, height = 6)


comboplot = grid.arrange(breastmilk.plot + ggtitle("A) Breastmilk diet"),
                         cowsmilk.plot + ggtitle("B) Cows-milk diet"),
                         ncol=2, widths=c(4,4))

ggsave(plot = comboplot, filename = "bmizplot_diet.jpg", width = 16, height = 6)




#########################
## Run for BMI z-score ##
## In overall cohort   ##
#########################

## Run gaussian model for Centered bmi z-score with ageinweeks as the time variable.
# We allow for random intercepts and random effects and run stratified by infant sex
# We initialize using the gridsearch function based on a 1-group model

set.seed(333)

# Cubic B-splines #
overallm1.cs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=long2)
overallm1.cs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=3) , random =~1, subject = "ifedid", data=long2)
overallm1.cs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=3) , random =~1, subject = "ifedid", data=long2)

overallm1.cs.4$pred$resid_ss_sq <- overallm1.cs.4$pred$resid_ss^2
overallm1.cs.4.sse <- sum(overallm1.cs.4$pred$resid_ss_sq) 

overallm1.cs.3$pred$resid_ss_sq <- overallm1.cs.3$pred$resid_ss^2
overallm1.cs.3.sse <- sum(overallm1.cs.3$pred$resid_ss_sq) 

overallm1.cs.2$pred$resid_ss_sq <- overallm1.cs.2$pred$resid_ss^2
overallm1.cs.2.sse <- sum(overallm1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
overallm1.qs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=2) , random =~1, subject = "ifedid", data=long2)
overallm1.qs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=2) , random =~1, subject = "ifedid", data=long2)
overallm1.qs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=2) , random =~1, subject = "ifedid", data=long2)

overallm1.qs.4$pred$resid_ss_sq <- overallm1.qs.4$pred$resid_ss^2
overallm1.qs.4.sse <- sum(overallm1.qs.4$pred$resid_ss_sq) 

overallm1.qs.3$pred$resid_ss_sq <- overallm1.qs.3$pred$resid_ss^2
overallm1.qs.3.sse <- sum(overallm1.qs.3$pred$resid_ss_sq) 

overallm1.qs.2$pred$resid_ss_sq <- overallm1.qs.2$pred$resid_ss^2
overallm1.qs.2.sse <- sum(overallm1.qs.2$pred$resid_ss_sq) 

# Compare SSEs #
overallm1.cs.4.sse
overallm1.cs.3.sse
overallm1.cs.2.sse

overallm1.qs.4.sse
overallm1.qs.3.sse
overallm1.qs.2.sse
## Lowest SSE is from: cubic model with 4 knots##

## Determine Number of Groups ##
overallm1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=long2)
overallm2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=2, B=overallm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=overallm1.s)
overallm3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=3, B=overallm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=overallm1.s)
overallm4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=4, B=overallm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=overallm1.s)
overallm5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=5, B=overallm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=overallm1.s)
overallm6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=6, B=overallm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=overallm1.s)
overalltab.s<-as.data.frame(summarytable(overallm1.s, overallm2.s, overallm3.s, overallm4.s, overallm5.s, overallm6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(overallm1.s, overallm2.s, overallm3.s, overallm4.s, overallm5.s, overallm6.s, which = c("AIC", "BIC", "entropy"))  
plot(overallm4.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE) 
write.csv(overalltab.s, "tabforGMM_overall.csv")

## Number of participants/group?
wide2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(wide, overallm4.s$pprob))
table(wide2$class)

## Plots
# Output predictions by class
overallpred <- data.frame(ageinweeks=seq(0,36))
overallfit <- predictY(overallm4.s,draws=TRUE,newdata=overallpred,var.time="ageinweeks")
overallfit.pred <- cbind(overallfit$time, overallfit$pred)
overallfit <- overallfit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
overallfit <- overallfit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
overallfit$norm_bmiz = overallfit$Ypred
overallfit$lower.Ypred <- ifelse(overallfit$lower.Ypred > -2.5, overallfit$lower.Ypred, NA) ## This was added due to weird values for the lower CI in class 4
overallfit$label <- ifelse(overallfit$class == 'class1', "Stable, n=82",
                               ifelse(overallfit$class =='class2', "Rapid Increase, n=47",
                                      ifelse(overallfit$class =='class3', "Delayed Increase, n=27","Decrease, n=43")))

# Overlay with original data
long3 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(long2, overallm4.s$pprob))
long3$class = str_c("Class ", long3$class)
long3$label <- ifelse(long3$class == 'Class 1', "Stable, n=82",
                          ifelse(long3$class =='Class 2', "Rapid Increase, n=47",
                                 ifelse(long3$class =='Class 3', "Delayed Increase, n=27","Decrease, n=43")))

# Plot
overall.plot <- ggplot(long3, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=overallfit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=overallfit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() + ggtitle("A) ", subtitle = " ") +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        legend.position = "left", legend.text=element_text(size=14), legend.title = element_text(size=16)) 

ggsave(plot = overall.plot, filename = "bmizplot_overall.jpg", width = 8, height = 6)



###################################################
## Run for BMI z-score                           ##
## In overall cohort                             ##
## But limit to visits after exposure collection ##
## at the last study visit (12+ weeks)           ##
###################################################
table(long2$visit) 
limit <- subset(long2, visit >= 12) ## 1065 observations
length(table(limit$ifedid)) ## n=192
quantile(table(limit$ifedid)) ## Median (IQR) number of visits: 5 (5, 7)

quantile(limit$ageinweeks,c(.20,.40,.60,.80))

## Test splines models: b-splines ##
## 2 knots: 16.857, 24.857 weeks split the age distribution at 33 and 66th percentiles ##
## 3 knots: 15.857, 21, 27.857 weeks split the age distribution at 25, 50, 75th percentiles ##
## 4 knots: 15.11429 19.51429 24.00000 28.42857  weeks split the age distribution at 20, 40, 60, and 80th percentiles ##

## Run gaussian model for Centered bmi z-score with ageinweeks as the time variable.
# We allow for random intercepts and random effects and run stratified by infant sex
# We initialize using the gridsearch function based on a 1-group model
set.seed(333)

# Cubic B-splines #
limit.m1.cs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3) , random =~1, subject = "ifedid", data=limit)
limit.m1.cs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.857, 21, 27.857), degree=3) , random =~1, subject = "ifedid", data=limit)
limit.m1.cs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(16.857, 24.857), degree=3) , random =~1, subject = "ifedid", data=limit)

limit.m1.cs.4$pred$resid_ss_sq <- limit.m1.cs.4$pred$resid_ss^2
limit.m1.cs.4.sse <- sum(limit.m1.cs.4$pred$resid_ss_sq) 

limit.m1.cs.3$pred$resid_ss_sq <- limit.m1.cs.3$pred$resid_ss^2
limit.m1.cs.3.sse <- sum(limit.m1.cs.3$pred$resid_ss_sq) 

limit.m1.cs.2$pred$resid_ss_sq <- limit.m1.cs.2$pred$resid_ss^2
limit.m1.cs.2.sse <- sum(limit.m1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
limit.m1.qs.4 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=2) , random =~1, subject = "ifedid", data=limit)
limit.m1.qs.3 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.857, 21, 27.857), degree=2) , random =~1, subject = "ifedid", data=limit)
limit.m1.qs.2 <- hlme(norm_bmiz ~ bs(ageinweeks, knots=c(16.857, 24.857), degree=2) , random =~1, subject = "ifedid", data=limit)

limit.m1.qs.4$pred$resid_ss_sq <- limit.m1.qs.4$pred$resid_ss^2
limit.m1.qs.4.sse <- sum(limit.m1.qs.4$pred$resid_ss_sq) 

limit.m1.qs.3$pred$resid_ss_sq <- limit.m1.qs.3$pred$resid_ss^2
limit.m1.qs.3.sse <- sum(limit.m1.qs.3$pred$resid_ss_sq) 

limit.m1.qs.2$pred$resid_ss_sq <- limit.m1.qs.2$pred$resid_ss^2
limit.m1.qs.2.sse <- sum(limit.m1.qs.2$pred$resid_ss_sq) 

# Compare SSEs #
limit.m1.cs.4.sse
limit.m1.cs.3.sse
limit.m1.cs.2.sse

limit.m1.qs.4.sse
limit.m1.qs.3.sse
limit.m1.qs.2.sse
## Lowest SSE is from: cubic model with 4 knots##


## Determine Number of Groups ##
limit.m1.s<-hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3) , random =~1, subject = "ifedid", data=limit)
limit.m2.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3), random =~1, subject = "ifedid", data=limit, ng=2, B=limit.m1.s, mixture=~bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3)), rep=20, maxiter=20, minit=limit.m1.s)
limit.m3.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3), random =~1, subject = "ifedid", data=limit, ng=3, B=limit.m1.s, mixture=~bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3)), rep=20, maxiter=20, minit=limit.m1.s)
limit.m4.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3), random =~1, subject = "ifedid", data=limit, ng=4, B=limit.m1.s, mixture=~bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3)), rep=20, maxiter=20, minit=limit.m1.s)
limit.m5.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3), random =~1, subject = "ifedid", data=limit, ng=5, B=limit.m1.s, mixture=~bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3)), rep=20, maxiter=20, minit=limit.m1.s)
limit.m6.s<-gridsearch(hlme(norm_bmiz ~ bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3), random =~1, subject = "ifedid", data=limit, ng=6, B=limit.m1.s, mixture=~bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3)), rep=20, maxiter=20, minit=limit.m1.s)
limit.tab.s<-as.data.frame(summarytable(limit.m1.s, limit.m2.s, limit.m3.s, limit.m4.s, limit.m5.s, limit.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(limit.m1.s, limit.m2.s, limit.m3.s, limit.m4.s, limit.m5.s, limit.m6.s, which = c("AIC", "BIC", "entropy"))  
## None of these have a good entropy, and only the 2-group solution has >10% within one group
plot(limit.m2.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE) 
plot(limit.m4.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE)  #n=7 in Rapid Increase 
write.csv(limit.tab.s, "tabforGMM_limit.csv")


### 2-group solution
## Number of participants/group?
wide2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(wide, limit.m2.s$pprob))
table(wide2$class)

## Plots
# Output predictions by class
limit.pred <- data.frame(ageinweeks=seq(12,36))
limit.fit <- predictY(limit.m2.s,draws=TRUE,newdata=limit.pred,var.time="ageinweeks")
limit.fit.pred <- cbind(limit.fit$time, limit.fit$pred)
limit.fit <- limit.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
limit.fit <- limit.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
limit.fit$norm_bmiz = limit.fit$Ypred
limit.fit$lower.Ypred <- ifelse(limit.fit$lower.Ypred > -2.5, limit.fit$lower.Ypred, NA) 
limit.fit$label <- ifelse(limit.fit$class == 'class1', "Delayed Increase, n=92", "Stable, n=100")

# Overlay with original data
long3 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(limit, limit.m2.s$pprob))
long3$class = str_c("Class ", long3$class)
long3$label <- ifelse(long3$class == 'Class 1', "Delayed Increase, n=92", "Stable, n=100")

# Plot
limit.plot <- ggplot(long3, aes(x=ageinweeks,y=norm_bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=limit.fit, aes(x=ageinweeks, y=norm_bmiz, col=as.factor(label))) + 
  geom_ribbon(data=limit.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered BMI z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 

ggsave(plot = limit.plot, filename = "bmizplot_limit2.jpg", width = 8, height = 5)



########################################
## Run for BMI z-score NOT Centered ##
## In overall cohort                  ##
########################################

## Run gaussian model for Centered bmi z-score with ageinweeks as the time variable.
# We allow for random intercepts and random effects and run stratified by infant sex
# We initialize using the gridsearch function based on a 1-group model

set.seed(333)

# Cubic B-splines #
long2.m1.cs.4 <- hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=long2)
long2.m1.cs.3 <- hlme(bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=3) , random =~1, subject = "ifedid", data=long2)
long2.m1.cs.2 <- hlme(bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=3) , random =~1, subject = "ifedid", data=long2)

long2.m1.cs.4$pred$resid_ss_sq <- long2.m1.cs.4$pred$resid_ss^2
long2.m1.cs.4.sse <- sum(long2.m1.cs.4$pred$resid_ss_sq) 

long2.m1.cs.3$pred$resid_ss_sq <- long2.m1.cs.3$pred$resid_ss^2
long2.m1.cs.3.sse <- sum(long2.m1.cs.3$pred$resid_ss_sq) 

long2.m1.cs.2$pred$resid_ss_sq <- long2.m1.cs.2$pred$resid_ss^2
long2.m1.cs.2.sse <- sum(long2.m1.cs.2$pred$resid_ss_sq) 

# Quadratic B-splines #
long2.m1.qs.4 <- hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=2) , random =~1, subject = "ifedid", data=long2)
long2.m1.qs.3 <- hlme(bmiz ~ bs(ageinweeks, knots=c(4, 11.4286, 23), degree=2) , random =~1, subject = "ifedid", data=long2)
long2.m1.qs.2 <- hlme(bmiz ~ bs(ageinweeks, knots=c(5.85714, 19), degree=2) , random =~1, subject = "ifedid", data=long2)

long2.m1.qs.4$pred$resid_ss_sq <- long2.m1.qs.4$pred$resid_ss^2
long2.m1.qs.4.sse <- sum(long2.m1.qs.4$pred$resid_ss_sq) 

long2.m1.qs.3$pred$resid_ss_sq <- long2.m1.qs.3$pred$resid_ss^2
long2.m1.qs.3.sse <- sum(long2.m1.qs.3$pred$resid_ss_sq) 

long2.m1.qs.2$pred$resid_ss_sq <- long2.m1.qs.2$pred$resid_ss^2
long2.m1.qs.2.sse <- sum(long2.m1.qs.2$pred$resid_ss_sq) 

# Compare SSEs #
long2.m1.cs.4.sse
long2.m1.cs.3.sse
long2.m1.cs.2.sse

long2.m1.qs.4.sse
long2.m1.qs.3.sse
long2.m1.qs.2.sse
## Lowest SSE is from: cubic model with 4 knots##

## Determine Number of Groups ##
long2.m1.s<-hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=long2)
long2.m2.s<-gridsearch(hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=2, B=long2.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=long2.m1.s)
long2.m3.s<-gridsearch(hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=3, B=long2.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=long2.m1.s)
long2.m4.s<-gridsearch(hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=4, B=long2.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=long2.m1.s)
long2.m5.s<-gridsearch(hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=5, B=long2.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=long2.m1.s)
long2.m6.s<-gridsearch(hlme(bmiz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=6, B=long2.m1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=long2.m1.s)
long2.tab.s<-as.data.frame(summarytable(long2.m1.s, long2.m2.s, long2.m3.s, long2.m4.s, long2.m5.s, long2.m6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(long2.m1.s, long2.m2.s, long2.m3.s, long2.m4.s, long2.m5.s, long2.m6.s, which = c("AIC", "BIC", "entropy"))  
plot(long2.m4.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE) 
write.csv(long2.tab.s, "tabforGMM_overall_nonorm.csv")

## Number of participants/group?
wide2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(wide, long2.m4.s$pprob))
table(wide2$class)

## Plots
# Output predictions by class
long2.pred <- data.frame(ageinweeks=seq(0,36))
long2.fit <- predictY(long2.m4.s,draws=TRUE,newdata=long2.pred,var.time="ageinweeks")
long2.fit.pred <- cbind(long2.fit$time, long2.fit$pred)
long2.fit <- long2.fit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
long2.fit <- long2.fit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
long2.fit$bmiz = long2.fit$Ypred
long2.fit$lower.Ypred <- ifelse(long2.fit$lower.Ypred > -2.5, long2.fit$lower.Ypred, NA) 
long2.fit$label <- ifelse(long2.fit$class == 'class1', "Delayed Increase, n=25",
                          ifelse(long2.fit$class =='class2', "Stable, n=80",
                                 ifelse(long2.fit$class =='class3', "Decrease, n=49", "Rapid Increase, n=45")))

# Overlay with original data
long3 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(long2, long2.m4.s$pprob))
long3$class = str_c("Class ", long3$class)
long3$label <- ifelse(long3$class == 'Class 1', "Delayed Increase, n=25",
                          ifelse(long3$class =='Class 2', "Stable, n=80",
                                 ifelse(long3$class =='Class 3', "Decrease, n=49", "Rapid Increase, n=45")))

# Plot
long2.plot <- ggplot(long3, aes(x=ageinweeks,y=bmiz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=long2.fit, aes(x=ageinweeks, y=bmiz, col=as.factor(label))) + 
  geom_ribbon(data=long2.fit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "BMI z-score") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 

ggsave(plot = long2.plot, filename = "bmizplot_overall_nonorm.jpg", width = 8, height = 5)



############################
## Run for weight z-score ##
## In overall cohort      ##
############################

## Run gaussian model for Centered bmi z-score with ageinweeks as the time variable.
# We allow for random intercepts and random effects and run stratified by infant sex
# We initialize using the gridsearch function based on a 1-group model

set.seed(333)

## Determine Number of Groups ##
norm_wazm1.s<-hlme(norm_waz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) , random =~1, subject = "ifedid", data=long2)
norm_wazm2.s<-gridsearch(hlme(norm_waz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=2, B=norm_wazm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=norm_wazm1.s)
norm_wazm3.s<-gridsearch(hlme(norm_waz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=3, B=norm_wazm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=norm_wazm1.s)
norm_wazm4.s<-gridsearch(hlme(norm_waz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=4, B=norm_wazm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=norm_wazm1.s)
norm_wazm5.s<-gridsearch(hlme(norm_waz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=5, B=norm_wazm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=norm_wazm1.s)
norm_wazm6.s<-gridsearch(hlme(norm_waz ~ bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3), random =~1, subject = "ifedid", data=long2, ng=6, B=norm_wazm1.s, mixture=~bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3)), rep=20, maxiter=20, minit=norm_wazm1.s)
norm_waztab.s<-as.data.frame(summarytable(norm_wazm1.s, norm_wazm2.s, norm_wazm3.s, norm_wazm4.s, norm_wazm5.s, norm_wazm6.s, which = c('G', "AIC", "BIC", "entropy", "%class")))
summaryplot(norm_wazm1.s, norm_wazm2.s, norm_wazm3.s, norm_wazm4.s, norm_wazm5.s, norm_wazm6.s, which = c("AIC", "BIC", "entropy"))  
plot(norm_wazm4.s, which="fit", var.time="ageinweeks", marg=FALSE, shades = TRUE) 
write.csv(norm_waztab.s, "tabforGMM_norm_waz.csv")


### Four-group solution
## Number of participants/group?
wide2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(wide, norm_wazm4.s$pprob))
table(wide2$class)

## Plots
# Output predictions by class
norm_wazpred <- data.frame(ageinweeks=seq(0,36))
norm_wazfit <- predictY(norm_wazm4.s,draws=TRUE,newdata=norm_wazpred,var.time="ageinweeks")
norm_wazfit.pred <- cbind(norm_wazfit$time, norm_wazfit$pred)
norm_wazfit <- norm_wazfit.pred %>% pivot_longer(
  cols = -ageinweeks,
  names_to = c('name1', 'class'),
  names_pattern = "(.+)_(.+$)")
norm_wazfit <- norm_wazfit %>% pivot_wider(
  id_cols = c('ageinweeks','class'),
  names_from = name1,
  values_from = value
)
norm_wazfit$norm_waz = norm_wazfit$Ypred
norm_wazfit$lower.Ypred <- ifelse(norm_wazfit$lower.Ypred > -2.5, norm_wazfit$lower.Ypred, NA)
norm_wazfit$label <- ifelse(norm_wazfit$class == 'class1', "Rapid Increase, n=30",
                          ifelse(norm_wazfit$class =='class2', "Decrease, n=66",
                                 ifelse(norm_wazfit$class =='class3', "Stable, n=88", "Delayed Increase, n=15")))

# Overlay with original data
long3 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(long2, norm_wazm4.s$pprob))
long3$class = str_c("Class ", long3$class)
long3$label <- ifelse(long3$class == 'Class 1', "Rapid Increase, n=30",
                      ifelse(long3$class =='Class 2', "Decrease, n=66",
                             ifelse(long3$class =='Class 3', "Stable, n=88", "Delayed Increase, n=15")))

# Plot
norm_wazplot <- ggplot(long3, aes(x=ageinweeks,y=norm_waz,col=as.factor(label),shape=as.factor(label),linetype=as.factor(label))) + geom_point(alpha=0.2,size=2) +
  geom_line(data=norm_wazfit, aes(x=ageinweeks, y=norm_waz, col=as.factor(label))) + 
  geom_ribbon(data=norm_wazfit, aes(x=ageinweeks, ymin=lower.Ypred, ymax=upper.Ypred, col=as.factor(label), fill=as.factor(label)), linetype=0, alpha=0.2) +
  labs(x = "Age (weeks)", y = "Centered weight z-score") + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), limits=c(-3,3)) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette = "Dark2", 
                    guide = guide_legend(title="",reverse=TRUE)) +
  scale_shape_manual(values = c(15,16,17,18),
                     guide = guide_legend(title="",reverse=TRUE)) +
  scale_linetype_manual(values = c("solid",'twodash','longdash','dotted'),
                        guide = guide_legend(title="",reverse=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) 


ggsave(plot = norm_wazplot, filename = "norm_wazplot_overall.jpg", width = 8, height = 5)


comboplot = grid.arrange(long2.plot + ggtitle("A) Non-Centered BMI z-score"),
                         norm_wazplot + ggtitle("B) Centered Weight-for-Age Z-score"),
                         ncol=2, widths=c(4,4))

ggsave(plot = comboplot, filename = "bmiz_waz_plot.jpg", width = 16, height = 6)



#####################################################################
## Save posterior probabilities and class assignments into dataset ##
#####################################################################
malepprob <- male.m4.s$pprob
colnames(malepprob)[2:6] <- paste('male', colnames(malepprob)[2:6], sep='.')
femalepprob <- female.m4.s$pprob
colnames(femalepprob)[2:6] <- paste('female', colnames(femalepprob)[2:6], sep='.')
female28pprob <- female28.m4.s$pprob
colnames(female28pprob)[2:6] <- paste('female28', colnames(female28pprob)[2:6], sep='.')
breastmilkpprob <- breastmilk.m4.s$pprob
colnames(breastmilkpprob)[2:6] <- paste('breastmilk', colnames(breastmilkpprob)[2:6], sep='.')
cowsmilkpprob4 <- cowsmilk.m4.s$pprob
colnames(cowsmilkpprob4)[2:6] <- paste('cowsmilk4', colnames(cowsmilkpprob4)[2:6], sep='.')
overallpprob <- overallm4.s$pprob
colnames(overallpprob)[2:6] <- paste('overall', colnames(overallpprob)[2:6], sep='.')
# limitpprob <- limit.m2.s$pprob
# colnames(limitpprob)[2:4] <- paste('limit', colnames(limitpprob)[2:4], sep='.')
nonormpprob <- long2.m4.s$pprob
colnames(nonormpprob)[2:6] <- paste('nonorm', colnames(nonormpprob)[2:6], sep='.')
wazpprob <- norm_wazm4.s$pprob
colnames(wazpprob)[2:6] <- paste('waz', colnames(wazpprob)[2:6], sep='.')


outcomes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(malepprob, femalepprob, female28pprob, breastmilkpprob, cowsmilkpprob4,
                                                              overallpprob, nonormpprob, wazpprob)) 

write.csv(outcomes, "J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\GMM_Dataset_Wide.csv")




################################################################################
## Additional plot which shows the BMI z-scores over time by group assignment ##
################################################################################
outcomes <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\GMM_Dataset_Wide.csv")

dat.w.outcomes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(long, outcomes))
dat.bmiz <- dat.w.outcomes[, c('ifedid','visit','ageinweeks','bmiz',"n_bmiz","overall.class")]
dat.bmiz <- dat.bmiz[complete.cases(dat.bmiz),]

## Average number of measures of BMI by group?
data_table <- data.table(dat.bmiz)  
n_bmiz <- data_table[ , .(count = length(unique(bmiz))), by = ifedid]
dat.bmiz <- Reduce(function(x, y) merge(x, y, all=TRUE), list(dat.bmiz, n_bmiz))
dat.bmiz$count # By ID

vis1 <- subset(dat.bmiz, visit == 1)
CreateTableOne(data=vis1,
               vars='count',
               factorVars='overall.class',
               strata="overall.class")  # By group - all have about 10 measures


## 4 panel spaghetti plot by group with counts included
# Group 1
g1 <- subset(dat.bmiz, overall.class == 1)
g1.bmizplot <- ggplot(data = g1, aes(x = ageinweeks, y = bmiz, group = ifedid)) + 
  geom_line(color = "#E7298A") + xlab("Age (weeks)\n ") + ylab(" \nBMI z-score") + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -1.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 1.96, linetype = 'dashed', color = 'grey') +
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) +
  ggtitle("Stable", subtitle = "82 infants with a mean of 9.8 visits")
g1.bmizplot

g2 <- subset(dat.bmiz, overall.class == 2)
g2.bmizplot <- ggplot(data = g2, aes(x = ageinweeks, y = bmiz, group = ifedid)) + 
  geom_line(color = "#7570B3") + xlab("Age (weeks)\n ") + ylab(" \nBMI z-score") + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -1.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 1.96, linetype = 'dashed', color = 'grey') +
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) +
  ggtitle("Rapid Increase", subtitle = "47 infants with a mean of 10.3 visits")
g2.bmizplot

g3 <- subset(dat.bmiz, overall.class == 3)
g3.bmizplot <- ggplot(data = g3, aes(x = ageinweeks, y = bmiz, group = ifedid)) + 
  geom_line(color = "#D95F02") + xlab("Age (weeks)\n ") + ylab(" \nBMI z-score") + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -1.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 1.96, linetype = 'dashed', color = 'grey') +
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) +
  ggtitle("Delayed Increase", subtitle = "27 infants with a mean of 9.4 visits")
g3.bmizplot

g4 <- subset(dat.bmiz, overall.class == 4)
g4.bmizplot <- ggplot(data = g4, aes(x = ageinweeks, y = bmiz, group = ifedid)) + 
  geom_line(color = "#1B9E77") + xlab("Age (weeks)\n ") + ylab(" \nBMI z-score") + 
  scale_x_continuous(breaks=c(0,2,4,6,8,12,16,20,24,28,32,36)) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), limits=c(-4,4)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = -1.96, linetype = 'dashed', color = 'grey') +
  geom_hline(yintercept = 1.96, linetype = 'dashed', color = 'grey') +
  theme_bw() +  
  theme(text = element_text(size=19),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.text=element_text(size=14), legend.title = element_text(size=16)) +
  ggtitle("Decrease", subtitle = "43 infants with a mean of 10.4 visits")
g4.bmizplot


comboplot = grid.arrange(g1.bmizplot,g2.bmizplot,g3.bmizplot,g4.bmizplot, ncol=2, widths = c(4,4))
comboplot
ggsave(plot = comboplot, filename = "combo_bmizplot_bygroup.jpg", width = 16, height = 12, dpi=600)


overall.plot <- overall.plot + theme(title = element_blank())
ggsave(plot = overall.plot, filename = "overall_bmizplot.jpg", width = 8, height = 5, dpi=600)

overall.plot <- overall.plot + ggtitle("A) ") + theme(title = element_text(size = 16))
comboplot = grid.arrange(g1.bmizplot,g2.bmizplot,g3.bmizplot,g4.bmizplot, ncol=2, widths = c(4,4), top=textGrob("B) ", gp = gpar(fontsize = 19), hjust = 0, x = 0))
comboplot

layout <- "
AA
##
BB
BB"
comboplot2 <- overall.plot  / comboplot  + 
  plot_layout(design = layout, heights = c(5, 2, 10))

comboplot2
ggsave(plot = comboplot2, filename = "Combined BMIZ Trajectory Plots.jpg", width = 16, height = 20, dpi=900)






## Difference in BMI z at visit 1 by group?
dat.w.outcomes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(wide, outcomes))
dat.w.outcomes$bmiz_1
dat.bmiz <- dat.w.outcomes[, c('ifedid','bmiz_1',"overall.class")]
dat.bmiz <- dat.bmiz[complete.cases(dat.bmiz),]
tab1 <- CreateTableOne(data=dat.bmiz,
               vars='bmiz_1',
               factorVars='overall.class',
               strata="overall.class")
tab1 <- print(tab1)
write.csv(tab1, "J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Initial BMIz by Group Assignment.csv")
#                     1 (Stable)    2 (Rapid)      3 (Delayed)    4 (Decrease)
# bmiz_1 (mean (SD)) -0.18 (0.70)   -0.57 (0.82)   -0.06 (0.82)   0.39 (0.79) <0.001 

## Comparisons of stable to rapid
dat.bmiz.limit <- subset(dat.bmiz, overall.class == 1 | overall.class == 2)
CreateTableOne(data=dat.bmiz.limit,
               vars='bmiz_1',
               factorVars='overall.class',
               strata="overall.class")
## This is different as well



## Difference in BMI z at last visit by group?
dat.w.outcomes <- Reduce(function(x, y) merge(x, y, all=TRUE), list(long, outcomes))
dat.bmiz <- dat.w.outcomes[, c('ifedid','visit','ageinweeks','bmiz',"n_bmiz","overall.class")]
dat.bmiz <- dat.bmiz[complete.cases(dat.bmiz),]
data_table <- data.table(dat.bmiz)  
dat.last.bmiz <- data_table[, .SD[which.max(ageinweeks), ], by=ifedid] 
tab1 <- CreateTableOne(data=dat.last.bmiz,
                       vars='bmiz',
                       factorVars='overall.class',
                       strata="overall.class")
tab1 <- print(tab1)
write.csv(tab1, "J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results\\Last BMIz by Group Assignment.csv")
#                     1 (Stable)    2 (Rapid)      3 (Delayed)    4 (Decrease)
#  n                  82            47             27             43       
# bmiz (mean (SD))    0.17 (0.80)   1.37 (0.85)    0.15 (1.01)    -0.56 (0.84)
# p-value of difference in bmiz (mean (SD)) <0.001   
