# Programmer: Danielle Stevens
# Last updated on: 11/17/2026
# 
# Purpose of code: This code reads in the imputed IFED exposure, covariate, and confounder data,
#                  Runs exposure distribution tables and figures,
#                  And investigates associations between phenols and growth outcomes (long and wide form) in infancy


# install.packages('tidyverse') 
# install.packages('dplyr')
# install.packages('htmlTable') 
# install.packages('haven')
# install.packages('foreign')
# install.packages('mde')
# install.packages('mice') 
# install.packages('mitools')
# install.packages('miceadds')
# install.packages('broom.mixed')
# install.packages('nnet')
# install.packages('splines')
# install.packages('tableone')
# install.packages('pheatmap')
# install.packages('dichromat')
# install.packages('qgcomp')
# install.packages('qgcompint')
# install.packages('ggplot2')
# install.packages('lme4')
# install.packages('nlme')
# install.packages('irr')


#### Libraries
library(tidyverse) 
library(dplyr)
library(htmlTable) 
library(haven)
library(foreign)
library(mde)
library(mice) 
library(mitools)
library(miceadds)
library(broom.mixed)
library(nnet)
library(splines)
library(tableone)
library(pheatmap)
library(dichromat)
library(qgcomp)
library(qgcompint)
library(ggplot2)
library(lme4)
library(nlme)


##################
# Variable lists #
##################

#set working directory
setwd("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Results")


## Define variables of interest
# Exposures
exposure.labels <- c('Bisphenol A', 'Benzophenone-3', '2,4 Dichlorophenol', '2,5 Dichlorophenol', 
                     'Methyl Paraben', 'Propyl Paraben', 'Triclosan')
exposures <- c('bpa', 'bp_3', 'dcp_24', 'dcp_25', 'm_pb', 'p_pb', 'tcs')
exposures.obrien.v68 <- paste0(exposures,"_obrien_v68")
exposures.obrien.v12 <- paste0(exposures,"_obrien_v12")
exposures.obrien <- paste0(exposures,"_obrien_gm")
exposures.obrien.ln <- paste0(exposures,"_obrien_gm_ln")
exposures.gm.ln <- paste0(exposures,"_gm_ln")

# Confounders
confounders.long <- c('ageinweeks', 'feedinggroup', 'mage', 'meduc', 'childsgender', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')
confounders.wide <- c('feedinggroup', 'mage', 'meduc', 'childsgender', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')

###############################################################################################################################################################################################################
## Read in datasets
imp.dat.wide <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\Imputed dataset with outcomes wide.csv") ## THIS IS TOO LONG FOR A WIDEFORM DATAFRAME ##

## Flip to wide for weights analysis
imp.dat.weighted = imp.dat.wide %>% 
  pivot_longer(cols = contains('prob'),
               names_to = c(".value","group."),
               names_sep = "\\.")


################################################
## Table 2: Table of exposure characteristics ##
################################################

## I am just going to take the mean of the medians and the mean of the first and third quantiles for these tables based on:
## https://stackoverflow.com/questions/65379516/how-to-proceed-with-descriptive-statistics-median-iqr-frequencies-proportion 

# Get dataset

med.iqr <- data.frame()
funct.exp.dist <- function(inputvars,vartype){
  nomiss <- subset(imp.dat.wide[,c(".imp", inputvars)], .imp!=0)
  
  exposure.med <- aggregate(. ~ .imp, nomiss, median)
  mean.exposure.med <- data.frame(colMeans(exposure.med))
  exposure.q1 <- aggregate(. ~ .imp, nomiss, FUN='quantile', probs=c(0.25))
  mean.exposure.q1 <- data.frame(colMeans(exposure.q1))
  exposure.q3 <- aggregate(. ~ .imp, nomiss, FUN='quantile', probs=c(0.75))
  mean.exposure.q3 <- data.frame(colMeans(exposure.q3))
  
  med.iqr <- data.frame(cbind(vartype,mean.exposure.med,mean.exposure.q1,mean.exposure.q3))
  colnames(med.iqr)<-c("Variable","Median","Q1","Q3")
  med.iqr$MedIQR1<- paste0(round(med.iqr$Median,1)," (",round(med.iqr$Q1,1),", ",round(med.iqr$Q3,1),")")
  med.iqr$MedIQR2<- paste0(round(med.iqr$Median,2)," (",round(med.iqr$Q1,2),", ",round(med.iqr$Q3,2),")")
  med.iqr<-med.iqr[-1,]
  }
gm.dist <- funct.exp.dist(exposures.obrien,"GM")
v68.dist <- funct.exp.dist(exposures.obrien.v68,"V68")
v12.dist <- funct.exp.dist(exposures.obrien.v12,"V12")

exp.dist <- rbind(gm.dist,v68.dist,v12.dist)
exp.dist$Label <- rep(exposures,3)
exp.dist.wide = exp.dist %>% pivot_wider(id_cols = c("Label"),
                                         names_from = c("Variable"),
                                             values_from = c("MedIQR1"))

write.csv(exp.dist.wide, "Table 2 - Exposure with Obrien Correction Table.csv")




################################################################################
## Figure S2: Correlation matrix of exposures                                 ##
## We run into issues using the micombine.cor because we cannot use as.mids() ##
## after all the merging & variable editing we've done. So we're just going   ##
## to loop and average "naive correlation"                                    ##
################################################################################

#Running loop over imputed datasets
cormat = matrix(0, nrow = 7, ncol = 7)
for(j in 1:10) {
  #grab the jth imputed dataset
  imputation <- imp.dat.wide %>% filter(.imp==j)
  
  #create the design effect object for that dataset
  imputation.cor <- cor(imputation[,exposures.obrien], use = "complete.obs", method = "spearman")
  
  #add to the previous correlation matrix
  cormat = cormat + imputation.cor
}

#Averaging the output from the loop to create the pooled correlation matrix
cormat.avg = round(cormat / 10, digits = 2)
min(cormat.avg)
colnames(cormat.avg) <- exposure.labels
rownames(cormat.avg) <- exposure.labels
write.csv(cormat.avg,"Correlation Matrix.csv")

#plot
p <- pheatmap(cormat.avg,
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
              fontsize=25)


png(height=1800, width=1800, file="Correlation Matrix.png", type = "cairo")
p
dev.off()


########################################
## Predictors of exposure & outcomes  ##
## Just do ones Mandy did not include ##
########################################

imp = subset(imp.dat.wide, .imp ==0)
CreateTableOne(vars=exposures.obrien,
               factorVars='smkpreg',
               strata='smkpreg',
               data=imp) #Significant differences in TCS

CreateTableOne(vars=exposures.obrien,
               factorVars='alcpregany',
               strata='alcpregany',
               data=imp) #Significant differences in DCPs

CreateTableOne(vars=exposures.obrien,
               factorVars='gestationalage',
               strata='gestationalage',
               data=imp) #No significant differences but note that methyl paraben declines with increasing GA


##################################################################
## Did Gestational age at delievery impact growth trajectories? ##
##################################################################
imp$overall <- ifelse(imp$overall.class == '1', "Decrease",
                      ifelse(imp$overall.class =='2', "Delayed Increase",
                             ifelse(imp$overall.class =='3', "Rapid Growth","Stable Growth")))


CreateTableOne(vars='overall',
               factorVars='gestationalage',
               strata='gestationalage',
               data=imp) #No significant differences but note that methyl paraben declines with increasing GA



####################################################
## Table S2: Agreement between outcome categories ##
####################################################
imp$female <- ifelse(imp$female.class == '1', "Stable Growth",
                     ifelse(imp$female.class =='2', "Delayed Increase",
                            ifelse(imp$female.class =='3',"Rapid Growth", "Decrease")))
imp$female28 <- ifelse(imp$female28.class == '1', "Rapid Growth",
                       ifelse(imp$female28.class =='2', "Delayed Increase",
                              ifelse(imp$female28.class =='3', "Decrease", "Stable Growth")))
imp$male <- ifelse(imp$male.class == '1', "Stable Growth",
                   ifelse(imp$male.class =='2', "Delayed Increase",
                          ifelse(imp$male.class =='3', "Decrease","Rapid Growth")))
imp$breastmilk <- ifelse(imp$breastmilk.class == '1', "Decrease",
                         ifelse(imp$breastmilk.class =='2', "Rapid Growth",
                                ifelse(imp$breastmilk.class =='3', "Stable Growth","Delayed Increase")))
imp$cowsmilk <- ifelse(imp$cowsmilk4.class == '1', "Stable Growth",
                       ifelse(imp$cowsmilk4.class =='2', "Decrease",
                              ifelse(imp$cowsmilk4.class == "3", "Delayed Increase", "Rapid Growth")))
imp$nonorm <- ifelse(imp$nonorm.class == '1', "Rapid Growth",
                     ifelse(imp$nonorm.class =='2', "Stable Growth",
                            ifelse(imp$nonorm.class =='3', "Delayed Increase","Decrease")))
imp$waz <- ifelse(imp$waz.class == '1', "Stable Growth",
                  ifelse(imp$waz.class =='2', "Decrease",
                         ifelse(imp$waz.class =='3', "Rapid Growth", "Delayed Increase")))

library(Publish)
library(irr)
female.female28 <- kappa2(imp[,c('female','female28')])$value #.76
tabs1<-CreateTableOne(vars=c('female28'),
                      factorVars=c('female28','female'),
                      data=imp,
                      strata='female',
                      test=FALSE)
print.tabs1<- print(tabs1)
write.csv(print.tabs1,"Table of Overlap Between GMM Groupings for Female Trajectories.csv")

overall.female <- kappa2(imp[,c('overall','female')])$value #0.93
overall.female28 <- kappa2(imp[,c('overall','female28')])$value #0.78
overall.male <- kappa2(imp[,c('overall','male')])$value #0.82
overall.breastmilk <- kappa2(imp[,c('overall','breastmilk')])$value #0.82
overall.cowsmilk <- kappa2(imp[,c('overall','cowsmilk')])$value #0.75
overall.nonorm <- kappa2(imp[,c('overall','nonorm')])$value #0.87
overall.waz <- kappa2(imp[,c('overall','waz')])$value #0.59

tabs1<-CreateTableOne(vars=c('overall','female','female28','male','breastmilk','cowsmilk','nonorm','waz'),
                      factorVars=c('overall','female','female28','male','breastmilk','cowsmilk','nonorm','waz'),
                      data=imp,
                      strata='overall',
                      test=FALSE)
print.tabs1<- print(tabs1)
write.csv(print.tabs1,"Table of Overlap Between GMM Groupings.csv")



###########################################
## Calculate IQR of exposures for models ##
###########################################

## IQR-estimates for generating standardized results
iqr.all = data.frame()
for(i in 1:length(exposures)){
  term = exposures.obrien.ln[i]
  
  q = quantile(imp[,exposures.obrien.ln[i]], quantiles = c(0.25, 0.75), na.rm=T)
  iqr = round(q[[4]] - q[[2]],2)
  iqr = cbind(term, iqr)
  
  iqr.all = rbind(iqr.all, iqr)
}
iqr.all

lngm.dist <- funct.exp.dist(exposures.obrien.ln,"LN")
lngm.dist

########################
########################
## Multinomial Models ##
########################
########################


###########################################
## Primary analysis in overall dataframe ##
###########################################

## Edit dataframe for analysis of overall trajectories
imp.dat.weighted$class <- ifelse(imp.dat.weighted$overall.class == '1', "Decrease, n=43",
                           ifelse(imp.dat.weighted$overall.class =='2', "Delayed Increase, n=27",
                                  ifelse(imp.dat.weighted$overall.class =='3', "Rapid Growth, n=47","Stable Growth, n=82")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=82")
table(imp.dat.weighted$overall.class)
table(imp.dat.weighted$class)/11 #Check by dividing n's by 10 to get the n= for each group

## Create imputation list
list.df = list()
for(i in 1:10){
  list.df[[i]] = imp.dat.weighted[which(imp.dat.weighted$.imp==i),]
}
imp.list.df = imputationList(list.df)

## Test model
results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c('bpa_obrien_gm_ln', confounders.wide), collapse = "+"))), 
                                     weights = overall, na.omit=T))

tidied_list <- lapply(results, function(mod) {
  tidy(mod) %>%
    mutate(term = paste(y.level, term, sep = ": "))
})

pooled <- pool(results)
p.results = summary(pooled)
p.results$term2 <- tidied_list[[1]]$term
p.results <- p.results %>%
  separate(
    term2,
    into = c("y.level", "term"),
    sep = ": "
  )
p.results.term = p.results %>% filter(term == paste('bpa_obrien_gm_ln'))


#Unadjusted models
unadj.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", exposures.obrien.ln[i])), na.omit=T))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  unadj.results = rbind(unadj.results, p.results.term)
}

#Adjusted models
adj.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide), collapse = "+"))), 
                                       weights = overall, na.omit=T))

    tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )  

  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  adj.results = rbind(adj.results, p.results.term)
}

unadj.results2 = merge(unadj.results, iqr.all, by = "term", all = T)
unadj.results2 = unadj.results2 %>% mutate(iqr = as.numeric(iqr))
unadj.results2 = unadj.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
unadj.results2$OR = exp(unadj.results2$estimate*unadj.results2$iqr)
unadj.results2$LCL = exp(unadj.results2$estimate*unadj.results2$iqr - 1.96*(unadj.results2$std.error*unadj.results2$iqr))
unadj.results2$UCL = exp(unadj.results2$estimate*unadj.results2$iqr + 1.96*(unadj.results2$std.error*unadj.results2$iqr))
unadj.results2$ORCI = str_c(round(unadj.results2$OR,2), " (", round(unadj.results2$LCL,2), ", ", round(unadj.results2$UCL,2), ")")
unadj.results3 = unadj.results2 %>% pivot_wider(id_cols = c("term"),
                                                names_from = "y.level",
                                                values_from = c("ORCI"))
write.csv(unadj.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Unadjusted Multinomial Model Estimates with Weights.csv")
write.csv(unadj.results3, "Table of unadjusted Multinomial Model Estimates with Weights.csv")


adj.results2 = merge(adj.results, iqr.all, by = "term", all = T)
adj.results2 = adj.results2 %>% mutate(iqr = as.numeric(iqr))
adj.results2 = adj.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
adj.results2$OR = exp(adj.results2$estimate*adj.results2$iqr)
adj.results2$LCL = exp(adj.results2$estimate*adj.results2$iqr - 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$UCL = exp(adj.results2$estimate*adj.results2$iqr + 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$ORCI = str_c(round(adj.results2$OR,2), " (", round(adj.results2$LCL,2), ", ", round(adj.results2$UCL,2), ")")
adj.results3 = adj.results2 %>% pivot_wider(id_cols = c("term"),
                                                names_from = "y.level",
                                                values_from = c("ORCI"))
write.csv(adj.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Adjusted Multinomial Model Estimates with Weights.csv")
write.csv(adj.results3, "Table of Adjusted Multinomial Model Estimates with Weights.csv")



## Re-do this main model with different adjustments ##
21*length(confounders.wide)


#Adjusted models
adj.results <- data.frame()
lab.results <- data.frame()
all.results <- data.frame()
for(j in 1:length(confounders.wide)){
  for(i in 1:length(exposures)){
    results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide[-j]), collapse = "+"))), 
                                       weights = overall, na.omit=T))
    tidied_list <- lapply(results, function(mod) {
      tidy(mod) %>%
        mutate(term = paste(y.level, term, sep = ": "))
    })
    
    pooled <- pool(results)
    p.results = summary(pooled)
    p.results$term2 <- tidied_list[[1]]$term
    p.results <- p.results %>%
      separate(
        term2,
        into = c("y.level", "term"),
        sep = ": "
      )

    p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
    adj.results = rbind(adj.results, data.frame(p.results.term, not.adjusted.for = confounders.wide[j]))
  }}
adj.results

adj.results2 = merge(adj.results, iqr.all, by = "term", all = T)
adj.results2 = adj.results2 %>% mutate(iqr = as.numeric(iqr))
adj.results2 = adj.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
adj.results2$OR = exp(adj.results2$estimate*adj.results2$iqr)
adj.results2$LCL = exp(adj.results2$estimate*adj.results2$iqr - 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$UCL = exp(adj.results2$estimate*adj.results2$iqr + 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$ORCI = str_c(round(adj.results2$OR,2), " (", round(adj.results2$LCL,2), ", ", round(adj.results2$UCL,2), ")")
adj.results2 = adj.results2[,c("term","not.adjusted.for","y.level","ORCI")]
adj.results3 = adj.results2 %>% pivot_wider(.,
                                            names_from = "y.level",
                                            values_from = "ORCI")
write.csv(adj.results3, "Table of Adjusted Multinomial Model Estimates with Weights for Different Adjustments.csv")


###############################
## Secondary analysis by sex ##
###############################

#### Males
## Edit dataframe for analysis of sex-specific trajectories in males
imp.dat.weighted$class <- ifelse(imp.dat.weighted$male.class == '1', "Stable Growth, n=57",
                                 ifelse(imp.dat.weighted$male.class =='2', "Delayed Increase, n=8",
                                        ifelse(imp.dat.weighted$male.class =='3', "Decrease, n=21","Rapid Growth, n=15")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=57")
table(imp.dat.weighted$overall.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Limit to males only for analysis
male.mice = filter(imp.dat.weighted, imp.dat.weighted$childsgender == 2)
list.m.df = list()
for(i in 1:10){
  list.m.df[[i]] = male.mice[which(male.mice$.imp==i),]
}
imp.list.m.df = imputationList(list.m.df)

## We are doing adjusted models only for now
#Same covariate list as used in preliminary models - minus child's sex assigned at birth
confounders.wide.sex <- c('feedinggroup', 'mage', 'meduc', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')

male.results = data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.m.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide.sex), collapse="+"))),
                                         weights = male, na.action = na.omit))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  male.results = rbind(male.results, p.results.term)
}

male.results2 = merge(male.results, iqr.all, by = "term", all = T)
male.results2 = male.results2 %>% mutate(iqr = as.numeric(iqr))
male.results2 = male.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
male.results2$OR = exp(male.results2$estimate*male.results2$iqr)
male.results2$LCL = exp(male.results2$estimate*male.results2$iqr - 1.96*(male.results2$std.error*male.results2$iqr))
male.results2$UCL = exp(male.results2$estimate*male.results2$iqr + 1.96*(male.results2$std.error*male.results2$iqr))
male.results2$ORCI = str_c(round(male.results2$OR,2), " (", round(male.results2$LCL,2), ", ", round(male.results2$UCL,2), ")")
male.results3 = male.results2 %>% pivot_wider(id_cols = c("term"),
                                              names_from = "y.level",
                                              values_from = c("ORCI"))

write.csv(male.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Male Multinomial Model Estimates with Weights.csv")
write.csv(male.results3, "Table of Male Adjusted Multinomial Model Estimates with Weights.csv")



#### Females
## Edit dataframe for analysis of sex-specific trajectories in females
imp.dat.weighted$class <- ifelse(imp.dat.weighted$female.class == '1', "Stable Growth, n=30",
                                 ifelse(imp.dat.weighted$female.class =='2', "Delayed Increase, n=14",
                                        ifelse(imp.dat.weighted$female.class =='3', "Rapid Growth, n=30", "Decrease, n=24")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=30")
table(imp.dat.weighted$female.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Limit to females only for analysis
female.mice = filter(imp.dat.weighted, imp.dat.weighted$childsgender == 1)
list.f.df = list()
for(i in 1:10){
  list.f.df[[i]] = female.mice[which(female.mice$.imp==i),]
}
imp.list.f.df = imputationList(list.f.df)

## We are doing adjusted models only for now
#Same covariate list as used in preliminary models - minus child's sex assigned at birth
confounders.wide.sex <- c('feedinggroup', 'mage', 'meduc', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')
imp.list.f.df$imputations[[1]]
female.results = data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.f.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide.sex), collapse="+"))), 
                                         weights = female, na.action = na.omit))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  female.results = rbind(female.results, p.results.term)
}

female.results2 = merge(female.results, iqr.all, by = "term", all = T)
female.results2 = female.results2 %>% mutate(iqr = as.numeric(iqr))
female.results2 = female.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
female.results2$OR = exp(female.results2$estimate*female.results2$iqr)
female.results2$LCL = exp(female.results2$estimate*female.results2$iqr - 1.96*(female.results2$std.error*female.results2$iqr))
female.results2$UCL = exp(female.results2$estimate*female.results2$iqr + 1.96*(female.results2$std.error*female.results2$iqr))
female.results2$ORCI = str_c(round(female.results2$OR,2), " (", round(female.results2$LCL,2), ", ", round(female.results2$UCL,2), ")")
female.results3 = female.results2 %>% pivot_wider(id_cols = c("term"),
                                                  names_from = "y.level",
                                                  values_from = c("ORCI"))

write.csv(female.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Female Multinomial Model Estimates with Weights.csv")
write.csv(female.results3, "Table of Female Adjusted Multinomial Model Estimates with Weights.csv")



#### Females sensitivity analysis limited to before 28 weeks
## Edit dataframe for analysis of sex-specific trajectories in females
imp.dat.weighted$class <- ifelse(imp.dat.weighted$female28.class == '1', "Rapid Growth, n=22",
                                 ifelse(imp.dat.weighted$female28.class =='2', "Delayed Increase, n=12",
                                        ifelse(imp.dat.weighted$female28.class =='3', "Decrease, n=30", "Stable Growth, n=34")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=34")
table(imp.dat.weighted$female28.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Limit to females only for analysis
female.mice = filter(imp.dat.weighted, imp.dat.weighted$childsgender == 1)
list.f.df = list()
for(i in 1:10){
  list.f.df[[i]] = female.mice[which(female.mice$.imp==i),]
}
imp.list.f.df = imputationList(list.f.df)


## We are doing adjusted models only for now
#Same covariate list as used in preliminary models - minus child's sex assigned at birth
confounders.wide.sex <- c('feedinggroup', 'mage', 'meduc', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')
imp.list.f.df$imputations[[1]]
female.results = data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.f.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide.sex), collapse="+"))),
                                         weights = female28, na.action = na.omit))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  female.results = rbind(female.results, p.results.term)
}

female.results2 = merge(female.results, iqr.all, by = "term", all = T)
female.results2 = female.results2 %>% mutate(iqr = as.numeric(iqr))
female.results2 = female.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
female.results2$OR = exp(female.results2$estimate*female.results2$iqr)
female.results2$LCL = exp(female.results2$estimate*female.results2$iqr - 1.96*(female.results2$std.error*female.results2$iqr))
female.results2$UCL = exp(female.results2$estimate*female.results2$iqr + 1.96*(female.results2$std.error*female.results2$iqr))
female.results2$ORCI = str_c(round(female.results2$OR,2), " (", round(female.results2$LCL,2), ", ", round(female.results2$UCL,2), ")")
female.results3 = female.results2 %>% pivot_wider(id_cols = c("term"),
                                                  names_from = "y.level",
                                                  values_from = c("ORCI"))

write.csv(female.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Female Multinomial Model Estimates Limited to 28 Weeks with Weights.csv")
write.csv(female.results3, "Table of Female Adjusted Multinomial Model Estimates Limited to 28 Weeks with Weights.csv")



####################################
## Run stratified by child's diet ##
####################################

#### Cowsmilk
## Edit dataframe for analysis of sex-specific trajectories in cowsmilks
imp.dat.weighted$class <- ifelse(imp.dat.weighted$cowsmilk4.class == '1', "Stable Growth, n=45",
                                 ifelse(imp.dat.weighted$cowsmilk4.class =='2', "Decrease, n=21",
                                        ifelse(imp.dat.weighted$cowsmilk4.class == "3", "Delayed Increase, n=24", "Rapid Growth, n=19")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=45")
table(imp.dat.weighted$cowsmilk4.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Limit to cowsmilks only for analysis
cowsmilk.mice = filter(imp.dat.weighted, imp.dat.weighted$feedinggroup == 2)
list.c.df = list()
for(i in 1:10){
  list.c.df[[i]] = cowsmilk.mice[which(cowsmilk.mice$.imp==i),]
}
imp.list.c.df = imputationList(list.c.df)

## We are doing adjusted models only for now
#Same covariate list as used in preliminary models - minus child's diet
confounders.wide.diet <- c('childsgender', 'mage', 'meduc', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')

cowsmilk.results = data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.c.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide.diet), collapse="+"))),
                                         weights = cowsmilk4, na.action = na.omit))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  cowsmilk.results = rbind(cowsmilk.results, p.results.term)
}

cowsmilk.results2 = merge(cowsmilk.results, iqr.all, by = "term", all = T)
cowsmilk.results2 = cowsmilk.results2 %>% mutate(iqr = as.numeric(iqr))
cowsmilk.results2 = cowsmilk.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
cowsmilk.results2$OR = exp(cowsmilk.results2$estimate*cowsmilk.results2$iqr)
cowsmilk.results2$LCL = exp(cowsmilk.results2$estimate*cowsmilk.results2$iqr - 1.96*(cowsmilk.results2$std.error*cowsmilk.results2$iqr))
cowsmilk.results2$UCL = exp(cowsmilk.results2$estimate*cowsmilk.results2$iqr + 1.96*(cowsmilk.results2$std.error*cowsmilk.results2$iqr))
cowsmilk.results2$ORCI = str_c(round(cowsmilk.results2$OR,2), " (", round(cowsmilk.results2$LCL,2), ", ", round(cowsmilk.results2$UCL,2), ")")
cowsmilk.results3 = cowsmilk.results2 %>% pivot_wider(id_cols = c("term"),
                                                      names_from = "y.level",
                                                      values_from = c("ORCI"))

write.csv(cowsmilk.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Cowsmilk Multinomial Model Estimates with Weights.csv")
write.csv(cowsmilk.results3, "Table of Cowsmilk Adjusted Multinomial Model Estimates with Weights.csv")


#### breastmilk
## Edit dataframe for analysis of sex-specific trajectories in breastmilks
imp.dat.weighted$class <- ifelse(imp.dat.weighted$breastmilk.class == '1', "Decrease, n=22",
                                 ifelse(imp.dat.weighted$breastmilk.class =='2', "Rapid Growth, n=10",
                                        ifelse(imp.dat.weighted$breastmilk.class =='3', "Stable Growth, n=49","Delayed Increase, n=9")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=49")
table(imp.dat.weighted$breastmilk.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Limit to breastmilks only for analysis
breastmilk.mice = filter(imp.dat.weighted, imp.dat.weighted$feedinggroup == 1)
list.b.df = list()
for(i in 1:10){
  list.b.df[[i]] = breastmilk.mice[which(breastmilk.mice$.imp==i),]
}
imp.list.b.df = imputationList(list.b.df)

## We are doing adjusted models only for now
#Same covariate list as used in preliminary models - minus child's diet
confounders.wide.diet <- c('childsgender', 'mage', 'meduc', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')

breastmilk.results = data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.b.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide.diet), collapse="+"))), 
                                         weights = breastmilk, na.action = na.omit))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  breastmilk.results = rbind(breastmilk.results, p.results.term)
}

breastmilk.results2 = merge(breastmilk.results, iqr.all, by = "term", all = T)
breastmilk.results2 = breastmilk.results2 %>% mutate(iqr = as.numeric(iqr))
breastmilk.results2 = breastmilk.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
breastmilk.results2$OR = exp(breastmilk.results2$estimate*breastmilk.results2$iqr)
breastmilk.results2$LCL = exp(breastmilk.results2$estimate*breastmilk.results2$iqr - 1.96*(breastmilk.results2$std.error*breastmilk.results2$iqr))
breastmilk.results2$UCL = exp(breastmilk.results2$estimate*breastmilk.results2$iqr + 1.96*(breastmilk.results2$std.error*breastmilk.results2$iqr))
breastmilk.results2$ORCI = str_c(round(breastmilk.results2$OR,2), " (", round(breastmilk.results2$LCL,2), ", ", round(breastmilk.results2$UCL,2), ")")
breastmilk.results3 = breastmilk.results2 %>% pivot_wider(id_cols = c("term"),
                                                          names_from = "y.level",
                                                          values_from = c("ORCI"))

write.csv(breastmilk.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Breastmilk Multinomial Model Estimates with Weights.csv")
write.csv(breastmilk.results3, "Table of Breastmilk Adjusted Multinomial Model Estimates with Weights.csv")





######################################
## Sensitivity analysis using BMI Z ##
######################################

## Edit dataframe for analysis of nonorm trajectories
imp.dat.weighted$class <- ifelse(imp.dat.weighted$nonorm.class == '1', "Rapid Growth, n=45",
                                 ifelse(imp.dat.weighted$nonorm.class =='2', "Stable Growth, n=80",
                                        ifelse(imp.dat.weighted$nonorm.class =='3', "Delayed Increase, n=25", "Decrease, n=49")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=80")
table(imp.dat.weighted$nonorm.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Create imputation list
list.df = list()
for(i in 1:10){
  list.df[[i]] = imp.dat.weighted[which(imp.dat.weighted$.imp==i),]
}
imp.list.df = imputationList(list.df)

#nonormusted models
nonorm.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide), collapse = "+"))), 
                                       weights = nonorm, na.omit=T))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  nonorm.results = rbind(nonorm.results, p.results.term)
}

nonorm.results2 = merge(nonorm.results, iqr.all, by = "term", all = T)
nonorm.results2 = nonorm.results2 %>% mutate(iqr = as.numeric(iqr))
nonorm.results2 = nonorm.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
nonorm.results2$OR = exp(nonorm.results2$estimate*nonorm.results2$iqr)
nonorm.results2$LCL = exp(nonorm.results2$estimate*nonorm.results2$iqr - 1.96*(nonorm.results2$std.error*nonorm.results2$iqr))
nonorm.results2$UCL = exp(nonorm.results2$estimate*nonorm.results2$iqr + 1.96*(nonorm.results2$std.error*nonorm.results2$iqr))
nonorm.results2$ORCI = str_c(round(nonorm.results2$OR,2), " (", round(nonorm.results2$LCL,2), ", ", round(nonorm.results2$UCL,2), ")")
nonorm.results3 = nonorm.results2 %>% pivot_wider(id_cols = c("term"),
                                                  names_from = "y.level",
                                                  values_from = c("ORCI"))
write.csv(nonorm.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Non-normalized Multinomial Model Estimates with Weights.csv")
write.csv(nonorm.results3, "Table of Non-normalized Multinomial Model Estimates with Weights.csv")





#################################################
## Sensitivity analysis using weight-for-age Z ##
#################################################

## Edit dataframe for analysis of waz trajectories
imp.dat.weighted$class <- ifelse(imp.dat.weighted$waz.class == '1', "Stable Growth, n=88",
                                 ifelse(imp.dat.weighted$waz.class =='2', "Decrease, n=66",
                                        ifelse(imp.dat.weighted$waz.class =='3', "Rapid Growth, n=30", "Delayed Increase, n=15")))


imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=88")
table(imp.dat.weighted$waz.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## Create imputation list
list.df = list()
for(i in 1:10){
  list.df[[i]] = imp.dat.weighted[which(imp.dat.weighted$.imp==i),]
}
imp.list.df = imputationList(list.df)

#waz models
waz.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposures.obrien.ln[i], confounders.wide), collapse = "+"))), 
                                       weights = waz, na.omit=T))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.obrien.ln[i]))
  waz.results = rbind(waz.results, p.results.term)
}

waz.results2 = merge(waz.results, iqr.all, by = "term", all = T)
waz.results2 = waz.results2 %>% mutate(iqr = as.numeric(iqr))
waz.results2 = waz.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
waz.results2$OR = exp(waz.results2$estimate*waz.results2$iqr)
waz.results2$LCL = exp(waz.results2$estimate*waz.results2$iqr - 1.96*(waz.results2$std.error*waz.results2$iqr))
waz.results2$UCL = exp(waz.results2$estimate*waz.results2$iqr + 1.96*(waz.results2$std.error*waz.results2$iqr))
waz.results2$ORCI = str_c(round(waz.results2$OR,2), " (", round(waz.results2$LCL,2), ", ", round(waz.results2$UCL,2), ")")
waz.results3 = waz.results2 %>% pivot_wider(id_cols = c("term"),
                                            names_from = "y.level",
                                            values_from = c("ORCI"))
write.csv(waz.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Weight-for-Age Multinomial Model Estimates with Weights.csv")
write.csv(waz.results3, "Table of Weight-for-Age Multinomial Model Estimates with Weights.csv")



##########################################################
## Run models associating phenols with group assignment ##
## Use Visit 6-8 only                                   ##
##########################################################

## Edit dataframe for analysis of overall trajectories
imp.dat.weighted$class <- ifelse(imp.dat.weighted$overall.class == '1', "Decrease, n=43",
                                 ifelse(imp.dat.weighted$overall.class =='2', "Delayed Increase, n=27",
                                        ifelse(imp.dat.weighted$overall.class =='3', "Rapid Growth, n=47","Stable Growth, n=82")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=82")
table(imp.dat.weighted$overall.class)
table(imp.dat.weighted$class) #Check by dividing n's by 10 to get the n= for each group

## log-transform exposure
exposure.obrien.v68.ln <- paste0(exposures.obrien.v68, '_ln')
exposure.obrien.v12.ln <- paste0(exposures.obrien.v12, '_ln')

for(i in 1:length(exposures)){ 
  imp.dat.weighted[exposure.obrien.v68.ln[i]] = log(imp.dat.weighted[exposures.obrien.v68[i]])
  imp.dat.weighted[exposure.obrien.v12.ln[i]] = log(imp.dat.weighted[exposures.obrien.v12[i]])
}

## Fix the few inf and NaNs that this causes
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
imp.dat.weighted <- replace(imp.dat.weighted, is.nan(imp.dat.weighted), NA)
imp.dat.weighted[imp.dat.weighted==Inf] = NA
imp.dat.weighted[imp.dat.weighted==-Inf] = NA


## Create imputation list
list.df = list()
for(i in 1:10){
  list.df[[i]] = imp.dat.weighted[which(imp.dat.weighted$.imp==i),]
}
imp.list.df = imputationList(list.df)

## Test model
results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c('bpa_obrien_v68_ln', confounders.wide), collapse = "+"))), 
                                     weights = overall, na.omit=T))
tidied_list <- lapply(results, function(mod) {
  tidy(mod) %>%
    mutate(term = paste(y.level, term, sep = ": "))
})

pooled <- pool(results)
p.results = summary(pooled)
p.results$term2 <- tidied_list[[1]]$term
p.results <- p.results %>%
  separate(
    term2,
    into = c("y.level", "term"),
    sep = ": "
  )
p.results.term = p.results %>% filter(term == paste('bpa_obrien_v68_ln'))

#Adjusted models at 6-8 weeks
adj.v68.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposure.obrien.v68.ln[i], confounders.wide), collapse = "+"))), 
                                       weights = overall, na.omit=T))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposure.obrien.v68.ln[i]))
  adj.v68.results = rbind(adj.v68.results, p.results.term)
}


#Adjusted models at 12 weeks
adj.v12.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposure.obrien.v12.ln[i], confounders.wide), collapse = "+"))), 
                                       weights = overall, na.omit=T))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposure.obrien.v12.ln[i]))
  adj.v12.results = rbind(adj.v12.results, p.results.term)
}

iqr.all$variable = exposure.obrien.v68.ln
adj.v68.results$variable = adj.v68.results$term
adj.results2 = merge(adj.v68.results, iqr.all, by = "variable", all = T)
adj.results2 = adj.results2 %>% mutate(iqr = as.numeric(iqr))
adj.results2 = adj.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
adj.results2$OR = exp(adj.results2$estimate*adj.results2$iqr)
adj.results2$LCL = exp(adj.results2$estimate*adj.results2$iqr - 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$UCL = exp(adj.results2$estimate*adj.results2$iqr + 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$ORCI = str_c(round(adj.results2$OR,2), " (", round(adj.results2$LCL,2), ", ", round(adj.results2$UCL,2), ")")
adj.results3 = adj.results2 %>% pivot_wider(id_cols = c("variable"),
                                            names_from = "y.level",
                                            values_from = c("ORCI"))
write.csv(adj.results2 %>% select(c("variable", "y.level", "OR", "LCL", "UCL")), "Adjusted Multinomial Model Estimates at the 68 Week Visit with Weights.csv")
write.csv(adj.results3, "Table of Adjusted Multinomial Model Estimates at the 68 Week Visit with Weights.csv")


iqr.all$variable = exposure.obrien.v12.ln
adj.v12.results$variable = adj.v12.results$term
adj.results2 = merge(adj.v12.results, iqr.all, by = "variable", all = T)
adj.results2 = adj.results2 %>% mutate(iqr = as.numeric(iqr))
adj.results2 = adj.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
adj.results2$OR = exp(adj.results2$estimate*adj.results2$iqr)
adj.results2$LCL = exp(adj.results2$estimate*adj.results2$iqr - 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$UCL = exp(adj.results2$estimate*adj.results2$iqr + 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$ORCI = str_c(round(adj.results2$OR,2), " (", round(adj.results2$LCL,2), ", ", round(adj.results2$UCL,2), ")")
adj.results3 = adj.results2 %>% pivot_wider(id_cols = c("variable"),
                                            names_from = "y.level",
                                            values_from = c("ORCI"))
write.csv(adj.results2 %>% select(c("variable", "y.level", "OR", "LCL", "UCL")), "Adjusted Multinomial Model Estimates at the 12 Week Visit with Weights.csv")
write.csv(adj.results3, "Table of Adjusted Multinomial Model Estimates at the 12 Week Visit with Weights.csv")




##########################################################
## Run models associating phenols with group assignment ##
## Use Non-Creatinine Standardized Measures             ##
##########################################################
imp.dat.weighted$class <- ifelse(imp.dat.weighted$overall.class == '1', "Decrease, n=43",
                                 ifelse(imp.dat.weighted$overall.class =='2', "Delayed Increase, n=27",
                                        ifelse(imp.dat.weighted$overall.class =='3', "Rapid Growth, n=47","Stable Growth, n=82")))

imp.dat.weighted$class <- relevel(as.factor(imp.dat.weighted$class), ref = "Stable Growth, n=82")
table(imp.dat.weighted$overall.class)
table(imp.dat.weighted$class)/11/4 #Check by dividing n's by 11 to get the n= for each group

## Create imputation list
list.df = list()
for(i in 1:10){
  list.df[[i]] = imp.dat.weighted[which(imp.dat.weighted$.imp==i),]
}
imp.list.df = imputationList(list.df)

## Test model
results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c('bpa_gm_ln', confounders.wide), collapse = "+"))), 
                                     weights = overall, na.omit=T))
tidied_list <- lapply(results, function(mod) {
  tidy(mod) %>%
    mutate(term = paste(y.level, term, sep = ": "))
})

pooled <- pool(results)
p.results = summary(pooled)
p.results$term2 <- tidied_list[[1]]$term
p.results <- p.results %>%
  separate(
    term2,
    into = c("y.level", "term"),
    sep = ": "
  )
p.results.term = p.results %>% filter(term == paste('bpa_gm_ln'))

#Unadjusted models
unadj.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", exposures.gm.ln[i])), na.omit=T))
  
  tidied_list <- lapply(results, function(mod) {
    tidy(mod) %>%
      mutate(term = paste(y.level, term, sep = ": "))
  })
  
  pooled <- pool(results)
  p.results = summary(pooled)
  p.results$term2 <- tidied_list[[1]]$term
  p.results <- p.results %>%
    separate(
      term2,
      into = c("y.level", "term"),
      sep = ": "
    )
  p.results.term = p.results %>% filter(term == paste(exposures.gm.ln[i]))
  unadj.results = rbind(unadj.results, p.results.term)
}

#Adjusted models
adj.results <- data.frame()
for(i in 1:length(exposures)){
  results = with(imp.list.df, multinom(as.formula(paste("class", "~", paste(c(exposures.gm.ln[i], confounders.wide), collapse = "+"))), 
                                       weights = overall, na.omit=T))
  
 tidied_list <- lapply(results, function(mod) {
  tidy(mod) %>%
    mutate(term = paste(y.level, term, sep = ": "))
})

pooled <- pool(results)
p.results = summary(pooled)
p.results$term2 <- tidied_list[[1]]$term
p.results <- p.results %>%
  separate(
    term2,
    into = c("y.level", "term"),
    sep = ": "
  )  
p.results.term = p.results %>% filter(term == paste(exposures.gm.ln[i]))
  adj.results = rbind(adj.results, p.results.term)
}


iqr.all.2 <- iqr.all
iqr.all.2$term <- exposures.gm.ln

adj.results2 = merge(adj.results, iqr.all.2, by = "term", all = T)
adj.results2 = adj.results2 %>% mutate(iqr = as.numeric(iqr))
adj.results2 = adj.results2 %>% mutate(iqr = ifelse(is.na(iqr), 1, iqr))
adj.results2$OR = exp(adj.results2$estimate*adj.results2$iqr)
adj.results2$LCL = exp(adj.results2$estimate*adj.results2$iqr - 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$UCL = exp(adj.results2$estimate*adj.results2$iqr + 1.96*(adj.results2$std.error*adj.results2$iqr))
adj.results2$ORCI = str_c(round(adj.results2$OR,2), " (", round(adj.results2$LCL,2), ", ", round(adj.results2$UCL,2), ")")
adj.results3 = adj.results2 %>% pivot_wider(id_cols = c("term"),
                                            names_from = "y.level",
                                            values_from = c("ORCI"))
write.csv(adj.results2 %>% select(c("term", "y.level", "OR", "LCL", "UCL")), "Adjusted Multinomial Model Estimates For Non-Creatinine Adjusted Phenols with Weights.csv")
write.csv(adj.results3, "Table of Adjusted Multinomial Model Estimates For Non-Creatinine Adjusted Phenols with Weights.csv")




###################
###################
## QGCOMP Models ##
###################
###################


###########################################
## Primary analysis in overall dataframe ##
###########################################

## Convert estimates by dividing by the IQR
exposures.obrien.iqr <- paste0(exposures, '_obrien_iqr')
exposures.obrien.iqr.v68 <- paste0(exposures.obrien.iqr, '_v68')
exposures.obrien.iqr.v12 <- paste0(exposures.obrien.iqr, '_v12')
exposures.iqr <- paste0(exposures, '_iqr')

for(i in 1:length(exposures)){
  print(exposures.obrien.ln[i])
  print(iqr.all[i,2])
  # blank exposure i
  imp.dat.weighted[[exposures.obrien.iqr[i]]] <- NA
  imp.dat.weighted[[exposures.obrien.iqr.v68[i]]] <- NA
  imp.dat.weighted[[exposures.obrien.iqr.v12[i]]] <- NA
  imp.dat.weighted[[exposures.iqr[i]]] <- NA
  # input values for .imp j
  for(j in 0:10){ 
    print(paste0("imputation = ", j))
    # divide by IQR 
    imp.dat.weighted[which(imp.dat.weighted$.imp==j),c(exposures.obrien.iqr[i])] <- 
      imp.dat.weighted[which(imp.dat.weighted$.imp==j), exposures.obrien.ln[i]]/as.numeric(iqr.all[i,2])
    # divide by IQR - use same iqr as main analysis and not visit-specific to keep consistent
    imp.dat.weighted[which(imp.dat.weighted$.imp==j),c(exposures.obrien.iqr.v68[i])] <- 
      imp.dat.weighted[which(imp.dat.weighted$.imp==j), exposure.obrien.v68.ln[i]]/as.numeric(iqr.all[i,2])
    # divide by IQR - use same iqr as main analysis and not visit-specific to keep consistent
    imp.dat.weighted[which(imp.dat.weighted$.imp==j),c(exposures.obrien.iqr.v12[i])] <- 
      imp.dat.weighted[which(imp.dat.weighted$.imp==j), exposure.obrien.v12.ln[i]]/as.numeric(iqr.all[i,2])
    # divide by IQR - use same iqr as main analysis and not the non-creatinine adjusted to keep consistent
    imp.dat.weighted[which(imp.dat.weighted$.imp==j),c(exposures.iqr[i])] <- 
      imp.dat.weighted[which(imp.dat.weighted$.imp==j), exposures.gm.ln[i]]/as.numeric(iqr.all[i,2])
  }
}
summary(imp.dat.weighted[exposures.obrien.iqr])
view(imp.dat.weighted[1:30,c('bpa_obrien_gm_ln','bpa_obrien_iqr',
                             'bpa_obrien_v68_ln','bpa_obrien_iqr_v68',
                             'bpa_obrien_v12_ln','bpa_obrien_iqr_v12',
                             'bpa_gm_ln','bpa_iqr')])

## No. of imputations 
imp.length <- 10  

## Outcomes
## Edit dataframe for analysis of overall trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$overall.class == '2', 1,
                                ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$overall.class == '1', 1,
                                ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$overall.class == '3', 1,
                                ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Test function
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.boot(f = as.formula(paste('decrease', "~", paste(exposures.obrien.iqr, collapse = "+"))),
                                       expnms = exposures.obrien.iqr, 
                                       family = binomial(), 
                                       weights = overall,
                                       q = NULL, bayes = TRUE, seed=125, rr=TRUE,
                                       data = imp.dat.weighted.decrease[which(imp.dat.weighted.decrease$.imp==x),])))

qg.fit.imp #Gives overview of results for each imputation

# Pool across imputations using mitools package (on my computer, pool() function won't work)
# We want to multiply the estimate and the SEs by 2 to approximate the IQR estimates from our single-pollutant models
MIcombineRes<-MIcombine(qg.fit.imp$analyses)

MIcombinevals <- function(MIcombineRes,outcome) {
  est <- MIcombineRes$coefficients
  se <- sqrt(diag(MIcombineRes$variance))
  x <- 1 - 0.05/2
  tStat <- MIcombineRes$coefficients/sqrt(diag(MIcombineRes$variance))
  pval <- 2*pt(-abs(tStat),df=MIcombineRes$df)
  lcl <- (MIcombineRes$coefficients) - qt(x,df=MIcombineRes$df)*(se)
  ucl <- (MIcombineRes$coefficients) + qt(x,df=MIcombineRes$df)*(se)
  est = exp(est)
  lcl = exp(lcl)
  ucl = exp(ucl)
  estci1 <- str_c(round(est,1), " (", round(lcl,1), ", ", round(ucl,1), ")")
  estci2 <- str_c(round(est,2), " (", round(lcl,2), ", ", round(ucl,2), ")")
  label <- outcome
  results <- cbind(outcome, est, lcl, ucl, estci1, estci2, pval)
  results[2,]
} #Calculate the lcl and ucl based on https://bookdown.org/mwheymans/bookmi/rubins-rules.html#degrees-of-freedom-and-p-values
#Pval calculation taken from existing f'n at https://thestatsgeek.com/2020/11/05/p-values-after-multiple-imputation-using-mitools-in-r/ 
qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),'outcome_mod')
qg.psi

## Set-up Function: unadjusted models
myqgcompfit <- function(dat, outcome_mod){
  
  dat <- dat
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = as.formula(paste(outcome_mod, "~", paste(exposures.obrien.iqr, collapse = "+"))),
                                              expnms = exposures.obrien.iqr, 
                                              family = binomial(), 
                                              weights = overall,
                                              q = NULL, bayes = TRUE, seed=125, rr=TRUE,
                                              data = dat[which(dat$.imp==x),])))
  
  #Output parameter estimates for psi1
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),outcome_mod)
  return(qg.psi)
}


delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc")
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease") 
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid") 

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Unadjusted QGCOMP RR Model Estimates with Weights.csv")
write.csv(out.results.2, "Table of Unadjusted QGCOMP RR Model Estimates with Weights.csv")

## Set-up Function: adjusted models
myqgcompfit <- function(dat, outcome_mod, covariate_mod, exposure_mod, weight_mod){
  
  dat <- dat
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = as.formula(paste(outcome_mod, "~", paste(c(exposure_mod, covariate_mod), collapse = "+"))),
                                              data = dat[which(dat$.imp==x),],
                                              expnms = exposure_mod, 
                                              family = binomial(), 
                                              weights = as.matrix(dat[which(dat$.imp==x), weight_mod]),
                                              q = NULL, bayes = TRUE, seed=125, rr=TRUE)))
  
  #Output parameter estimates for psi1
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),outcome_mod)
  return(qg.psi)
}

dat.weights <- as.matrix(imp.dat.weighted.delayinc[which(imp.dat.weighted.delayinc$.imp==1), 'overall'])


delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide, exposures.obrien.iqr, 'overall') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide, exposures.obrien.iqr, 'overall') 
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide, exposures.obrien.iqr, 'overall') 

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Adjusted QGCOMP RR Model Estimates with Weights.csv")
write.csv(out.results.2, "Table of Adjusted QGCOMP RR Model Estimates with Weights.csv")


####################################################################
## Sensitivity analysis using non-creatinine standardized phenols ##
####################################################################
dat.weights <- as.matrix(imp.dat.weighted.delayinc[which(imp.dat.weighted.delayinc$.imp==1), 'overall'])

delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide, exposures.iqr, 'overall') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide, exposures.iqr, 'overall') 
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide, exposures.iqr, 'overall') 

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Adjusted QGCOMP RR Model Estimates for Non-Creatinine Standardized Phenols with Weights.csv")
write.csv(out.results.2, "Table of Adjusted QGCOMP RR Model Estimates for Non-Creatinine Standardized Phenols with Weights.csv")


###############################
## Secondary analysis by sex ##
###############################

#### Males
## Edit dataframe for analysis of overall trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$male.class == '2', 1,
                              ifelse(imp.dat.weighted$male.class =='1', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$male.class == '3', 1,
                                ifelse(imp.dat.weighted$male.class =='1', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$male.class == '4', 1,
                             ifelse(imp.dat.weighted$male.class =='1', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide.sex, exposures.obrien.iqr, 'male') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide.sex, exposures.obrien.iqr, 'male') 
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide.sex, exposures.obrien.iqr, 'male') 

m.out.results.1 = rbind(c("Males",decrease),c("Males",delayinc),c("Males",rapid))
m.out.results.2 = t(m.out.results.1)[c('outcome','estci2'),]
write.csv(m.out.results.1, "Adjusted QGCOMP RR Model Estimates in Males with Weights.csv")
write.csv(m.out.results.2, "Table of Adjusted QGCOMP RR Model Estimates in Males with Weights.csv")


#### Females
## Edit dataframe for analysis of overall trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$female.class == '2', 1,
                              ifelse(imp.dat.weighted$female.class =='1', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$female.class == '4', 1,
                                ifelse(imp.dat.weighted$female.class =='1', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$female.class == '3', 1,
                             ifelse(imp.dat.weighted$female.class =='1', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide.sex, exposures.obrien.iqr, 'female') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide.sex, exposures.obrien.iqr, 'female')  
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide.sex, exposures.obrien.iqr, 'female') 

f.out.results.1 = rbind(c("Females",decrease),c("Females",delayinc),c("Females",rapid))
f.out.results.2 = t(f.out.results.1)[c('outcome','estci2'),]
write.csv(f.out.results.1, "Adjusted QGCOMP RR Model Estimates in Females with Weights.csv")
write.csv(f.out.results.2, "Table of Adjusted QGCOMP RR Model Estimates in Females with Weights.csv")


#### Females Sensitivity Analysis for <28 weeks
## Edit dataframe for analysis of overall trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$female28.class == '2', 1,
                                ifelse(imp.dat.weighted$female28.class =='4', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$female28.class == '3', 1,
                                ifelse(imp.dat.weighted$female28.class =='4', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$female28.class == '1', 1,
                             ifelse(imp.dat.weighted$female28.class =='4', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide.sex, exposures.obrien.iqr, 'female28') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide.sex, exposures.obrien.iqr, 'female28')
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide.sex, exposures.obrien.iqr, 'female28')

f.out.results.1 = rbind(c("Females",decrease),c("Females",delayinc),c("Females",rapid))
f.out.results.2 = t(f.out.results.1)[c('outcome','estci2'),]
write.csv(f.out.results.1, "Adjusted QGCOMP RR Model Estimates in Females Sensitivity Analysis for 28 Weeks with Weights.csv")
write.csv(f.out.results.2, "Table of Adjusted QGCOMP RR Model Estimates in Females Sensitivity Analysis for 28 Weeks with Weights.csv")



########################################
## Secondary analysis by Child's diet ##
########################################

#### Cowsmilk
## Edit dataframe for analysis of trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$cowsmilk4.class == '3', 1,
                              ifelse(imp.dat.weighted$cowsmilk4.class =='1', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$cowsmilk4.class == '2', 1,
                                ifelse(imp.dat.weighted$cowsmilk4.class =='1', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$cowsmilk4.class == '4', 1,
                             ifelse(imp.dat.weighted$cowsmilk4.class =='1', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide.diet, exposures.obrien.iqr, 'cowsmilk4') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide.diet, exposures.obrien.iqr, 'cowsmilk4')  
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide.diet, exposures.obrien.iqr, 'cowsmilk4') 

c.out.results.1 = rbind(c("Cows-milk",decrease),c("Cows-milk",delayinc),c("Cows-milk",rapid))
c.out.results.2 = t(c.out.results.1)[c('outcome','estci2'),]
write.csv(c.out.results.1, "Adjusted QGCOMP RR Model Estimates in Cows-milk Fed Infants with Weights.csv")
write.csv(c.out.results.2, "Table of Adjusted QGCOMP RR Model Estimates in Cows-milk Fed Infants with Weights.csv")



#### breastmilk
## Edit dataframe for analysis of trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$breastmilk.class == '4', 1,
                              ifelse(imp.dat.weighted$breastmilk.class =='3', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$breastmilk.class == '1', 1,
                                ifelse(imp.dat.weighted$breastmilk.class =='3', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$breastmilk.class == '2', 1,
                             ifelse(imp.dat.weighted$breastmilk.class =='3', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide.diet, exposures.obrien.iqr, 'breastmilk') 
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide.diet, exposures.obrien.iqr, 'breastmilk') 
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide.diet, exposures.obrien.iqr, 'breastmilk') 

b.out.results.1 = rbind(c("Breastmilk",decrease),c("Breastmilk",delayinc),c("Breastmilk",rapid))
b.out.results.2 = t(b.out.results.1)[c('outcome','estci2'),]
write.csv(b.out.results.1, "Adjusted QGCOMP RR Model Estimates in Breastmilk Fed Infants with Weights.csv")
write.csv(b.out.results.2, "Table of Adjusted QGCOMP RR Model Estimates in Breastmilk Fed Infants with Weights.csv")



######################################
## Sensitivity analysis using BMI Z ##
######################################
## Edit dataframe for analysis of trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$nonorm.class == '3', 1,
                              ifelse(imp.dat.weighted$nonorm.class =='2', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$nonorm.class == '4', 1,
                                ifelse(imp.dat.weighted$nonorm.class =='2', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$nonorm.class == '1', 1,
                             ifelse(imp.dat.weighted$nonorm.class =='2', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide, exposures.obrien.iqr, 'nonorm')
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide, exposures.obrien.iqr, 'nonorm')
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide, exposures.obrien.iqr, 'nonorm')

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Adjusted QGCOMP RR Model Estimates for BMI Z-score with Weights.csv")
write.csv(out.results.2, "Table of Adjusted QGCOMP RR Model Estimates for BMI Z-score with Weights.csv")




#################################################
## Sensitivity analysis using weight-for-age Z ##
#################################################
## Edit dataframe for analysis of trajectories
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$waz.class == '4', 1,
                              ifelse(imp.dat.weighted$waz.class =='1', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$waz.class == '2', 1,
                                ifelse(imp.dat.weighted$waz.class =='1', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$waz.class == '3', 1,
                             ifelse(imp.dat.weighted$waz.class =='1', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

table(imp.dat.weighted$delayinc)
table(imp.dat.weighted$decrease)
table(imp.dat.weighted$rapid)

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposures.obrien.iqr)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposures.obrien.iqr)]),] ## complete dataset needed

## Run models
delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide, exposures.obrien.iqr, 'waz')
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide, exposures.obrien.iqr, 'waz')
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide, exposures.obrien.iqr, 'waz')

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Adjusted QGCOMP RR Model Estimates for Weight-for-Age Z-score with Weights.csv")
write.csv(out.results.2, "Table of Adjusted QGCOMP RR Model Estimates for Weight-for-Age Z-score with Weights.csv")





#####################################################
## Sensitivity analysis for exposures at 6-8 weeks ##
#####################################################
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$overall.class == '2', 1,
                                    ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$overall.class == '1', 1,
                                    ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$overall.class == '3', 1,
                                 ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposure.obrien.v68.ln)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposure.obrien.v68.ln)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposure.obrien.v68.ln)]),] ## complete dataset needed

## Sample size?
table(subset(imp.dat.weighted.delayinc$delayinc, imp.dat.weighted.delayinc$.imp>0)) #27
table(subset(imp.dat.weighted.decrease$decrease, imp.dat.weighted.decrease$.imp>0)) #34
table(subset(imp.dat.weighted.rapid$rapid, imp.dat.weighted.rapid$.imp>0)) #38
# Ref (Stable Growth): 77

delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide, exposure.obrien.v68.ln, 'overall')
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide, exposure.obrien.v68.ln, 'overall') 
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide, exposure.obrien.v68.ln, 'overall')

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Adjusted QGCOMP RR Model Estimates for Exposures at 6-8 Weeks with Weights.csv")
write.csv(out.results.2, "Table of Adjusted QGCOMP RR Model Estimates for Exposures at 6-8 Weeks with Weights.csv")



####################################################
## Sensitivity analysis for exposures at 12 weeks ##
####################################################
imp.dat.weighted$delayinc <- ifelse(imp.dat.weighted$overall.class == '2', 1,
                                    ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is Delayed Increase, 0 if group is Stable Growth
imp.dat.weighted$decrease <- ifelse(imp.dat.weighted$overall.class == '1', 1,
                                    ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is decrease, 0 if group is Stable Growth
imp.dat.weighted$rapid <- ifelse(imp.dat.weighted$overall.class == '3', 1,
                                 ifelse(imp.dat.weighted$overall.class =='4', 0, NA)) ## 1 if group is Rapid Growth, 0 if group is Stable Growth

imp.dat.weighted.delayinc <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('delayinc',exposure.obrien.v12.ln)]),] ## complete dataset needed
imp.dat.weighted.decrease <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('decrease', exposure.obrien.v12.ln)]),] ## complete dataset needed
imp.dat.weighted.rapid <- imp.dat.weighted[complete.cases(imp.dat.weighted[,c('rapid', exposure.obrien.v12.ln)]),] ## complete dataset needed

## Sample size?
table(subset(imp.dat.weighted.delayinc$delayinc, imp.dat.weighted.delayinc$.imp>0)) #21
table(subset(imp.dat.weighted.decrease$decrease, imp.dat.weighted.decrease$.imp>0)) #34
table(subset(imp.dat.weighted.rapid$rapid, imp.dat.weighted.rapid$.imp>0)) #37
# Ref (Stable Growth): 66

delayinc   <- myqgcompfit(imp.dat.weighted.delayinc, "delayinc", confounders.wide, exposure.obrien.v12.ln, 'overall')
decrease   <- myqgcompfit(imp.dat.weighted.decrease, "decrease", confounders.wide, exposure.obrien.v12.ln, 'overall')
rapid   <- myqgcompfit(imp.dat.weighted.rapid, "rapid", confounders.wide, exposure.obrien.v12.ln, 'overall')

out.results.1 = rbind(decrease,delayinc,rapid)
out.results.2 = t(out.results.1)[c('outcome','estci2'),]
write.csv(out.results.1, "Adjusted QGCOMP RR Model Estimates for Exposures at 12 Weeks with Weights.csv")
write.csv(out.results.2, "Table of Adjusted QGCOMP RR Model Estimates for Exposures at 12 Weeks with Weights.csv")




#####################################################
## Run models associating phenols with mean growth ##
#####################################################

############################################
## Single-pollutant models are LME models ##
## with random intercepts                 ##
## and a cubic b-spline with 4 knots      ##
## These are run in SAS                   ##
## with filename "IFED Cleaning_Final"    ##
############################################



###################
###################
## QGCOMP Models ##
###################
###################


###########################################
## Primary analysis in overall dataframe ##
###########################################
## Setting parameters for main models
imp.length <- 10  #No. of imputations 

## Dataframe
imp.dat.long = read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\Imputed dataset with outcomes long.csv")

## Convert estimates by dividing by the IQR
exposures.obrien.iqr <- paste0(exposures, '_iqr')
exposures.obrien.iqr.v68 <- paste0(exposures.obrien.iqr, '_v68')
exposures.obrien.iqr.v12 <- paste0(exposures.obrien.iqr, '_v12')

for(i in 1:length(exposures)){
  print(exposures.obrien.ln[i])
  print(iqr.all[i,2])
  # blank exposure i
  imp.dat.long[[exposures.obrien.iqr[i]]] <- NA
  imp.dat.long[[exposures.obrien.iqr.v68[i]]] <- NA
  imp.dat.long[[exposures.obrien.iqr.v12[i]]] <- NA
  # input values for .imp j
  for(j in 0:10){ 
    print(paste0("imputation = ", j))
    # divide by IQR 
    imp.dat.long[which(imp.dat.long$.imp==j),c(exposures.obrien.iqr[i])] <- 
      imp.dat.long[which(imp.dat.long$.imp==j), exposures.obrien.ln[i]]/as.numeric(iqr.all[i,2])
    # divide by IQR - use same iqr as main analysis and not visit-specific to keep consistent
    imp.dat.long[which(imp.dat.long$.imp==j),c(exposures.obrien.iqr.v68[i])] <- 
      imp.dat.long[which(imp.dat.long$.imp==j), exposure.obrien.v68.ln[i]]/as.numeric(iqr.all[i,2])
    # divide by IQR - use same iqr as main analysis and not visit-specific to keep consistent
    imp.dat.long[which(imp.dat.long$.imp==j),c(exposures.obrien.iqr.v12[i])] <- 
      imp.dat.long[which(imp.dat.long$.imp==j), exposure.obrien.v12.ln[i]]/as.numeric(iqr.all[i,2])
  }
}
summary(imp.dat.long[exposures.obrien.iqr])
view(imp.dat.long[1:30,c('bpa_obrien_gm_ln','bpa_iqr',
                             'bpa_obrien_v68_ln','bpa_iqr_v68',
                             'bpa_obrien_v12_ln','bpa_iqr_v12')])


## Set up function for combining qgcomp estimates
MIcombinevals <- function(MIcombineRes,outcome) {
  est <- MIcombineRes$coefficients
  se <- sqrt(diag(MIcombineRes$variance))
  x <- 1 - 0.05/2
  tStat <- MIcombineRes$coefficients/sqrt(diag(MIcombineRes$variance))
  pval <- 2*pt(-abs(tStat),df=MIcombineRes$df)
  lcl <- (MIcombineRes$coefficients) - qt(x,df=MIcombineRes$df)*(se)
  ucl <- (MIcombineRes$coefficients) + qt(x,df=MIcombineRes$df)*(se)
  estci1 <- str_c(round(est,1), " (", round(lcl,1), ", ", round(ucl,1), ")")
  estci2 <- str_c(round(est,2), " (", round(lcl,2), ", ", round(ucl,2), ")")
  label <- outcome
  results <- cbind(outcome, est, lcl, ucl, estci1, estci2, pval)
  results[2,]
} #Calculate the lcl and ucl based on https://bookdown.org/mwheymans/bookmi/rubins-rules.html#degrees-of-freedom-and-p-values
#Pval calculation taken from existing f'n at https://thestatsgeek.com/2020/11/05/p-values-after-multiple-imputation-using-mitools-in-r/ 


## Getting data ready for analysis 
dat1 <- imp.dat.long[complete.cases(imp.dat.long[,c('.imp','ifedid','bmiz','ageinweeks', confounders.wide, exposures.obrien.iqr)]),] ## complete dataset needed

## Test function
## Ran 3 times:
### 1) exposures.obrien.iqr - iqr-standardized exposures - with q=NULL
### 2) exposures.obrien.iqr - iqr-standardized exposures - with q=4
### 3) exposures.obrien.ln - log-transformed exposures - with q=4

# ## One
# qgformula <- as.formula(paste('bmiz', "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
#                               paste( c(exposures.obrien.iqr, confounders.wide), collapse = "+")))
# qg.fit.imp <- list(
#   analyses = lapply(1:imp.length, 
#                     function(x) qgcomp.boot(f = qgformula,
#                                             expnms = exposures.obrien.iqr, 
#                                             family = gaussian(), 
#                                             id = 'ifedid',
#                                             q = NULL, seed=125, 
#                                             data = dat1[which(dat1$.imp==x),])))
# 
# qg.fit.imp 
# qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),'outcome_mod')
# 
# ## Two
# qgformula <- as.formula(paste('bmiz', "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
#                               paste( c(exposures.obrien.iqr, confounders.wide), collapse = "+")))
# qg.fit.imp.4 <- list(
#   analyses = lapply(1:imp.length, 
#                     function(x) qgcomp.boot(f = qgformula,
#                                             expnms = exposures.obrien.iqr, 
#                                             family = gaussian(), 
#                                             id = 'ifedid',
#                                             q = 4, seed=125, 
#                                             data = dat1[which(dat1$.imp==x),])))
# 
# qg.fit.imp.4
# qg.psi.4 <- MIcombinevals(MIcombine(qg.fit.imp.4$analyses),'outcome_mod') # Used exposures.obrien.iqr and q=4
# 
# ## Three
# dat1 <- imp.dat.long[complete.cases(imp.dat.long[,c('.imp','ifedid','bmiz','ageinweeks', confounders.wide, exposures.obrien.ln)]),] ## complete dataset needed
# 
# qgformula <- as.formula(paste('bmiz', "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
#                               paste( c(exposures.obrien.ln, confounders.wide), collapse = "+")))
# qg.fit.imp.4.og <- list(
#   analyses = lapply(1:imp.length, 
#                     function(x) qgcomp.boot(f = qgformula,
#                                             expnms = exposures.obrien.ln, 
#                                             family = gaussian(), 
#                                             id = 'ifedid',
#                                             q = 4, seed=125, 
#                                             data = dat1[which(dat1$.imp==x),])))
# 
# qg.fit.imp.4.og
# qg.psi.4.og <- MIcombinevals(MIcombine(qg.fit.imp.4.og$analyses),'outcome_mod') # Used exposures.obrien.ln and q=4
# 
# ## Compare all 3 approaches
# qg.psi
# qg.psi.4
# qg.psi.4.og
# ## The first - which is what Alex recommended - is slightly different. So I think this is ok to move forward with. It is frustratingly less precise :/ ##


## Set-up Function: unadjusted models
myqgcompfit <- function(dat, outcome_mod, covariate_mod, exposure_mod){
  
  dat1 <- dat
  dat2 <- dat1[complete.cases(dat1[,c('.imp','ifedid',outcome_mod,'ageinweeks', covariate_mod, exposure_mod)]),] 
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = as.formula(paste(outcome_mod, "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
                                                                   paste(exposure_mod, collapse = "+"))),
                                              expnms = exposure_mod, 
                                              family = gaussian(), 
                                              id = 'ifedid',
                                              q = NULL, bayes = TRUE, seed=125, 
                                              data = dat2[which(dat2$.imp==x),])))
  
  
  #Output parameter estimates for psi1
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),outcome_mod)
  return(qg.psi)
}

u.bmiz   <- myqgcompfit(imp.dat.long, "bmiz", confounders.wide, exposures.obrien.iqr) 
u.bmiz$label <- "Primary Analysis: Unadjusted"
u.bmiz <- data.frame(u.bmiz)



## Adjusted models
myqgcompfit <- function(dat, outcome_mod, covariate_mod, exposure_mod){
  
  dat1 <- dat
  dat2 <- dat1[complete.cases(dat1[,c('.imp','ifedid',outcome_mod,'ageinweeks', covariate_mod, exposure_mod)]),] 
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp(f = as.formula(paste(outcome_mod, "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
                                                              paste(c(exposure_mod, covariate_mod), collapse = "+"))),
                                         expnms = exposure_mod, 
                                         family = gaussian(), 
                                         id = 'ifedid',
                                         q = NULL, bayes = TRUE, seed=125, 
                                         data = dat2[which(dat2$.imp==x),])))
  
  
  #Output parameter estimates for psi1
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),outcome_mod)
  return(qg.psi)
}

a.bmiz   <- myqgcompfit(imp.dat.long, "bmiz", confounders.wide, exposures.obrien.iqr) 
a.bmiz$label <- "Primary Analysis: Adjusted"
a.bmiz <- data.frame(a.bmiz)




## Adjusted models - with interaction term
# Create function that will pull out the pvalue  of the interaction term
MIcombinevals.int <- function(input) {
  MIcombineRes <- MIcombine(input)
  est <- MIcombineRes$coefficients*2
  tStat <- MIcombineRes$coefficients/sqrt(diag(MIcombineRes$variance))
  pval <- 2*pt(-abs(tStat),df=MIcombineRes$df)[[4]]
  pval} 

#Investigate interactions: Sex at all ages
dat <- imp.dat.long[complete.cases(imp.dat.long[,c('.imp','ifedid',"bmiz",'ageinweeks', confounders.wide, exposures.obrien.iqr)]),] 
qgformula <- as.formula(paste("bmiz", "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
                 paste(c(exposures.obrien.iqr, confounders.wide), collapse = "+")))
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.emm.noboot(f = qgformula,
                                                  expnms = exposures.obrien.iqr, 
                                                  family = gaussian(), 
                                                  id = 'ifedid',
                                                  emmvar = 'childsgender',
                                                  q = NULL, bayes = TRUE,
                                                  data = dat[which(dat$.imp==x),])))

qg.fit.imp
MIcombinevals.int(qg.fit.imp$analyses) #0.02

#Investigate interactions: Sex for females before 28 weeks
imp.dat.long$newage <- ifelse(imp.dat.long$childsgender==1 & imp.dat.long$ageinweeks > 28, NA, imp.dat.long$ageinweeks)
summary(imp.dat.long$newage)
f <- subset(imp.dat.long, childsgender==1)
summary(f$newage)

dat1 <- imp.dat.long[complete.cases(imp.dat.long[,c('.imp','ifedid',"bmiz",'newage', confounders.wide, exposures.obrien.iqr)]),] 
qgformula <- as.formula(paste("bmiz", "~", "bs(newage, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
                              paste(c(exposures.obrien.iqr, confounders.wide), collapse = "+")))
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.emm.noboot(f = qgformula,
                                                  expnms = exposures.obrien.iqr, 
                                                  family = gaussian(), 
                                                  id = 'ifedid',
                                                  emmvar = 'childsgender',
                                                  q = NULL, bayes = TRUE,
                                                  data = dat1[which(dat1$.imp==x),])))

qg.fit.imp
MIcombinevals.int(qg.fit.imp$analyses) #0.11

#Investigate interactions: infant diet
dat <- imp.dat.long[complete.cases(imp.dat.long[,c('.imp','ifedid',"bmiz",'ageinweeks', confounders.wide, exposures.obrien.iqr)]),] 
qgformula <- as.formula(paste("bmiz", "~", "bs(ageinweeks, knots=c(3.42857, 7.85714, 15.8571, 24.2857), degree=3) +",
                              paste(c(exposures.obrien.iqr, confounders.wide), collapse = "+")))
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.emm.noboot(f = qgformula,
                                                  expnms = exposures.obrien.iqr, 
                                                  family = gaussian(), 
                                                  id = 'ifedid',
                                                  emmvar = 'feedinggroup',
                                                  q = NULL, bayes = TRUE,
                                                  data = dat[which(dat$.imp==x),])))

qg.fit.imp
MIcombinevals.int(qg.fit.imp$analyses) #0.77




#Stratified model - all females
imp.dat.long.female <- subset(imp.dat.long, imp.dat.long$childsgender==1)
f.bmiz   <- myqgcompfit(imp.dat.long.female, "bmiz", confounders.wide.sex, exposures.obrien.iqr) 
f.bmiz$label = "Secondary Analysis: Females"
f.bmiz <- data.frame(f.bmiz)

## Sensitivity analysis for females limited to <= 28 weeks to mimic the male distributions ##
imp.dat.long.female <- subset(imp.dat.long, imp.dat.long$childsgender==1 & imp.dat.long$ageinweeks<=28)
f2.bmiz   <- myqgcompfit(imp.dat.long.female, "bmiz", confounders.wide.sex, exposures.obrien.iqr) 
f2.bmiz$label = "Sensitivity Analysis 28 Weeks: Females"
f2.bmiz <- data.frame(f2.bmiz)
write.csv(f2.bmiz, "Longitudinal QGCOMP Estimates for Females Before 28 Weeks.csv")

#Stratified model - all males
imp.dat.long.male <- subset(imp.dat.long, imp.dat.long$childsgender==2)
m.bmiz   <- myqgcompfit(imp.dat.long.male, "bmiz", confounders.wide.sex, exposures.obrien.iqr) 
m.bmiz$label = "Secondary Analysis: Males"
m.bmiz <- data.frame(m.bmiz)

#Stratified model - breastmilk
imp.dat.long.breast <- subset(imp.dat.long, imp.dat.long$feedinggroup==1)
b.bmiz   <- myqgcompfit(imp.dat.long.breast, "bmiz", confounders.wide.diet, exposures.obrien.iqr) 
b.bmiz$label = "Secondary Analysis: Breastmilk"
b.bmiz <- data.frame(b.bmiz)

#Stratified model - cowsmilk
imp.dat.long.cow <- subset(imp.dat.long, imp.dat.long$feedinggroup==2)
c.bmiz   <- myqgcompfit(imp.dat.long.cow, "bmiz", confounders.wide.diet, exposures.obrien.iqr) 
c.bmiz$label = "Secondary Analysis: Cows-milk"
c.bmiz <- data.frame(c.bmiz)

#Models for weight-for-age
waz   <- myqgcompfit(imp.dat.long, "waz", confounders.wide, exposures.obrien.iqr) 
waz$label = "Sensitivity Analysis: Weight-for-Age Z-score"
waz <- data.frame(waz)

#Models for visit 6-8
v68.bmiz   <- myqgcompfit(imp.dat.long, "bmiz", confounders.wide, exposures.obrien.iqr.v68) 
v68.bmiz$label <- "Sensitivity Analysis: Exposure at 6-8Weeks"
v68.bmiz <- data.frame(v68.bmiz)

#Models for visit 12
v12.bmiz   <- myqgcompfit(imp.dat.long, "bmiz", confounders.wide, exposures.obrien.iqr.v12) 
v12.bmiz$label <- "Sensitivity Analysis: Exposure at 12Weeks"
v12.bmiz <- data.frame(v12.bmiz)

#Models for non-creatinine adjusted exposure
d.bmiz   <- myqgcompfit(imp.dat.long, "bmiz", confounders.wide, exposures.gm.ln) 
d.bmiz$label <- "Sensitivity Analysis: Non-Creatinine Adjusted Exposure"
d.bmiz <- data.frame(d.bmiz)


## Sensitivity analysis for 12+ weeks
myqgcompfit <- function(dat, outcome_mod, covariate_mod, exposure_mod){
  
  dat1 <- dat
  dat2 <- dat1[complete.cases(dat1[,c('.imp','ifedid',outcome_mod,'ageinweeks', covariate_mod, exposure_mod)]),] 
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp(f = as.formula(paste(outcome_mod, "~", "bs(ageinweeks, knots=c(15.11429, 19.51429, 24, 28.42857), degree=3) +",
                                                              paste(c(exposure_mod, covariate_mod), collapse = "+"))),
                                         expnms = exposure_mod, 
                                         family = gaussian(), 
                                         id = 'ifedid',
                                         q = NULL, bayes = TRUE, seed=125, 
                                         data = dat2[which(dat2$.imp==x),])))
  
  
  #Output parameter estimates for psi1
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),outcome_mod)
  return(qg.psi)
}

imp.dat.long.limit <- subset(imp.dat.long, ageinweeks >= 12)
limit.bmiz   <- myqgcompfit(imp.dat.long.limit, "bmiz", confounders.wide, exposures.obrien.iqr) 
limit.bmiz$label <- "Sensitivity Analysis: 12+ Weeks"
limit.bmiz <- data.frame(limit.bmiz)


## Combine results ##
qgcomp.results.1 = rbind(u.bmiz,a.bmiz,f.bmiz,f2.bmiz,m.bmiz,b.bmiz,c.bmiz,waz,limit.bmiz,v68.bmiz,v12.bmiz,d.bmiz)
qgcomp.results.2 = t(qgcomp.results.1)[c("label","estci2"),] 
write.csv(qgcomp.results.1, "Longitudinal QGCOMP Estimates.csv")
write.csv(qgcomp.results.2, "Table of Longitudinal QGCOMP Estimates.csv")


