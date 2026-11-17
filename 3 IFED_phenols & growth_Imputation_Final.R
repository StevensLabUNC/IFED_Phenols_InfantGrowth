# Programmer: Danielle Stevens
# Last updated on: 11/11/2026
# 
# Purpose of code: This code reads in the IFED exposure, covariate, and confounder data,
#                  and imputes values below LOD for exposurs and creatinine using mice and qgcomp packages
#                  as well as missingness in the batch and smoking variables.
#                  Exposure variables are then cleaned and the O'Brien method of standardization is employed to correct for urinary dilution                            
#                  


############
## Set up ##
############

## Libraries
library(stringr)
# library(tidyverse)
# library(plyr)
# library(tidyr)
# library(dplyr) 
library(mice)
library(miceadds)
library(splines)
library(qgcomp)
library(mde)
library(haven)
library(tidyverse)
# packageVersion('qgcomp') #2.9.0
# packageVersion('mice') #3.14.0


#set working directory
library(here)
datdir <- "J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data"

## Define variables of interest
exposures <- c('bpa', 'bp_3', 'dcp_24', 'dcp_25', 'm_pb', 'p_pb', 'tcs')
exposures.v68 <- paste0(exposures,"_v68")
exposures.v12 <- paste0(exposures,"_v12")
covariates <- c('season', 'ampm', 'ageinweeks', 'batch')
covariates.v68 <- paste0(covariates,"_v68")
covariates.v12 <- paste0(covariates,"_v12")
confounders <- c('feedinggroup', 'mage', 'meduc', 'childsgender', 'crace', 'gestationalage', 'smkpreg', 'alcpregany')


## Create a dataset for the imputation ##
exposures.dat <- read_sas(file.path(datdir, "exposures.sas7bdat"))
confounders.dat <- read_sas(file.path(datdir, "ifed_wide.sas7bdat"))
names(exposures.dat) <- tolower(names(exposures.dat))
names(confounders.dat) <- tolower(names(confounders.dat))
conf2 <- confounders.dat[,c('ifedid', confounders)]

imp.dat <- Reduce(function(x, y) merge(x, y, all=TRUE), 
                  list(confounders.dat[,c('ifedid', confounders)],exposures.dat)) 
# okay
# %>%
  # mutate(type1 = 1)
# imp.dat2 <- full_join(exposures.dat, conf2)%>%
#   mutate(type2 = 1)
# impall <- full_join(imp.dat, imp.dat2)

#############
# Exposures #
#############

# vector of LODs
# Add a 1 for creatinine as well
exposures <- c(exposures,'creat')
lod <- rep(c(0.1,0.2,0.1,0.1,1,0.1,1,1),2)
exposures.list <- c(exposures.v68,"creat_v68",exposures.v12,"creat_v12")
exposures.flag <- paste0(exposures.list, '_flag')

myflag0 <- exposures.dat[, c("ifedid", exposures.list)]
myflag <- pivot_longer(myflag0, -ifedid) %>%
  mutate(name2 = str_replace(name, "p_", "p"),
         name3 = str_replace(name, "m_", "m"),
         nameF = ifelse(substr(name, 1, 4) %in% c("m_pb"), name3, name2)) %>%
  separate(nameF, into = c("name1", "time"), sep = "_") %>%
  dplyr::filter(name1 != "creat") %>%
  mutate(na1 = ifelse(is.na(value), 1, 0)) %>%  
  group_by(ifedid, time) %>%
  summarize(sum1 = sum(na1)) %>%
  pivot_wider(names_from = "time", values_from = "sum1") %>%
  mutate(myflag = ifelse(v68 == 7, 1, ifelse(v12 == 7, 2, 0))) 

# Flag if all missing
imp.dat$flag_exp <- ifelse(is.na(imp.dat$bpa_v68)== TRUE & 
                             is.na(imp.dat$bp_3_v68)== TRUE &  
                             is.na(imp.dat$dcp_24_v68)== TRUE &  
                             is.na(imp.dat$dcp_25_v68)== TRUE &  
                             is.na(imp.dat$m_pb_v68)== TRUE &  
                             is.na(imp.dat$p_pb_v68)== TRUE &  
                             is.na(imp.dat$tcs_v68), 1, 
                           ifelse(is.na(imp.dat$bpa_v12)== TRUE &  
                                    is.na(imp.dat$bp_3_v12)== TRUE &  
                                    is.na(imp.dat$dcp_24_v12)== TRUE &  
                                    is.na(imp.dat$dcp_25_v12)== TRUE &  
                                    is.na(imp.dat$m_pb_v12)== TRUE &  
                                    is.na(imp.dat$p_pb_v12)== TRUE &  
                                    is.na(imp.dat$tcs_v12), 2, 0))

# add my missing
# imp.dat <- full_join(imp.dat, myflag)

# Flag if missing due to <LOD 
table(imp.dat$flag_exp)

exposures.list.v68 <- c(exposures.v68,"creat_v68")
exposures.list.v12 <- c(exposures.v12,"creat_v12")
exposures.flag.v68 <- paste0(exposures.list.v68, '_flag')
exposures.flag.v12 <- paste0(exposures.list.v12, '_flag')

# Create indicator of missing exposure & flag based on type of missingness, 0 = not missing, 1 = <LOD, 2 = missing/not measured
# 1 is all missing 6/8, 2 is all missing 12 from earlier
for(i in 1:length(exposures.list.v68)){
  imp.dat[exposures.flag.v68[i]]<-ifelse(is.na(imp.dat[exposures.list.v68[i]]) & (imp.dat$flag_exp == 0 | imp.dat$flag_exp == 2), 1, 
                                         ifelse(imp.dat$flag_exp == 1, 2, 0))
  print(table(imp.dat[exposures.flag.v68[i]], exclude=F))
}

for(i in 1:length(exposures.list.v12)){
  imp.dat[exposures.flag.v12[i]]<-ifelse(is.na(imp.dat[exposures.list.v12[i]]) & (imp.dat$flag_exp == 0 | imp.dat$flag_exp == 1), 1, 
                                         ifelse(imp.dat$flag_exp == 2, 2, 0))
  print(table(imp.dat[exposures.flag.v12[i]], exclude=F))
}


############################################
## Order & clean variables for imputation ##
############################################
# Data cleaning - get rid of formats and ensure it is a dataframe
imp.dat[] <- lapply(imp.dat, as.vector)
imp.dat <- data.frame(imp.dat)
str(imp.dat)

# Variable formats
imp.dat <- imp.dat %>% mutate(ifedid = as.integer(ifedid),
                              season_v68 = as.factor(season_v68),
                              season_v12 = as.factor(season_v12),
                              ampm_v68 = as.factor(ampm_v68),
                              ampm_v12 = as.factor(ampm_v12),
                              batch_v68 = as.factor(batch_v68),
                              batch_v12 = as.factor(batch_v12),
                              feedinggroup = as.factor(feedinggroup),
                              mage = as.factor(mage),
                              meduc = as.factor(meduc),
                              childsgender = as.factor(childsgender),
                              crace = as.factor(crace),
                              gestationalage = as.factor(gestationalage),
                              smkpreg = as.factor(smkpreg),
                              alcpregany = as.factor(alcpregany))
str(imp.dat)

# Sort by missingness
sortmiss<-sort_by_missingness(imp.dat[,c(confounders, exposures.list, covariates.v68, covariates.v12)],sort_by='counts',descending=FALSE)
miss.order <- c('ifedid', 'flag_exp', exposures.flag, paste(sortmiss$variable))

imp.dat[,miss.order]
df.mice <- imp.dat %>% select(all_of(c(miss.order)))
head(names(df.mice), n=60)


## Functions -------------------------------------------------------------

## Note: These functions are just used to help check changes to data 
#        later in the code.

#### Summary by list elements
simple.summary <- function(mydf, mylist){
  # blank matrix 
  mat1 <- matrix(NA, 
                 nrow = length(mylist), ncol = 4,
                 dimnames = list(mylist, c('min','med','max','NA'))) # row/col names
  # loop through list elements #
  for(i in 1:length(mylist)){
    
    mymin <- min(mydf[[mylist[i]]], na.rm = T)
    mymed <- quantile(mydf[[mylist[i]]],probs = .5, na.rm = T)
    mymax <- max(mydf[[mylist[i]]], na.rm = T)
    myna <- sum(is.na(mydf[[mylist[i]]]))
    
    mat1[i, 1] <- round(mymin, 3)
    mat1[i, 2] <- round(mymed, 3)
    mat1[i, 3] <- round(mymax, 3)
    mat1[i, 4] <- myna
  }
  # output matrix #
  mat1
}

#### TRUE/FALSE TABLE SUMMARY (Used for summary of df.where)
table.truefalse <- function(mydf, mylist){
  # blank matrix 
  mat1 <- matrix(NA, 
                 nrow = length(mylist), ncol = 2,
                 dimnames = list(mylist, c('FALSE','TRUE'))) # row/col names
  # loop through list elements #
  for(i in 1:length(mylist)){
    # save values
    my.false <- table(df.where[mylist[i]])['FALSE']
    my.true <- table(df.where[mylist[i]])['TRUE']
    # input into matrix
    mat1[i, 1] <- my.false
    mat1[i, 2] <- my.true
  }
  # output matrix #
  mat1
}

################
# Where matrix #
################

# TRUE = impute; FALSE = don't impute
# Only impute when exposures.flag == 1 b/c == 2 means it is missing for everyone at that timepoint
df.where <- df.mice 

## Exposure 
for(i in 1:length(exposures.list)){
  df.where[[exposures.list[i]]] <- FALSE # all false to start 
  # Below LOD
  df.where[[exposures.list[i]]][which(df.where[[exposures.flag[i]]] == 1)] <- TRUE}

# check total to be imputed
table.truefalse(df.where, exposures.list) # TRUE = will be imputed (<LOD_impute + missing)

## Covariates: impute 1 NA smoking
# variables not collected longitudinally alongside exposures
for(i in 1:length(confounders)){
  # Make true/false based on NA #

  df.where[[confounders[i]]] <- ifelse(test = is.na(df.where[confounders[i]]),
                                       yes = "TRUE",
                                       no = "FALSE")
}
table.truefalse(df.where, confounders) # NA = no missing values 


# variables  collected longitudinally alongside exposures
for(i in 1:length(covariates.v68)){
  # Make true/false based on NA #
  df.where[[covariates.v68[i]]] <- ifelse(test = is.na(df.where[covariates.v68[i]]),
                                          yes = "TRUE",
                                          no = "FALSE")
  # then edit based on missing that study visit #
  df.where[[covariates.v68[i]]][which(df.where[['flag_exp']] == 1)] <- FALSE
}

table.truefalse(df.where, covariates.v68) # NA = no missing values 

for(i in 1:length(covariates.v12)){
  # Make true/false based on NA #
  df.where[[covariates.v12[i]]] <- ifelse(test = is.na(df.where[covariates.v12[i]]),
                                          yes = "TRUE",
                                          no = "FALSE")
  # then edit based on missing that study visit #
  df.where[[covariates.v12[i]]][which(df.where[['flag_exp']] == 2)] <- FALSE
}

table.truefalse(df.where, covariates.v12) # NA = no missing values 



## Make all other locations = FALSE (ie don't impute)
# select by name
nonimpute <- df.where %>% 
  select(-all_of(c(exposures.list,confounders,covariates.v68,covariates.v12))) %>% 
  names() 
# change
for(i in 1:length(nonimpute)){
  df.where[[nonimpute[i]]] <- FALSE
}
# Check changes # 
table.truefalse(df.where, nonimpute) # everything should be FALSE

df.where #Impute values where TRUE and do not impute where FALSE


##################
# Methods matrix #
##################

## BLANK METHOD MATRIX 
df.method <- matrix("", nrow = 1, ncol = ncol(df.mice), # same cols as original
                    dimnames = list("", names(df.mice))) 

## EXPOSURES - using the left censored log-normal distribution
for(i in 1:length(exposures.list)){
  df.method[[1, exposures.list[i]]] <- "leftcenslognorm" 
}  

## Covariates
cov2 <- c(confounders,covariates.v68,covariates.v12)
for(i in 1:length(cov2)){
  # indicator for binary vs multinomial factor # 
  factorlength <- length(levels(df.mice[[cov2[i]]]))
  # input method based on numeric, binary, or multinomial #
  df.method[[1, cov2[i]]] <- ifelse(
    test = is.numeric(df.mice[[cov2[i]]]),
    yes = "pmm", # numeric
    no = ifelse(test = factorlength==2,
                yes = "logreg", # binary factor:input specific function below
                no = "polyreg")) # categorical factor:input specific function below
  rm(factorlength)
}

# check methods 
df.method



#################################################
# Left-censoring matrix of values for exposures #
#################################################

# For <LOD variables, use LOD values. 
# For NA exposure variables, use larger range because it helps reduce model assumptions (see below)
df.lod.values <- df.mice 

### CENSOR POINTS for NA Exposures 
# Per Alex Keil, can use large value censor points for missing exposures due to machine error/inadequate sample
# Censor at data max or 95th percentile values?
range(df.lod.values[[exposures.list[1]]], na.rm = T)
quantile(df.lod.values[[exposures.list[1]]], probs = .95, na.rm = T) # 95th percentile 
# Decision: 95th seems safer than maximum #

### ALTER LOD.DF FOR IMPUTATION
exposure.lod <- paste0(exposures.list, '_lod')
for(i in 1:length(exposures.list)){
  # INPUT LEFT-CENSOR POINT FOR MISSING EXPOSURE CONCENTRATIONS 
  # set everything to LOD
  df.lod.values[exposure.lod[i]]<-lod[i]
  # for those that you need to impute: not imputing below the LOD
  df.lod.values[which(df.lod.values[[exposures.flag[i]]]==1), exposure.lod[i]] <- 
    as.numeric(quantile(df.lod.values[[exposures.list[i]]], probs = .05, na.rm = T)) # 95th percentile 
  # replace normal variable (same name needed, so replace normal var with lod var) #
  df.lod.values[[exposures.list[i]]] <- df.lod.values[[exposure.lod[i]]] 
}

# CHECK changes (exposure.lod and exposure should match)
table(df.lod.values[[exposure.lod[1]]], exclude = T)
table(df.lod.values[[exposures.list[1]]], exclude = T)



#####################
# Predictors matrix #
#####################
# Set up matrix
df.predict = NA

# Include: Variables to use in every imputation model
includelist <- c(cov2, exposures.list) 
includelist 

# Exclude: Variables to not use in any imputations (all others besides above
excludelist <- df.mice %>% select(-all_of(includelist)) %>% names() 

# Specify matrix 
df.predict <- quickpred(df.mice, # automatic selection
                        include = includelist, # 1 = predictors (won't be used for self-prediction) 
                        exclude = excludelist)  # 0 = not used as predictors 

# Update predictor matrix to include visit-specific gest'l age and creatinine for exposures
df.predict[exposures.list.v12,covariates.v68]<- 0 #Not predicting exposures at 12 weeks from batch and timing of exposure assessment at 68 weeks
df.predict[exposures.list.v68,covariates.v12]<- 0 #Not predicting exposures at 6/8 weeks from batch and timing of exposure assessment at 12 weeks
df.predict[covariates.v68,exposures.list.v12]<- 0 #Not predicting covariates at 6/8 weeks from exposures at 12 weeks
df.predict[covariates.v12,exposures.list.v68]<- 0 #Not predicting covariates at 12 weeks from exposures at 68 weeks
# JRK FIX: too many missing predictors to use/causes missingness
df.predict[exposures.list.v12,exposures.list.v68]<- 0 #Not predicting exposures at 6/8 weeks from exposures at 68 weeks
df.predict[exposures.list.v68,exposures.list.v12]<- 0 #Not predicting exposures at 12 weeks from exposures at 12 weeks


## NOTE: I am allowing exposures at other timepoints in the model ##
df.predict


##############
## Run MICE ##
##############
my.seed <- 3456 # seed

my.seed <- 10598
start_time <- Sys.time()
mids.original <- mice(df.mice, 
                      m = 10, # Num. of imputed datasets
                      maxit = 20, # Num. of iterations per imputed dataset 
                      seed = my.seed,
                      donors = 10L,
                      # METHODS # 
                      method = df.method, 
                      # JRK FIX: df.lod.values above used 95%ile, imputing
                      # values across the range of concentrations
                      # leaving NULL uses smallest value observed in data
                    #  lod = df.lod.values[, 1 : 50], # left-censored(lod) points for phthalates
                      # WHERE/Predictors #
                      where = df.where, # what is imputed
                      predictorMatrix = df.predict, # what is used for imputation
                      # preventing errors # 
                      remove_collinear = F, # ensures dropping with collinearity
                      remove.constant  = F)

end_time <- Sys.time()
end_time - start_time # Time difference 

warnings()
View(mids.original$loggedEvents) # loggedEvents
mids.original$method # check methods 
mids.original$predictorMatrix
mids.original$formulas$smkpreg # example formula used
mids.original$formulas$batch_v68 

## Plots
library(grid)
pdf(file="Imputation_diagnostics.pdf")
my.plot <- plot(mids.original)
my.plot
dev.off()


#####################
## Clean exposures ##
#####################

## Make into dataframe for editing
imp.ifed = complete(mids.original, "long", include = T)

## O'Brien method for cleaning exposures
## Predict creatinine at each study visit from confounders and visit-specific covariates

# Predictors
creat.predictors.v68 <- c(confounders,covariates.v68)
creat.predictors.v12 <- c(confounders,covariates.v12)

## Run model for creatinine at each timepoint
a<-imp.ifed[,c(exposures.list, creat.predictors.v68, creat.predictors.v12, "ifedid", ".imp")]

## Run a few to examine predictors
imputation = subset(a, .imp == 2)
out = glm(data=imputation, as.formula(paste('creat_v12', "~", paste(creat.predictors.v12, collapse = "+"))))
summary(out) #ampm sign for 68
             # race, education, smoking, and batch sign for 12 

imputation = subset(a, .imp == 5)
out = glm(data=imputation, as.formula(paste('creat_v68', "~", paste(creat.predictors.v68, collapse = "+"))))
summary(out) #am pm sign for 68

imputation = subset(a, .imp == 9)
out = glm(data=imputation, as.formula(paste('creat_v12', "~", paste(creat.predictors.v12, collapse = "+"))))
summary(out) #ampm and batch sign for 68
             # education, race, smoking, and batch sign for 12

## Median exposure at each visit
median(a$creat_v68[which(a$.imp==0)], na.rm = T) #11.9
median(a$creat_v12[which(a$.imp==0)], na.rm = T) #12.65

## Cycle through each visit & imputation
exposures <- c("bpa","bp_3","dcp_24","dcp_25","m_pb","p_pb","tcs")
exposure.boeniger.v68 <- paste0(exposures, "_boeniger_v68")
exposure.boeniger.v12 <- paste0(exposures, "_boeniger_v12")
exposure.obrien.v68 <- paste0(exposures, "_obrien_v68")
exposure.obrien.v12 <- paste0(exposures, "_obrien_v12")

values.v68 = data.frame()
for(j in 0:10){
  imputation = subset(a, .imp == j)
  out = glm(data=imputation, as.formula(paste('creat_v68', "~", paste(creat.predictors.v68, collapse = "+"))))
  imputation$creatpred_v68 = predict(out, newdata = imputation, interval="none") 
  for(i in 1:length(exposures.v68)){
    imputation[exposure.boeniger.v68[i]] <- imputation[exposures.v68[i]]*((11.9)/(imputation$creat_v68))
    imputation[exposure.obrien.v68[i]] <- imputation[exposures.v68[i]]*((imputation$creatpred_v68)/(imputation$creat_v68))}
  values = imputation %>% select(all_of(c("ifedid", ".imp", "creatpred_v68", exposure.boeniger.v68, exposure.obrien.v68)))
  values.v68 = rbind(values.v68, values)
}

values.v12 = data.frame()
for(j in 0:10){
  imputation = subset(a, .imp == j)
  out = glm(data=imputation, as.formula(paste('creat_v12', "~", paste(creat.predictors.v12, collapse = "+"))))
  imputation$creatpred_v12 = predict(out, newdata = imputation, interval="none") 
  
  # take one negative value and make small positive
  wh1 <- which(imputation$creatpred_v12 < 0)
  imputation$creatpred_v12[wh1] <- 0.01
  
  for(i in 1:length(exposures.v12)){
    imputation[exposure.boeniger.v12[i]] <- imputation[exposures.v12[i]]*((11.9)/(imputation$creat_v12))
    imputation[exposure.obrien.v12[i]] <- imputation[exposures.v12[i]]*((imputation$creatpred_v12)/(imputation$creat_v12))}
  values = imputation %>% select(all_of(c("ifedid", ".imp", "creatpred_v12", exposure.boeniger.v12, exposure.obrien.v12)))
  values.v12 = rbind(values.v12, values)
}

head(values.v12)

df_list <- list(imp.ifed,values.v68,values.v12)
df_list <- df_list %>% reduce(full_join, by = c("ifedid", ".imp"))

summary(df_list$bpa_v68)
summary(df_list$bpa_boeniger_v68)
summary(df_list$bpa_obrien_v68)


## Log-transform
exposure.boeniger.ln.v68 <- paste0(exposure.boeniger.v68, '_ln')
exposure.boeniger.ln.v12 <- paste0(exposure.boeniger.v12, '_ln')
exposure.boeniger.ln <- c(exposure.boeniger.ln.v68, exposure.boeniger.ln.v12)
exposure.boeniger <- c(exposure.boeniger.v68, exposure.boeniger.v12)

exposure.obrien.ln.v68 <- paste0(exposure.obrien.v68, '_ln')
exposure.obrien.ln.v12 <- paste0(exposure.obrien.v12, '_ln')
exposure.obrien.ln <- c(exposure.obrien.ln.v68, exposure.obrien.ln.v12)
exposure.obrien <- c(exposure.obrien.v68, exposure.obrien.v12)

# only 1 negative for obrien 12
for(i in 1:length(exposure.boeniger)){
  print(exposure.boeniger[i])
  #print(c(min(df_list[exposure.boeniger[i]], na.rm = T),min(df_list[exposure.obrien[i]] , na.rm = T)))
  # obrien negative for some 12
  temp <- df_list[exposure.obrien[i]]
  print(temp[which(df_list[exposure.obrien[i]] < 0), ])
  df_list[exposure.boeniger.ln[i]] = log(df_list[exposure.boeniger[i]])
  df_list[exposure.obrien.ln[i]] = log(df_list[exposure.obrien[i]])
}

summary(df_list[exposure.obrien])
summary(df_list[exposure.obrien.ln])

## Geometric mean for pregnancy-averaged exposure, then log-transform
exposure.boeniger.gm <- paste0(exposures, '_boeniger_gm')
exposure.obrien.gm <- paste0(exposures, '_obrien_gm')
exposure.boeniger.gm.ln <- paste0(exposures, '_boeniger_gm_ln')
exposure.obrien.gm.ln <- paste0(exposures, '_obrien_gm_ln')

for(i in 1:length(exposures)){ 
  df_list[exposure.boeniger.gm.ln[i]] = rowMeans(df_list[,c(exposure.boeniger.ln.v68[i],exposure.boeniger.ln.v12[i])], na.rm=TRUE)
  df_list[exposure.obrien.gm.ln[i]] = rowMeans(df_list[,c(exposure.obrien.ln.v68[i],exposure.obrien.ln.v12[i])], na.rm=TRUE)
  
  df_list[exposure.boeniger.gm[i]] = exp(df_list[exposure.obrien.gm.ln[i]])
  df_list[exposure.obrien.gm[i]] = exp(df_list[exposure.obrien.gm.ln[i]])
}

summary(df_list[exposure.obrien.gm])
view(df_list[1:20,c('bpa_v68','bpa_v12','bpa_obrien_v68_ln','bpa_obrien_v12_ln','bpa_obrien_v68','bpa_obrien_v12',
                    'bpa_obrien_gm_ln','bpa_obrien_gm')])


# check new variables
dat2 <- dplyr::select(df_list, ifedid, flag_exp, .imp, creatpred_v68 :tcs_obrien_gm) %>% dplyr::filter(.imp >0)#tcs_obrien_gm)
for(i in 2 : ncol(dat2)) {
  # flag_exp
  dat2var <- dat2[, colnames(dat2)[i]]

  wh1 <- unique(c(which(is.na(dat2var)), which(is.nan(dat2var))))
  dat2var <- dat2[wh1, c("ifedid", "flag_exp", ".imp", colnames(dat2)[i])]
  dat2var <- dat2var %>% distinct() %>% 
    dplyr::filter(flag_exp ==0)
  #  group_by(ifedid, flag_exp) %>% count()
  if(nrow(dat2var) > 0) {# && dat2var[,".imp"] != 1 ) {
    print(colnames(dat2)[i])
    
  print(dat2var)
  }

}



## Fix the few inf and NaNs that this causes
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
# only replacing .imp = 0: not a problem
df_list2 <- replace(df_list, is.nan(df_list), NA)

df_list[df_list==Inf] = NA
df_list[df_list==-Inf] = NA


#write.csv(df_list, "U:\\IFED\\data\\Imputed dataset.csv")
save(df_list, file =file.path(datdir, "df_list.RData"))




#############################################
## Merge outcomes into this dataframe      ##
## Create two datasets:                    ##
## 1) long form dataframe for bmiz and waz ##
## 2) wide form dataframe for GMM results  ##
#############################################

## Imputed exposures and covariates dataframe
imp.dat <- df_list 
summary(imp.dat)
summary(imp.dat$.imp)
table(imp.dat$ifedid) #2189 or 11 rows per participant (0 = original dataframe, 1-10 = imputed dataframe)
11*199 #2189


## Long form outcomes dataframe
outcomes <- read_sas(file.path(datdir,"ifed_long.sas7bdat"))
names(outcomes) <- tolower(names(outcomes))
summary(outcomes)
table(outcomes$ifedid) #1991 observations (max 12 visits/participant)

## Merge
imp.dat.long <- Reduce(function(x, y) merge(x, y, all=TRUE), list(imp.dat, outcomes[,c("ifedid",'ageinweeks',"bmiz","waz")]))
table(imp.dat.long$ifedid) #Looks correct (x 11 from above)


# outcomes2 <- dplyr::select(outcomes, ifedid, ageinweeks, bmiz, waz) %>% unique()
# imp2 <- full_join(outcomes2, imp.dat) %>% 
#   mutate(keep  = 1)
# test <- full_join(imp2, imp.dat.long)
# table(test$keep)

# Check merge occurred correctly
table(imp.dat.long$.imp) #10 imputations with 1991 in each imputation 
view(imp.dat.long[1:100,c('.imp','ifedid','bpa_obrien_gm_ln',confounders,'ageinweeks','bmiz')])
view(imp.dat.long[15000:15100,c('.imp','ifedid','bpa_obrien_gm_ln',confounders,'ageinweeks','bmiz')])
view(imp.dat.long[2000:2100,c('.imp','ifedid','bpa_obrien_gm_ln',confounders,'ageinweeks','bmiz')])

## On 11/15/2023, I calculated the geometric mean of the non-creatinine adjusted measures and log-transformed for a sensitivity analysis
## Log-transform & calculate geometric mean for pregnancy-averaged exposure
exposures.all <- c(exposures.v68, exposures.v12)
exposure.ln <- paste0(exposures.all, "_ln")
exposure.ln.v68 <- paste0(exposures.v68, "_ln")
exposure.ln.v12 <- paste0(exposures.v12, "_ln")
exposure.gm <- paste0(exposures, '_gm')
exposure.gm.ln <- paste0(exposures, '_gm_ln')

for(i in 1:length(exposures.all)){
  imp.dat.long[exposure.ln[i]] = log(imp.dat.long[exposures.all[i]])}
summary(imp.dat.long[exposure.ln])

for(i in 1:length(exposures)){
   imp.dat.long[exposure.gm.ln[i]] = rowMeans(imp.dat.long[,c(exposure.ln.v68[i],exposure.ln.v12[i])], na.rm=TRUE)
   imp.dat.long[exposure.gm[i]] = exp(imp.dat.long[exposure.gm.ln[i]])
}
summary(imp.dat.long[exposure.gm.ln])
summary(imp.dat.long[exposure.gm])

view(imp.dat.long[1:20,c('bpa_v68','bpa_v12','bpa_v68_ln','bpa_v12_ln', 'bpa_gm_ln','bpa_gm')])


## Fix the few inf and NaNs that this causes
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

imp.dat.long <- replace(imp.dat.long, is.nan(imp.dat.long), NA)

imp.dat.long[imp.dat.long==Inf] = NA
imp.dat.long[imp.dat.long==-Inf] = NA

write.csv(imp.dat.long, "J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\Imputed dataset with outcomes long.csv")
save(imp.dat.long, file = file.path(datdir, "implong.RData"))





## Wide form
outcomes <- read.csv("J:\\StevensLab\\IFED Study\\IFED Phenols & Growth_DRS\\Analysis\\Derived data\\GMM_Dataset_Wide.csv")
colnames(outcomes)
table(outcomes$ifedid)


imp.dat.wide <- Reduce(function(x, y) merge(x, y, all=TRUE), list(imp.dat, outcomes[,2:42]))
table(imp.dat.wide$ifedid) #Looks correct (x 11 from above)


# outcomes2 <- outcomes[, 2 : 42] %>% unique()
# imp2 <- full_join(outcomes2, imp.dat) %>%
#   mutate(keep  = 1)
# test <- full_join(imp2, imp.dat)
# table(test$keep)


## Same as above - add the non-creatinine adjusted measures back in
for(i in 1:length(exposures.all)){
  imp.dat.wide[exposure.ln[i]] = log(imp.dat.wide[exposures.all[i]])}
summary(imp.dat.wide[exposure.ln])

for(i in 1:length(exposures)){
  imp.dat.wide[exposure.gm.ln[i]] = rowMeans(imp.dat.wide[,c(exposure.ln.v68[i],exposure.ln.v12[i])], na.rm=TRUE)
  imp.dat.wide[exposure.gm[i]] = exp(imp.dat.wide[exposure.gm.ln[i]])
}
summary(imp.dat.wide[exposure.gm.ln])
summary(imp.dat.wide[exposure.gm])

view(imp.dat.wide[1:20,c('bpa_v68','bpa_v12','bpa_v68_ln','bpa_v12_ln', 'bpa_gm_ln','bpa_gm')])

## Fix the few inf and NaNs that this causes
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

imp.dat.wide <- replace(imp.dat.wide, is.nan(imp.dat.wide), NA)
# no infinite
imp.dat.wide[imp.dat.wide==Inf] = NA
imp.dat.wide[imp.dat.wide==-Inf] = NA


write.csv(imp.dat.wide, file.path(datdir,"Imputed dataset with outcomes wide.csv"))


