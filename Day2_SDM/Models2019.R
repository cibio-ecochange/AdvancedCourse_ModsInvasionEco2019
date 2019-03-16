##############################################
#########  Advanced course BIODIV2019 ########
##############################################


#######################################
############  LOAD BIOMOD #############
#######################################

library(biomod2) 


#######################################
#### IMPORTING AND EXPLORING DATA #####
#######################################

# Set the proper working directory from the base project location
setwd("./Day2_SDM")


### importing and exploring data | Acacia dealbata ###
acadea <- read.table("./DATA/acadea.txt", h=T, sep="\t")

is.data.frame(acadea)
head(acadea)
dim(acadea)
acadea[1:10,1]
acadea[1,1:5]
names(acadea)
ls()


### importing and exploring data | Ruscus aculeatus ###
rusacu <- read.table("./DATA/rusacu.txt", h=T, sep="\t")

is.data.frame(rusacu)
head(rusacu)
dim(rusacu)
rusacu[1:10,1]
rusacu[1,1:5]
names(rusacu)
ls()


#######################################
#########  SEE DATA IN SPACE ##########
#######################################

### Acacia dealbata ###
level.plot(acadea[,"acadea"], acadea [,(2:3)], title= "Acacia dealbata - sampled data")

### Environmental Predictors | Acacia dealbata ###
level.plot(acadea[,(4)], acadea [,(2:3)], title = "Minimum Temperature of Coldest Month")
level.plot(acadea[,(5)], acadea [,(2:3)], title = "Precipitation of Coldest Quarter")
level.plot(acadea[,(6)], acadea [,(2:3)], title = "Distance to the Hydrographical Network")
level.plot(acadea[,(7)], acadea [,(2:3)], title = "Diversity Index of the Slope")
level.plot(acadea[,(8)], acadea [,(2:3)], title = "Percentage of Arenossoils")
level.plot(acadea[,(9)], acadea [,(2:3)], title = "Diversity Index of the Soils")
level.plot(acadea[,(10)], acadea [,(2:3)], title = "Percentage of Artificial Stands")
level.plot(acadea[,(11)], acadea [,(2:3)], title = "Percentage of Agriculture")

save.image()

### Ruscus aculeatus ###
level.plot(rusacu[,"rusacu"], rusacu [,(2:3)], title= "Ruscus aculeatus - sampled data")

### Environmental Predictors | Ruscus aculeatus ###
level.plot(rusacu[,(4)], rusacu [,(2:3)], title = "Minimum Temperature of Coldest Month")
level.plot(rusacu[,(5)], rusacu [,(2:3)], title = "Precipitation of Wettest Quarter")
level.plot(rusacu[,(6)], rusacu [,(2:3)], title = "Distance to Urban Areas")
level.plot(rusacu[,(7)], rusacu [,(2:3)], title = "Diversity Index of the Slope")
level.plot(rusacu[,(8)], rusacu [,(2:3)], title = "Percentage of Regossoils")
level.plot(rusacu[,(9)], rusacu [,(2:3)], title = "Percentage of Semi-Natural Habitats")
level.plot(rusacu[,(10)], rusacu [,(2:3)], title = "Percentage of Artificial Stands")
level.plot(rusacu[,(11)], rusacu [,(2:3)], title = "Hidrographical Network Density")

save.image()


#######################################
######### Acacia dealbata  ############
#######################################

## get data
acadea <- read.table ("./DATA/acadea.txt", h=T, sep="\t")


## separating the response data, coordinates and predictors
myResp <- as.numeric(acadea$acadea)
myRespCoord <- acadea [,c("X", "Y")] 
myRespName <- "acadea"
myExpl <- acadea [, 4:11]

# formating the data, without pseudoabsences and without independet evaluation data
myBiomodData_Ad <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     eval.resp.var = NULL,
                                     eval.expl.var = NULL,
                                     eval.resp.xy = NULL,
                                     PA.nb.rep = 0,
                                     PA.nb.absences = 0,
                                     PA.strategy = 'random',
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     na.rm = TRUE)


# check the formatted data
myBiomodData_Ad

# plott it
plot(myBiomodData_Ad)

### Models' Calibration ###
myBiomodModelOut_Ad <- BIOMOD_Modeling( 
                           myBiomodData_Ad, 
                           models = c('GLM','GAM','GBM','CTA','ANN','SRE',
                                      'FDA','MARS','RF', "MAXENT.Tsuruoka"), 
                           NbRunEval = 3, 
                           DataSplit = 80, 
                           Yweights = NULL, 
                           VarImport = 5, 
                           models.eval.meth = c('ROC'),
                           SaveObj = TRUE,
                           rescal.all.models = TRUE,
                           models.options = BIOMOD_ModelingOptions(GAM=list(k=2)))

## summary of the model results
myBiomodModelOut_Ad


## model evaluations, assigned to new variable
myBiomodModelEval_Ad <- get_evaluations(myBiomodModelOut_Ad)

## models' ROC scores
myBiomodModelEval_Ad["ROC","Testing.data",,,]

# checking variables' importance 
get_variables_importance(myBiomodModelOut_Ad)


#######################################################################
## ENSEMBLE MODELLING
## build consensus model

myBiomodEM_Ad <- BIOMOD_EnsembleModeling( 
                     modeling.output = myBiomodModelOut_Ad,
                     chosen.models = 'all', 
                     em.by = 'all',
                     eval.metric = c('ROC'),
                     eval.metric.quality.threshold = c(0.7),
                     prob.mean = TRUE,
                     prob.cv = TRUE,
                     prob.ci = TRUE,
                     prob.ci.alpha = 0.05,
                     prob.median = TRUE,
                     committee.averaging = TRUE,
                     prob.mean.weight = TRUE,
                     prob.mean.weight.decay = 'proportional' )

# print summary                     
myBiomodEM_Ad
                     
# get evaluation scores 
get_evaluations(myBiomodEM_Ad)

# getting model's variable importances
## we need to load ensemble models
ensemble_models_names <- BIOMOD_LoadModels(myBiomodEM_Ad)
ensemble_models_names

varImpWMean <- (variables_importance(model=get("acadea_EMwmeanByROC_mergedAlgo_mergedRun_mergedData" ), data=myExpl, method="full_rand", nb_rand=2))

varImpWMean

########################################################
###  Model's Spatial Projection | Current conditions - 2000 ###
acadea_2000 <- read.table ("./DATA/spatAd.txt", h=T, sep="\t")

head(acadea_2000)
dim(acadea_2000)

myNewRespCoord_Ad2000 <- acadea_2000[,c("X", "Y")]
myNewExpla_Ad2000 <- acadea_2000[,4:11]

myBiomodProj_Ad2000 <- BIOMOD_Projection(
                         modeling.output = myBiomodModelOut_Ad,
                         new.env = myNewExpla_Ad2000,
                         proj.name = 'acadea_2000',
                         xy.new.env= myNewRespCoord_Ad2000,
                         selected.models = 'all',
                         binary.meth = 'ROC',
                         filtered.meth = 'ROC',
                         compress = 'xz',
                         clamping.mask = FALSE)

####Plot Acacia dealbata current models
plot(myBiomodProj_Ad2000)


#ENSEMBLE FORECAST - Acacia dealbato 2000
myBiomodEF_Ad2000 <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProj_Ad2000,
                            EM.output = myBiomodEM_Ad,
                            total.consensus = TRUE,
                            binary.meth = 'ROC',
                            filtered.meth = 'ROC' )

plot(myBiomodEF_Ad2000)

# loading and saving results

load("acadea/proj_acadea_2000/proj_acadea_2000_acadea_ensemble.RData")

load("acadea/proj_acadea_2000/proj_acadea_2000_acadea_ensemble_ROCbin.RData")


########################################################
###  Model's Spatial Projection | Future conditions - 2020 ###
acadea_2020 <- read.table ("./DATA/spat2020Ad.txt", h=T, sep="\t")

head(acadea_2020)

myNewRespCoord_Ad2020 <- acadea_2020[,c("X", "Y")]
myNewExpla_Ad2020 <- acadea_2020[,4:11]

myBiomodProj_Ad2020 <- BIOMOD_Projection(
                         modeling.output = myBiomodModelOut_Ad,
                         new.env = myNewExpla_Ad2020,
                         proj.name = 'acadea_2020',
                         xy.new.env= myNewRespCoord_Ad2020,
                         selected.models = 'all',
                         binary.meth = 'ROC',
                         filtered.meth = 'ROC',
                         compress = 'xz',
                         clamping.mask = FALSE)

####Plot Acacia dealbata future models
plot(myBiomodProj_Ad2020)


#ENSEMBLE FORECAST - Acacia dealbato 2020
myBiomodEF_Ad2020 <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProj_Ad2020,
                            EM.output = myBiomodEM_Ad,
                            total.consensus = TRUE,
                            binary.meth = 'ROC',
                            filtered.meth = 'ROC' )

####Plot Acacia dealbata future ensemble forecasted models
plot(myBiomodEF_Ad2020)


# loading and saving results

load("acadea/proj_acadea_2020/proj_acadea_2020_acadea_ensemble.RData")
load("acadea/proj_acadea_2020/proj_acadea_2020_acadea_ensemble_ROCbin.RData")



#######################################
######### Ruscus aculeatus  ###########
#######################################

## get data
rusacu <- read.table ("./DATA/rusacu.txt", h=T, sep="\t")


## separating the response data, coordinates and predictors
myResp <- as.numeric(rusacu$rusacu)
myRespCoord <- rusacu [,c("X", "Y")] 
myRespName <- "rusacu"
myExpl <- rusacu [, 4:11]

# formating the data, without pseudoabsences and without independet evaluation data
myBiomodData_Ra <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     eval.resp.var = NULL,
                                     eval.expl.var = NULL,
                                     eval.resp.xy = NULL,
                                     PA.nb.rep = 0,
                                     PA.nb.absences = 0,
                                     PA.strategy = 'random',
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     na.rm = TRUE)


# check the formatted data
myBiomodData_Ra

# plott it
plot(myBiomodData_Ra)

### Models' Calibration ###
myBiomodModelOut_Ra <- BIOMOD_Modeling( 
                           myBiomodData_Ra, 
                           models = c('GLM','GAM','GBM','CTA','ANN','SRE','FDA',
                                      'MARS','RF', "MAXENT.Tsuruoka"), 
                           models.options = NULL, 
                           NbRunEval = 3, 
                           DataSplit = 80, 
                           Yweights = NULL, 
                           VarImport = 5, 
                           models.eval.meth = c('ROC'),
                           SaveObj = TRUE,
                           rescal.all.models = TRUE)

## summary of the model results
myBiomodModelOut_Ra

## model evaluations, assigned to new variable
myBiomodModelEval_Ra <- get_evaluations(myBiomodModelOut_Ra)

## models' ROC scores
myBiomodModelEval_Ra["ROC","Testing.data",,,]

# checking variables' importance 
get_variables_importance(myBiomodModelOut_Ra)

#######################################################################
## ENSEMBLE MODELLING
## build consensus model

myBiomodEM_Ra <- BIOMOD_EnsembleModeling( 
                     modeling.output = myBiomodModelOut_Ra,
                     chosen.models = 'all', 
                     em.by = 'all',
                     eval.metric = c('ROC'),
                     eval.metric.quality.threshold = c(0.7),
                     prob.mean = TRUE,
                     prob.cv = TRUE,
                     prob.ci = TRUE,
                     prob.ci.alpha = 0.05,
                     prob.median = TRUE,
                     committee.averaging = TRUE,
                     prob.mean.weight = TRUE,
                     prob.mean.weight.decay = 'proportional' )

# print summary                     
myBiomodEM_Ra
                     
# get evaluation scores 
get_evaluations(myBiomodEM_Ra)

# getting model's variable importances
## we need to load ensemble models
ensemble_models_names_Ra <- BIOMOD_LoadModels(myBiomodEM_Ra)
ensemble_models_names_Ra

varImpWMean <- (variables_importance(model=get("rusacu_EMwmeanByROC_mergedAlgo_mergedRun_mergedData" ), data=myExpl, method="full_rand", nb_rand=2))


varImpWMean

########################################################
###  Model's Spatial Projection | Current conditions - 2000 ###
rusacu_2000 <- read.table ("./DATA/spatRa.txt", h=T, sep="\t")

myNewRespCoord_Ra2000 <- rusacu_2000[,c("X", "Y")]
myNewExpla_Ra2000 <- rusacu_2000[,4:11]

myBiomodProj_Ra2000 <- BIOMOD_Projection(
                         modeling.output = myBiomodModelOut_Ra,
                         new.env = myNewExpla_Ra2000,
                         proj.name = 'rusacu_2000',
                         xy.new.env= myNewRespCoord_Ra2000,
                         selected.models = 'all',
                         binary.meth = 'ROC',
                         filtered.meth = 'ROC',
                         compress = 'xz',
                         clamping.mask = FALSE)
                         
                         
####Plot Ruscus aculeatus current models
#plot(myBiomodProj_Ra2000)
                         

#ENSEMBLE FORECAST - Ruscus aculeatus 2000
myBiomodEF_Ra2000 <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProj_Ra2000,
                            EM.output = myBiomodEM_Ra,
                            total.consensus = TRUE,
                            binary.meth = 'ROC',
                            filtered.meth = 'ROC' )


####Plot Ruscus aculeatus current ensemble forecasted models
#plot(myBiomodEF_Ra2000)

# loading and saving results

load("rusacu/proj_rusacu_2000/proj_rusacu_2000_rusacu_ensemble.RData")

load("rusacu/proj_rusacu_2000/proj_rusacu_2000_rusacu_ensemble_ROCbin.RData")



########################################################
###  Model's Spatial Projection | Future conditions - 2020 ###
rusacu_2020 <- read.table ("./DATA/spat2020Ra.txt", h=T, sep="\t")

myNewRespCoord_Ra2020 <- rusacu_2020[,c("X", "Y")]
myNewExpla_Ra2020 <- rusacu_2020[,4:11]

myBiomodProj_Ra2020 <- BIOMOD_Projection(
                         modeling.output = myBiomodModelOut_Ra,
                         new.env = myNewExpla_Ra2020,
                         proj.name = 'rusacu_2020',
                         xy.new.env= myNewRespCoord_Ra2020,
                         selected.models = 'all',
                         binary.meth = 'ROC',
                         filtered.meth = 'ROC',
                         compress = 'xz',
                         clamping.mask = FALSE)

####Plot Ruscus aculeatus future models
#plot(myBiomodProj_Ra2020)

#ENSEMBLE FORECAST - Ruscus aculeatus 2020
myBiomodEF_Ra2020 <- BIOMOD_EnsembleForecasting( projection.output = myBiomodProj_Ra2020,
                            EM.output = myBiomodEM_Ra,
                            total.consensus = TRUE,
                            binary.meth = 'ROC',
                            filtered.meth = 'ROC' )

####Plot Ruscus aculeatus future ensemble forecasted models
#plot(myBiomodEF_Ra2020)


# loading and saving results

load("rusacu/proj_rusacu_2020/proj_rusacu_2020_rusacu_ensemble.RData")
load("rusacu/proj_rusacu_2020/proj_rusacu_2020_rusacu_ensemble_ROCbin.RData")


#####################################
####### SPECIES RANGE CHANGE ########
#####################################

###  Compt.By.Species ###

#Stores the summary of range change for each species (sorted by rows). The first four columns are absolute values whereas the next 3 ones are relative values: 
#Disa represents the number of pixels predicted to be lost by the given species. 
#Stable0 is the number of pixels which are not currently occupied by the given species 
#and not predicted to be. 

#Stable1 represents the number of pixels currently occupied by the given species, and predicted to remain occupied into the future. 

#Gain represents the number of pixels which are currently not occupied by the given species but predicted to be into the future. 
#PercLoss corresponds to the percentage of currently occupied sites to be lost (Disa/(Disa+Stable1) 
#PercGain corresponds to the percentage of new sites considering the species??? current dis- tribution size (Gain/(Disa+Stable1). For example, if the there are 30 sites cur- rently occupied and 15 new sites are projected to be occupied in future, it makes PercGain=+50(%). 
#SpeciesRangeChange
#it is the overall projection outcome, equal to PercGain-PercLoss. It does not assess for any migration shifts as it strictly compares the range sizes between current and future states. 
#CurrentRangeSize
#represents the modelled current range size (number of pixels occupied) of the given species. 
#FutureRangeSize0Disp
#represents the future modelled range size assuming no migration of the given species. 
#FutureRangeSize1Disp
#Diff.By.Pixel#
#represents the future modelled range size assuming migration of the given species (depending on the datasets given in input, if Migration has been used or not). 
#the summary of range change for each species (sorted by columns and with the pixel in rows). For each species, a pixel could have four different values : -2 if the given pixel is predicted to be lost by the species. -1 if the given pixel is predicted to be stable for the species. 0 is the given pixel was not occupied, and will not be in the future. 1 if the given pixel was not occupied, and is predicted to be into the future. This table could be easily plotted into GIS software in order to represent the pattern of change for the selected species (or even with the level.plot() function). ####


### Acacia dealbata ###
load("acadea/proj_acadea_2020/proj_acadea_2020_acadea_ensemble_ROCbin.RData")
load("acadea/proj_acadea_2000/proj_acadea_2000_acadea_ensemble_ROCbin.RData")


AdBiomodRangeSize <- BIOMOD_RangeSize(CurrentPred = proj_acadea_2000_acadea_ensemble_ROCbin ,FutureProj = proj_acadea_2020_acadea_ensemble_ROCbin)

# value 0 - stable absence between 2000-2020
# value -1 - stable presence between 2000-2020
# value 1 - potencial colonisation between 2000-2020
# value -2 - potencial extinction between 2000-2020

#Plot the Species Range Change
AdBiomodRangeSize$Diff.By.Pixel
multiple.plot (AdBiomodRangeSize$Diff.By.Pixel, acadea_2000[,2:3])

AdBiomodRangeSize$Compt.By.Models


### Ruscus aculeatus ###
load("rusacu/proj_rusacu_2020/proj_rusacu_2020_rusacu_ensemble_ROCbin.RData")
load("rusacu/proj_rusacu_2000/proj_rusacu_2000_rusacu_ensemble_ROCbin.RData")


RaBiomodRangeSize <- BIOMOD_RangeSize(CurrentPred = proj_rusacu_2000_rusacu_ensemble_ROCbin ,FutureProj = proj_rusacu_2020_rusacu_ensemble_ROCbin)

# value 0 - stable absence between 2000-2020
# value -1 - stable presence between 2000-2020
# value 1 - potencial colonisation between 2000-2020
# value -2 - potencial extinction between 2000-2020

#Plot the Species Range Change
RaBiomodRangeSize$Diff.By.Pixel
multiple.plot (RaBiomodRangeSize$Diff.By.Pixel, rusacu_2000[,2:3])

RaBiomodRangeSize$Compt.By.Models



#####################################
#### CONFLICTS BETWEEN SPECIES ######
#####################################


### Acacia dealbata ###

load("acadea/proj_acadea_2020/proj_acadea_2020_acadea_ensemble_ROCbin.RData")
load("acadea/proj_acadea_2000/proj_acadea_2000_acadea_ensemble_ROCbin.RData")

totalareaAd <- read.table("./DATA/spatAd.txt", sep="\t", h=T)

proj.cur.acadea.mean.auc.xy <- data.frame(totalareaAd[,1:3], AcadeaBinAUC= proj_acadea_2000_acadea_ensemble_ROCbin[,"acadea_EMmeanByROC_mergedAlgo_mergedRun_mergedData"])

proj.fut.acadea.mean.auc.xy <- data.frame(totalareaAd[,1:3], AcadeaBinAUC= proj_acadea_2020_acadea_ensemble_ROCbin[,"acadea_EMmeanByROC_mergedAlgo_mergedRun_mergedData"])


### Ruscus aculeatus ###

load("rusacu/proj_rusacu_2020/proj_rusacu_2020_rusacu_ensemble_ROCbin.RData")
load("rusacu/proj_rusacu_2000/proj_rusacu_2000_rusacu_ensemble_ROCbin.RData")

totalareaRa <- read.table("./DATA/spatRa.txt", sep="\t", h=T)

proj.cur.rusacu.mean.auc.xy <- data.frame(totalareaRa[,1:3], RusacuBinAUC= proj_rusacu_2000_rusacu_ensemble_ROCbin[,"rusacu_EMmeanByROC_mergedAlgo_mergedRun_mergedData"])

proj.fut.rusacu.mean.auc.xy <- data.frame(totalareaRa[,1:3], RusacuBinAUC= proj_rusacu_2020_rusacu_ensemble_ROCbin[,"rusacu_EMmeanByROC_mergedAlgo_mergedRun_mergedData"])

### Preparing data to calculate conflicts ###

df.merge.conf.cur <- merge(proj.cur.acadea.mean.auc.xy,proj.cur.rusacu.mean.auc.xy,by.x=1:3,by.y=1:3,all=F)

df.merge.conf.fut <- merge(proj.fut.acadea.mean.auc.xy,proj.fut.rusacu.mean.auc.xy,by.x=1:3,by.y=1:3,all=F)


### Current conflict ###

conflict.cur <- data.frame(df.merge.conf.cur[,1:3], conflict.cur = as.vector(df.merge.conf.cur[,"RusacuBinAUC"],mode="numeric") + (-2 * as.vector(df.merge.conf.cur[,"AcadeaBinAUC"],mode="numeric")))


level.plot (conflict.cur[,(4)], totalareaRa[,2:3], title = "Current conflict")


# value 0 - both species absent
# value -1 - potential conflict
# value 1 - only Ruscus aculatus present
# value -2 - only Acacia dealbata present


Only.Acacia.cur <- sum(conflict.cur==-2)
Only.Acacia.cur

Only.Ruscus.cur <- sum(conflict.cur==1)
Only.Ruscus.cur

Acacia.Ruscus.cur <- sum(conflict.cur==-1)
Acacia.Ruscus.cur

No.species.cur <- sum(conflict.cur==-0)
No.species.cur


### Future conflict ###
conflict.fut <- data.frame(df.merge.conf.fut[,1:3], conflict.fut = as.vector(df.merge.conf.fut[,"RusacuBinAUC"],mode="numeric") + (-2 * as.vector(df.merge.conf.fut[,"AcadeaBinAUC"],mode="numeric")))


level.plot (conflict.fut[,(4)], totalareaRa[,2:3], title = "Future conflict")

Only.Acacia.fut <- sum(conflict.fut==-2)
Only.Acacia.fut

Only.Ruscus.fut <- sum(conflict.fut==1)
Only.Ruscus.fut

Acacia.Ruscus.fut <- sum(conflict.fut==-1)
Acacia.Ruscus.fut

No.species.fut <- sum(conflict.fut==-0)
No.species.fut





