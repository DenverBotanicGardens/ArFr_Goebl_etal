## A. Goebl - 
## Script updated 2025-11-02 
## Analysis of Artemisia frigida trait data from common garden experiment  


rm(list=ls())


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(tidyr)
library(corrplot)
library(dplyr)
library(stringr)
library(lme4)
library(car)
library(plotrix)
library(emmeans)

calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------





## SET WORKING DIRECTORY --------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
ARFR22 <- read.csv(file="20251031_ChatfieldDataClean2022_ARFR.csv", sep=",", header=TRUE, dec=".")
ARFR23 <- read.csv(file="20251031_ChatfieldDataClean2023_ARFR.csv", sep=",", header=TRUE, dec=".")
ARFR24 <- read.csv(file="20251031_ChatfieldDataClean2024_ARFR.csv", sep=",", header=TRUE, dec=".")

## Load PCA values from AE's analysis
pca_vals <- read.csv(file="20251006_pcaTableFromAE_ARFR.csv", sep=",", header=TRUE, dec=".")
## ----------------------------------------------------------------------------------------------




## ARFR - DATA RE-FORMAT AS NEEDED --------------------------------------------------------------
str(ARFR22)
str(ARFR23)
str(ARFR24)

ARFR22$Source <- as.factor(ARFR22$Source)
ARFR23$Source <- as.factor(ARFR23$Source)
ARFR24$Source <- as.factor(ARFR24$Source)

ARFR22$Treatment <- as.factor(ARFR22$Treatment)
ARFR23$Treatment <- as.factor(ARFR23$Treatment)
ARFR24$Treatment <- as.factor(ARFR24$Treatment)
## ----------------------------------------------------------------------------------------------




## ----------------------------------------------------------------------------------------------
## ARFR - PREPARE RESPONSE VARIABLES

## CHANGE DATA TO NA BASED ON VARIOUS CONDITIONS (E.G. EXCLUDE SURV DATA IN 2023-24 IF HARVESTED IN 2022)
ARFR22.ex <- ARFR22 %>% mutate(across(c(starts_with("Survival_"),starts_with("Length_")), 
             ~case_when(ExcludeBcNotReplaced=="Y" ~ as.numeric(NA), TRUE ~ as.numeric(.x))))
ARFR23.ex <- ARFR23 %>% mutate(across(c(starts_with("Survival_"), "Height_20230927"), 
             ~case_when(ARFR22.ExcludeBcNotReplaced=="Y" ~ as.numeric(NA), TRUE ~ as.numeric(.x))))
ARFR23.ex <- ARFR23.ex %>% mutate(across(c(starts_with("Survival_"), "Height_20230927"), 
             ~case_when(ARFR24.ExcludeBcHarvest=="Y" ~ as.numeric(NA), TRUE ~ as.numeric(.x))))
ARFR24.ex <- ARFR24 %>% mutate(across(c("Survival","LeafSurfaceArea_cm2","SLA_mm2permg","DryLeafMass_g",
             "InfBM2022smpls_HEADS_2024weigh","InfBM2022smpls_CHAFF_2024weigh"), 
            ~case_when(ARFR22.ExcludeBcNotReplaced=="Y" ~ as.numeric(NA), TRUE ~ as.numeric(.x))))
ARFR24.ex <- ARFR24.ex %>% mutate(across(c("Survival","LeafSurfaceArea_cm2","SLA_mm2permg","DryLeafMass_g",
             "InfBM2022smpls_HEADS_2024weigh","InfBM2022smpls_CHAFF_2024weigh"), 
            ~case_when(ARFR23.ExcludeSurvDueToInconsistData=="Y" ~ as.numeric(NA), TRUE ~ as.numeric(.x))))
ARFR24.ex <- ARFR24.ex %>% mutate(across(c("Survival","LeafSurfaceArea_cm2","SLA_mm2permg","DryLeafMass_g"), 
            ~case_when(ExcludeBcHarvest=="Y" ~ as.numeric(NA), TRUE ~ as.numeric(.x))))






## COMBINE AND ADD REPRODUCTION BIOMASS DATA --------------------------------------------------------------
## Combine flower head and chaff/seed weights + any missed samples from initial 2022 weights
## Change Chaff entries to zero (from NA) if no chaff weight, but heads were weighed
ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24.ex$InfBM2022smpls_HEADS_2024weigh) & is.na(ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh)] <- 0
## Add chaff and head weights together
ARFR24.ex$InfBM2022_2024updated <- ARFR24.ex$InfBM2022smpls_HEADS_2024weigh + ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh
## Add several individuals (524, 885, and 908) from 2022 weights that weren't available for 2024 re-weigh
ARFR24.ex$InfBM2022_2024updated[ARFR24.ex$ID==524] <- ARFR24.ex$InfBM2022_Wobag_g[ARFR24.ex$ID==524]
ARFR24.ex$InfBM2022_2024updated[ARFR24.ex$ID==885] <- ARFR24.ex$InfBM2022_Wobag_g[ARFR24.ex$ID==885]
ARFR24.ex$InfBM2022_2024updated[ARFR24.ex$ID==908] <- ARFR24.ex$InfBM2022_Wobag_g[ARFR24.ex$ID==908]



## CLEAN 2023 PLANT SIZE FIELD MEASUREMENTS -------------------------------------------------------
identical(ARFR23.ex$ExcludeSurvDueToInconsistData, ARFR23.ex$ExcludeSzDueToUncertainty)
nrow(ARFR23.ex[!is.na(ARFR23.ex$ExcludeSurvDueToInconsistData),])
nrow(ARFR23.ex[!is.na(ARFR23.ex$ExcludeSzDueToUncertainty),])

nrow(ARFR23.ex[!is.na(ARFR23.ex$Height_20230927>0),])
ARFR23.ex$Height_20230927[ARFR23.ex$ExcludeSzDueToUncertainty=="Y" & !is.na(ARFR23.ex$ExcludeSzDueToUncertainty)] <- NA
ARFR23.ex$Height_20230927[ARFR23.ex$ExcludeSurvDueToInconsistData=="Y" & !is.na(ARFR23.ex$ExcludeSurvDueToInconsistData)] <- NA
nrow(ARFR23[!is.na(ARFR23$Height_20230927>0),])




## COMBINE RELEVANT 2022, 2023, 2024 DATA
ARFR22.sel <- ARFR22.ex %>% dplyr::select(c("ID", "Length_cm_20220726","Survival_20220922"))               
ARFR23.sel <- ARFR23.ex %>% dplyr::select(c("ID","Height_20230927","Survival_20230801")) 
ARFR.sel <- left_join(ARFR24.ex, ARFR23.sel, by="ID") 
ARFR.sel <- left_join(ARFR.sel, ARFR22.sel, by="ID") 
ARFR.sel$Source <- as.factor(ARFR.sel$Source)
# --------------------------------------------------------------------------------------------------





## ARFR - VISUALIZE RAW DATA ---------------------------------------------------------------
AddnCols <- as.data.frame(cbind(as.character(ARFR.sel$PopAbbrev),as.character(ARFR.sel$PopCol),as.numeric(ARFR.sel$Lat),as.character(ARFR.sel$Source)))
colnames(AddnCols) <- c("PopAbbrev","PopCol","Lat","Source")


## Order populations for plotting 
## Order by average latitude (or other traits)
ARFR.latByMed <- with(ARFR.sel, reorder(Source, Lat, median, na.rm=TRUE))

ARFR.meds <- ARFR.sel %>% group_by(Source) %>% 
             dplyr::summarise(Height22_MD=median(Length_cm_20220726,na.rm=TRUE), AGB22_MD=median(AGB2022_MinusBag,na.rm=TRUE),
             ReproBMrw_MD=median(InfBM2022_2024updated,na.rm=TRUE), Height23_MD=median(Height_20230927,na.rm=TRUE),
             SLA_MD=median(SLA_mm2permg,na.rm=TRUE), Surv24_Count=length(na.omit(Survival)), Surv24_Sum=sum(Survival, na.rm=TRUE),
             Surv23_Count=length(na.omit(Survival_20230801)), Surv23_Sum=sum(Survival_20230801, na.rm=TRUE),
             Surv22_Count=length(na.omit(Survival_20220922)), Surv22_Sum=sum(Survival_20220922, na.rm=TRUE))

ARFR.meds <- left_join(ARFR.meds, AddnCols, by="Source")
ARFR.meds <- unique(ARFR.meds)

## Estimate survival each year
surv24.pop <- ARFR.meds$Surv24_Sum/ARFR.meds$Surv24_Count
surv23.pop <- ARFR.meds$Surv23_Sum/ARFR.meds$Surv23_Count
surv22.pop <- ARFR.meds$Surv22_Sum/ARFR.meds$Surv22_Count



## Boxplots of raw data 
ARFR.meds <- ARFR.meds[order(ARFR.meds$Lat),] #Order by latitude

par(mfrow=c(2,3))

## Size 2022
boxplot(Length_cm_20220726 ~ ARFR.latByMed, data=ARFR.sel,
        ylab="Height (cm)", xlab=NA, cex.lab=1.25, horizontal=FALSE,
        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2,
        main="FINAL SIZE 2022", cex.main=1.5, col=ARFR.meds$PopCol)


## Reproduction
boxplot(InfBM2022_2024updated ~ ARFR.latByMed, data=ARFR.sel, las=2,
        ylab="Reproductive biomass (g)", xlab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=FALSE, ylim=c(0,80),
        cex.main=1.5, col=ARFR.meds$PopCol, main="REPRODUCTIVE OUTPUT 2022")


## Size 2023
boxplot(Height_20230927 ~ ARFR.latByMed, data=ARFR.sel,
        ylab="Height (cm)", xlab=NA, cex.lab=1.25, horizontal=FALSE,
        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2, ylim=c(15,90),
        main="FINAL SIZE 2023", cex.main=1.5, col=ARFR.meds$PopCol)


## SLA 2024
boxplot(SLA_mm2permg ~ ARFR.latByMed, data=ARFR.sel, las=2,
        ylab="Specific leaf area (mm2/mg)", xlab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=FALSE, ylim=c(0,40),
        cex.main=1.5, col=ARFR.meds$PopCol, main="SPECIFIC LEAF AREA 2024")


## Survival 2024
barplot(surv24.pop, col=ARFR.meds$PopCol, ylim=c(0,1), cex.axis=0.99, names.arg=ARFR.meds$PopAbbrev,
        las=2, ylab="Survival rate", main="SURVIVAL 2022-2024", cex.main=1.5)
## ---------------------------------------------------






## MODEL TRAIT DATA ----------------------------------

## Re-order Source as factor before running models
AddnCols.unq <- unique(AddnCols)
AddnCols.unq <- AddnCols.unq[order(AddnCols.unq$Lat),] #Order by lat

ARFR.sel$Source <- factor(ARFR.sel$Source, levels=AddnCols.unq$Source)



## Plant size

## 2022
ARFR.sel$Block <- as.factor(ARFR.sel$Block)
sz22.mod <- lmer(Length_cm_20220726 ~ Source + (1|Block), data=ARFR.sel)
summary(sz22.mod)
Anova(sz22.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(sz22.mod, type="pearson")
dResid <- residuals(sz22.mod, type="deviance")
hist(pResid)                                          #Shape should be consistent with assumed error distribution (e.g. normal)
hist(dResid)
qqnorm(pResid)                                        #Points should roughly follow the diagonal line, even at tails
qqline(dResid)
plot(fitted(sz22.mod), pResid, abline(h=0,col="red")) #Residuals should be randomly scattered around 0 line

## Obtain model predicted values for response variables
predForSource <- as.data.frame(AddnCols.unq$Source) 
colnames(predForSource) <- "Source"
sz22.pred <- predict(sz22.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)


## 2023
sz23.mod <- lmer(Height_20230927 ~ Source + (1|Block), data=ARFR.sel)
summary(sz23.mod)
Anova(sz23.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(sz23.mod, type="pearson")
hist(pResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(sz23.mod), pResid, abline(h=0,col="red")) 

## Obtain model predicted values for response variables
sz23.pred <- predict(sz23.mod, newdata=predForSource, se.fit=TRUE, type="response", re.form=~0)



## SLA
hist(ARFR.sel$SLA_mm2permg)
hist(log(ARFR.sel$SLA_mm2permg))
sla.mod <- lmer(log(SLA_mm2permg) ~ Source + (1|Block), data=ARFR.sel)
summary(sla.mod)
Anova(sla.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(sla.mod, type="pearson")
dResid <- residuals(sla.mod, type="deviance")
hist(dResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(sla.mod), pResid, abline(h=0,col="red")) 

## Obtain model predicted values for response variables
sla.predLog <- predict(sla.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)
sla.predOrigFit <- exp(sla.predLog$fit)
sla.predOrigSE <- exp(sla.predLog$se.fit)



## Reproductive biomass
hist(ARFR.sel$InfBM2022_2024updated)
hist(log(ARFR.sel$InfBM2022_2024updated))

rbm.mod <- lmer(log(InfBM2022_2024updated) ~ Source + (1|Block), data=ARFR.sel)
summary(rbm.mod)
Anova(rbm.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(rbm.mod, type="pearson")
hist(pResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(rbm.mod), pResid, abline(h=0,col="red")) 

## Obtain model predicted values for response variables
rbm.predLog <- predict(rbm.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)
rbm.predOrigFit <- exp(rbm.predLog$fit)
rbm.predOrigSE <- exp(rbm.predLog$se.fit)



## Survival
surv24.mod <- glmer(Survival ~ Source + (1|Block), data=ARFR.sel, family=binomial(link="logit"))
summary(surv24.mod)
Anova(surv24.mod)

## Obtain model predicted values for response variables
surv24.pred <- predict(surv24.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)






## PLOT MODEL ESTIMATES AND SE
predForSource <- dplyr::left_join(predForSource, AddnCols, by="Source")
predForSource <- unique(predForSource)
preds <- cbind(predForSource, sz22.pred$fit, sz22.pred$se.fit, sz23.pred$fit, sz23.pred$se.fit,
               rbm.predOrigFit, rbm.predOrigSE, sla.predOrigFit, sla.predOrigSE, surv24.pred$fit, surv24.pred$se.fit)

par(mfrow=c(1,1))
plot(NA, NA, xlab="Seed source", ylab="Height (cm)",
     main="FINAL SIZE 2022", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(13.5,48.5))
arrows(1:11, preds$`sz22.pred$fit`+preds$`sz22.pred$se.fit`, 1:11, preds$`sz22.pred$fit`-preds$`sz22.pred$se.fit`,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$`sz22.pred$fit`, col="black", bg=preds$PopCol, pch=21, cex=1.45)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)


par(mfrow=c(2,2))
plot(NA, NA, xlab="Seed source", ylab="Reproductive  biomass (g)",
     main="REPRODUCTION 2022", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(0,40))
arrows(1:11, preds$rbm.predOrigFit+preds$rbm.predOrigSE, 1:11, preds$rbm.predOrigFit-preds$rbm.predOrigSE,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$rbm.predOrigFit, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)

plot(NA, NA, xlab="Seed source", ylab="Height (cm)",
     main="FINAL SIZE 2023", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(40,68))
arrows(1:11, preds$`sz23.pred$fit`+preds$`sz23.pred$se.fit`, 1:11, preds$`sz23.pred$fit`-preds$`sz23.pred$se.fit`,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$`sz23.pred$fit`, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)

plot(NA, NA, xlab="Seed source", ylab="Specific leaf area (mm2/mg)",
     main="SPECIFIC LEAF AREA 2024", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(10,16))
arrows(1:11, preds$sla.predOrigFit+preds$sla.predOrigSE, 1:11, preds$sla.predOrigFit-preds$sla.predOrigSE,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$sla.predOrigFit, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)

plot(NA, NA, xlab="Seed source", ylab="Survival rate",
     main="SURVIVAL 2022-2024", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(0,1.9))
arrows(1:11, preds$`surv24.pred$fit`+preds$`surv24.pred$se.fit`, 1:11, preds$`surv24.pred$fit`-preds$`surv24.pred$se.fit`,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$`surv24.pred$fit`, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)



## Look at correlation in trait estimates
plot(preds$`sz22.pred$fit`, preds$`sz23.pred$fit`)
plot(preds$`sz22.pred$fit`, preds$rbm.predOrigFit)
plot(preds$`sz23.pred$fit`, preds$rbm.predOrigFit)



## Look at pairwise differences in model estimates + sig between populations using emmeans 
sz23.pw <- emmeans(sz23.mod, specs = pairwise ~ Source, type="response")
rbm.pw <- emmeans(rbm.mod, specs = pairwise ~ Source, type="response")
sla.pw <- emmeans(sla.mod, specs = pairwise ~ Source, type="response")
surv.pw <- emmeans(surv24.mod, specs = pairwise ~ Source, type="response")
## ---------------------------------------------------------------------








## Trait PCA -----------------------------------------------------------
ARFR.traits <- ARFR.sel %>% dplyr::select(c("Survival","Length_cm_20220726","Height_20230927", "SLA_mm2permg",
                                           "InfBM2022_2024updated")) 
                                            
ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Length_cm_20220726),] #Remove indivs that died early & have no data
ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Survival),] #Remove rows without surv data

## Look at trait correlations
ARFR.traitsCor <- cor(ARFR.traits, use="pairwise.complete.obs")
corrplot(ARFR.traitsCor)

ARFR.traitsT <- t(ARFR.traits)

## Make covariance matrix and run pca
covMat.traits <- cov(ARFR.traitsT, use="pairwise.complete.obs")
pca.results <- prcomp(covMat.traits, center=TRUE)


## Get sample list with pop ID and colors
ARFR.indivPop <- ARFR.sel %>% dplyr::select(c("Source", "ID", "HexCode_Indv"))
ARFR.indivPop$ID <- as.factor(ARFR.indivPop$ID)
indivs.traitPCA <- as.factor(colnames(ARFR.traitsT))
indivs.traitPCA <- as.data.frame(indivs.traitPCA)
colnames(indivs.traitPCA) <- "ID"
indivs.traitPCA <- left_join(indivs.traitPCA, ARFR.indivPop, by="ID")



par(mfrow=c(1,1))
plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv, main="Trait PCA")
plot(x=pca.results$x[,2], y=pca.results$x[,3],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv)



## Calculate PC1 mean values for each source population
trait.PCscores <- as.data.frame(cbind(pca.results$x[,1], pca.results$x[,2], as.character(indivs.traitPCA$Source)))
colnames(trait.PCscores) <- c("PC1", "PC2", "Source")
trait.PCscores$PC1 <- as.numeric(trait.PCscores$PC1)
traitPC1.mean <- trait.PCscores %>% group_by(Source) %>% summarise(PC1mean = mean(PC1), n=n())

## Create color gradient and assign colors based on numeric continuous PC1 mean values
# Adapted from ChatGPT generated code
# Define a color gradient (e.g., from blue to red)
gradient_fn <- colorRamp(c("greenyellow",   "deeppink"))

# Normalize values to [0,1] scale
vals_norm <- (traitPC1.mean$PC1mean - min(traitPC1.mean$PC1mean)) / (max(traitPC1.mean$PC1mean) - min(traitPC1.mean$PC1mean))

# Get RGB colors (as integers 0-255)
rgb_matrix <- gradient_fn(vals_norm)

# Convert to hex color strings
colors.traitPC <- rgb(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3], maxColorValue = 255)

# Plot using colors
plot(traitPC1.mean$PC1mean, rep(1, length(traitPC1.mean$PC1mean)), col=colors.traitPC, pch=16, cex=1.75,
     xlab="Trait PC1 score", ylab=NA, main="Color representation of trait PC1 scores", yaxt='n')

traitPC1.mean$color <- colors.traitPC
## -----------------------------------------------------------------------------







### VCF table and genomic PCA  --------------------------------------------------------------------------
## Get list of sample names 
indvNames <- as.data.frame(as.character(pca_vals$sample.id))
colnames(indvNames) <- "Sample"

## Make column with ID using string replace 
indvNames$Temp <- str_replace(indvNames$Sample, "ARFR_", "")
indvNames$ID <- as.integer(str_replace(indvNames$Temp, "_sorted", ""))
## Join by ID to get source (pop ID)
indvNames <- left_join(indvNames, ARFR24, by="ID")




## PCA scores  
## All samples
pca_vals$Sample <- indvNames$Sample               
pca_vals$Source <- indvNames$Source                          #Add a column with pop ID
popNames <- unique(indvNames$PopID)


## Calculate mean PC1 values for each population
PC1.mean <- pca_vals %>% group_by(Source) %>% summarise(PC1mean = mean(EV1), n=n())
PC1.mean <- PC1.mean[1:11,]




## Create color gradient and assign colors based on numeric continuous PC1 mean values
# Adapted from ChatGPT generated code
# Define a color gradient 
gradient_fn <- colorRamp(c("greenyellow",   "deeppink"))

# Normalize your values to [0,1] scale
vals_norm <- (PC1.mean$PC1mean - min(PC1.mean$PC1mean)) / (max(PC1.mean$PC1mean) - min(PC1.mean$PC1mean))

# Get RGB colors (as integers 0-255)
rgb_matrix <- gradient_fn(vals_norm)

# Convert to hex color strings
colors <- rgb(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3], maxColorValue = 255)

# Plot using colors
plot(PC1.mean$PC1mean, rep(1, length(PC1.mean$PC1mean)), col=colors, pch=16, cex=1.75,
     xlab="Genomic PC1 score", ylab=NA, main="Color representation of genomic PC1 scores", yaxt='n')

PC1.mean$color <- colors

## Plot range of color gradient as a legend
genPCrange <- seq(from=min(PC1.mean$PC1mean), to=max(PC1.mean$PC1mean), by=0.01)
vals_normGenRange <- (genPCrange - min(genPCrange)) / (max(genPCrange) - min(genPCrange))
rgb_matrixGenRange <- gradient_fn(vals_normGenRange)
colors.genPCrange <- rgb(rgb_matrixGenRange[,1], rgb_matrixGenRange[,2], rgb_matrixGenRange[,3], maxColorValue = 255)
plot(genPCrange, rep(0.5, length(genPCrange)), col=colors.genPCrange, pch=15, cex=4)
## -------------------------------------------------------------------------------------------









