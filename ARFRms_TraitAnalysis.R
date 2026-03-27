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
library(FactoMineR)

calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------





## SET WORKING DIRECTORY --------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
ARFR22 <- read.csv(file="20251031_ChatfieldDataClean2022_ARFRms.csv", sep=",", header=TRUE, dec=".")
ARFR23 <- read.csv(file="20251031_ChatfieldDataClean2023_ARFRms.csv", sep=",", header=TRUE, dec=".")
ARFR24 <- read.csv(file="20251031_ChatfieldDataClean2024_ARFRms.csv", sep=",", header=TRUE, dec=".")

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






## COMBINE AND ADD REPRODUCTIVE BIOMASS DATA --------------------------------------------------------------
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

ARFR.sel$PopAbbrev <- as.character(str_replace_all(ARFR.sel$PopAbbrev, "\\.", ""))

AddnCols <- as.data.frame(cbind(as.character(ARFR.sel$PopAbbrev),as.character(ARFR.sel$PopCol),as.numeric(ARFR.sel$Lat),as.character(ARFR.sel$Source)))
colnames(AddnCols) <- c("PopAbbrev","PopCol","Lat","Source")
# --------------------------------------------------------------------------------------------------





## ARFR - VISUALIZE RAW DATA ---------------------------------------------------------------

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


## Survival 2024 Barplot
barplot(surv24.pop, col=ARFR.meds$PopCol, ylim=c(0,1), cex.axis=0.99, names.arg=ARFR.meds$PopAbbrev,
        las=2, ylab="Survival rate", main="SURVIVAL 2022-2024", cex.main=1.5)
## ---------------------------------------------------






## MODEL TRAIT DATA ----------------------------------
library(DHARMa)

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
hist(pResid)                                          #Shape should be consistent with assumed error distribution (e.g. normal)
qqnorm(pResid)                                        #Points should roughly follow the diagonal line, even at tails
qqline(pResid)
plot(fitted(sz22.mod), pResid, abline(h=0,col="red")) #Residuals should be randomly scattered around 0 line
sz22.simRes <- simulateResiduals(sz22.mod)
plot(sz22.simRes)

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
sz23.simRes <- simulateResiduals(sz23.mod)
plot(sz23.simRes)

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
hist(pResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(sla.mod), pResid, abline(h=0,col="red")) 
sla.simRes <- simulateResiduals(sla.mod)
plot(sla.simRes)

## Obtain model predicted values for response variables
sla.predLog <- predict(sla.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)
sla.predOrigFit <- exp(sla.predLog$fit)
sla.predOrigSE <- exp(sla.predLog$se.fit)



## Reproductive biomass
hist(ARFR.sel$InfBM2022_2024updated)
hist(log(ARFR.sel$InfBM2022_2024updated))

rbm.mod <- lmer(InfBM2022_2024updated ~ Source + (1|Block), data=ARFR.sel)

summary(rbm.mod)
Anova(rbm.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(rbm.mod, type="pearson")
hist(pResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(rbm.mod), pResid, abline(h=0,col="red")) 
rbm.simRes <- simulateResiduals(rbm.mod)
plot(rbm.simRes)

## Obtain model predicted values for response variables
rbm.pred <- predict(rbm.mod, newdata=predForSource, se.fit=TRUE, type="response", re.form=~0)


## Try model without zeros for reproduction
ARFR.sel$InfBM2022_2024updated[ARFR.sel$InfBM2022_2024updated == 0 & !is.na(ARFR.sel$InfBM2022_2024updated)] <- NA

rbm.modLog <- lmer(log(InfBM2022_2024updated) ~ Source + (1|Block), data=ARFR.sel)
rbm.modGam <- glmer(InfBM2022_2024updated ~ Source + (1|Block), family=Gamma(link="log"), data=InfBMno0)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(rbm.modLog, type="pearson")
hist(pResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(rbm.modLog), pResid, abline(h=0,col="red")) 
rbmLog.simRes <- simulateResiduals(rbm.modLog)
plot(rbmLog.simRes)

pResid <- residuals(rbm.modGam, type="pearson")
hist(pResid)                                          
qqnorm(pResid)                                        
qqline(pResid)
plot(fitted(rbm.modGam), pResid, abline(h=0,col="red")) 
rbmGam.simRes <- simulateResiduals(rbm.modGam)
plot(rbmGam.simRes)

## Obtain model predicted values for response variables
rbm.predLog <- predict(rbm.modLog, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)
rbm.predOrigFit <- exp(rbm.predLog$fit)
rbm.predOrigSE <- exp(rbm.predLog$se.fit)

rbm.predGam <- predict(rbm.modGam, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)





## Survival
surv24.mod <- glmer(Survival ~ Source + (1|Block), data=ARFR.sel, family=binomial(link="logit"))
summary(surv24.mod)
Anova(surv24.mod)

## Check for overdispersion (modified from ChatGPT) and other diagnostics
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type="pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  p <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  return(c(chisq = Pearson.chisq, ratio=ratio, rdf=rdf, p=p))
}

overdisp_fun(surv24.mod)

surv.simRes <- simulateResiduals(surv24.mod)
plot(surv.simRes)
testDispersion(surv.simRes)
testZeroInflation(surv.simRes)

VarCorr(surv24.mod) #Check random effects


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
ARFR.sel$Grwth <- ARFR.sel$Height_20230927 - ARFR.sel$Length_cm_20220726
ARFR.traits <- ARFR.sel %>% dplyr::select(c("Length_cm_20220726","Height_20230927","Grwth","Survival","SLA_mm2permg",
                                            "InfBM2022_2024updated","Source")) 

## Clean up and subset dataset due to missing values 
ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Length_cm_20220726),] #Remove indivs that died early & have no data
ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Survival) | !is.na(ARFR.traits$InfBM2022_2024updated),] #Remove rows without surv data except for plts w InfBM
ARFR.traits <- ARFR.traits[(!is.na(ARFR.traits$Height_20230927) | !is.na(ARFR.traits$SLA_mm2permg)) 
                           | !is.na(ARFR.traits$InfBM2022_2024updated),] #Keep rows that have InfBM and either 2023 sz or SLA

## Look at trait correlations
ARFR.traitsCor <- cor(ARFR.traits[,1:6], use="pairwise.complete.obs")
corrplot(ARFR.traitsCor)
## Exclude size in 2022 due to high correlation with growth and 2022 repro

## Calculate population means and use to impute missing values (modified from ChatGPT)
ARFR.traitsImp <- ARFR.traits[,2:7] %>% group_by(Source) %>%
  mutate(Grwth=ifelse(is.na(Grwth), mean(Grwth, na.rm=TRUE), Grwth),
  Sz23=ifelse(is.na(Height_20230927), mean(Height_20230927, na.rm=TRUE), Height_20230927),
  SLA=ifelse(is.na(SLA_mm2permg), mean(SLA_mm2permg, na.rm=TRUE), SLA_mm2permg),
  RBM=ifelse(is.na(InfBM2022_2024updated), mean(InfBM2022_2024updated, na.rm=TRUE), InfBM2022_2024updated)
  ) %>%
  ungroup()

ARFR.traitsImp <- ARFR.traitsImp %>% dplyr::select(c("Sz23","Grwth","Survival","SLA","RBM")) 

## Run PCA
pca.results <- PCA(ARFR.traitsImp, scale.unit=TRUE, graph=TRUE)


## Get sample list with pop ID and colors
ARFR.indivPop <- ARFR.sel %>% dplyr::select(c("Source", "ID", "HexCode_Indv"))
ARFR.indivPop$ID <- as.factor(ARFR.indivPop$ID)
indivs.traitPCA <- as.factor(rownames(ARFR.traits))
indivs.traitPCA <- as.data.frame(indivs.traitPCA)
colnames(indivs.traitPCA) <- "ID"
indivs.traitPCA <- left_join(indivs.traitPCA, ARFR.indivPop, by="ID")

pca.ind <- pca.results$ind
pca.load <- pca.results$var$coord
pca.eig <- pca.results$eig



## Plot
par(mfrow=c(1,1))
par(pty="s")

plot(x=pca.ind$coord[,1], y=pca.ind$coord[,2],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv, main="Trait PCA",
     xlab="PC1 (28.8% variance)", ylab="PC2 (26.1% variance)")
legend("topleft", AddnCols.unq$PopAbbrev, col=AddnCols.unq$PopCol, cex=0.95, pch=19)

## Add loadings to plot (modified from chatGPT)
rownames(pca.load)
traits.loadNames <- c("Size","Growth","Survival","SLA","Reproduction")
arrow_scale <- 4
arrows(0,0, pca.load[,1]*arrow_scale, pca.load[,2]*arrow_scale, length=0.1, col="black", lwd=2)

text(pca.load[,1]*arrow_scale, pca.load[,2]*arrow_scale,
     labels=traits.loadNames, pos=2, cex=0.9)


## Plot PC 3 and 4
plot(x=pca.ind$coord[,3], y=pca.ind$coord[,4],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv, main="Trait PCA")
## ------------



## Calculate mean PC values for each population
trait.PCscores <- as.data.frame(cbind(pca.ind$coord[,1], pca.ind$coord[,2], as.character(indivs.traitPCA$Source)))
colnames(trait.PCscores) <- c("PC1", "PC2", "Source")

trait.PCscores$PC1 <- as.numeric(trait.PCscores$PC1)
traitPC1.mean <- trait.PCscores %>% group_by(Source) %>% summarise(PC1mean = mean(PC1), n=n())

trait.PCscores$PC2 <- as.numeric(trait.PCscores$PC2)
traitPC2.mean <- trait.PCscores %>% group_by(Source) %>% summarise(PC2mean = mean(PC2), n=n())


## Trait PC1
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
dev.off()
plot(traitPC1.mean$PC1mean, rep(1, length(traitPC1.mean$PC1mean)), col=colors.traitPC, pch=16, cex=2,
     xlab="Trait PC1 score", ylab=NA, main="Color representation of trait PC1 scores", yaxt='n', cex.main=1.5)

traitPC1.mean$color <- colors.traitPC

## Plot range of color gradient as a legend
traitPCrange <- seq(from=min(traitPC1.mean$PC1mean), to=max(traitPC1.mean$PC1mean), by=0.01)
vals_normTraitRange <- (traitPCrange - min(traitPCrange)) / (max(traitPCrange) - min(traitPCrange))
rgb_matrixTraitRange <- gradient_fn(vals_normTraitRange)
colors.traitPCrange <- rgb(rgb_matrixTraitRange[,1], rgb_matrixTraitRange[,2], rgb_matrixTraitRange[,3], maxColorValue = 255)
plot(traitPCrange, rep(0.5, length(traitPCrange)), col=colors.traitPCrange, pch=15, cex=4)



## PC2
# Normalize values to [0,1] scale
vals_normTraitPC2 <- (traitPC2.mean$PC2mean - min(traitPC2.mean$PC2mean)) / (max(traitPC2.mean$PC2mean) - min(traitPC2.mean$PC2mean))

# Get RGB colors (as integers 0-255)
rgb_matrixTraitPC2 <- gradient_fn(vals_normTraitPC2)

# Convert to hex color strings
colorsTraitPC2 <- rgb(rgb_matrixTraitPC2[,1], rgb_matrixTraitPC2[,2], rgb_matrixTraitPC2[,3], maxColorValue = 255)

# Plot using colors
dev.off()
plot(traitPC2.mean$PC2mean, rep(1, length(traitPC2.mean$PC2mean)), col=colorsTraitPC2, pch=16, cex=1.75,
     xlab="Trait PC2 score", ylab=NA, main="Color representation of trait PC2 scores", yaxt='n')

traitPC2.mean$color <- colorsTraitPC2

## Plot range of color gradient as a legend
traitPC2range <- seq(from=min(traitPC2.mean$PC2mean), to=max(traitPC2.mean$PC2mean), by=0.01)
vals_normTraitPC2Range <- (traitPC2range - min(traitPC2range)) / (max(traitPC2range) - min(traitPC2range))
rgb_matrixTraitPC2Range <- gradient_fn(vals_normTraitPC2Range)
colors.traitPC2range <- rgb(rgb_matrixTraitPC2Range[,1], rgb_matrixTraitPC2Range[,2], rgb_matrixTraitPC2Range[,3], maxColorValue = 255)
plot(traitPC2range, rep(0.5, length(traitPC2range)), col=colors.traitPC2range, pch=15, cex=4)
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


## Calculate mean PC values for each population
PC1.mean <- pca_vals %>% group_by(Source) %>% summarise(PC1mean = mean(EV1), n=n())
PC1.mean <- PC1.mean[1:11,]

PC2.mean <- pca_vals %>% group_by(Source) %>% summarise(PC2mean = mean(EV2), n=n())
PC2.mean <- PC2.mean[1:11,]



## Create color gradient and assign colors based on numeric continuous PC1 mean values
# Adapted from ChatGPT generated code
# Define a color gradient 
gradient_fn <- colorRamp(c("greenyellow",   "deeppink"))

# Normalize values to [0,1] scale
vals_norm <- (PC1.mean$PC1mean - min(PC1.mean$PC1mean)) / (max(PC1.mean$PC1mean) - min(PC1.mean$PC1mean))

# Get RGB colors (as integers 0-255)
rgb_matrix <- gradient_fn(vals_norm)

# Convert to hex color strings
colors <- rgb(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3], maxColorValue = 255)

# Plot using colors
dev.off()
plot(PC1.mean$PC1mean, rep(1, length(PC1.mean$PC1mean)), col=colors, pch=16, cex=1.75,
     xlab="Genomic PC1 score", ylab=NA, main="Color representation of genomic PC1 scores", yaxt='n')

PC1.mean$color <- colors

## Plot range of color gradient as a legend
genPCrange <- seq(from=min(PC1.mean$PC1mean), to=max(PC1.mean$PC1mean), by=0.01)
vals_normGenRange <- (genPCrange - min(genPCrange)) / (max(genPCrange) - min(genPCrange))
rgb_matrixGenRange <- gradient_fn(vals_normGenRange)
colors.genPCrange <- rgb(rgb_matrixGenRange[,1], rgb_matrixGenRange[,2], rgb_matrixGenRange[,3], maxColorValue = 255)
plot(genPCrange, rep(0.5, length(genPCrange)), col=colors.genPCrange, pch=15, cex=4)



## PC2
# Normalize values to [0,1] scale
vals_normPC2 <- (PC2.mean$PC2mean - min(PC2.mean$PC2mean)) / (max(PC2.mean$PC2mean) - min(PC2.mean$PC2mean))

# Get RGB colors (as integers 0-255)
rgb_matrixPC2 <- gradient_fn(vals_normPC2)

# Convert to hex color strings
colorsPC2 <- rgb(rgb_matrixPC2[,1], rgb_matrixPC2[,2], rgb_matrixPC2[,3], maxColorValue = 255)

# Plot using colors
dev.off()
plot(PC2.mean$PC2mean, rep(1, length(PC2.mean$PC2mean)), col=colorsPC2, pch=16, cex=1.75,
     xlab="Genomic PC2 score", ylab=NA, main="Color representation of genomic PC2 scores", yaxt='n')

PC2.mean$color <- colorsPC2

## Plot range of color gradient as a legend
genPCrange <- seq(from=min(PC2.mean$PC2mean), to=max(PC2.mean$PC2mean), by=0.01)
vals_normGenRange <- (genPCrange - min(genPCrange)) / (max(genPCrange) - min(genPCrange))
rgb_matrixGenRange <- gradient_fn(vals_normGenRange)
colors.genPCrange <- rgb(rgb_matrixGenRange[,1], rgb_matrixGenRange[,2], rgb_matrixGenRange[,3], maxColorValue = 255)
plot(genPCrange, rep(0.5, length(genPCrange)), col=colors.genPCrange, pch=15, cex=4)
## -------------------------------------------------------------------------------------------









