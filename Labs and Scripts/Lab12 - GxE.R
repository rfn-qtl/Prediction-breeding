##########################
# AGRO 7075	- Prediction-based Breeding
# Lab 12 - GxE
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: April 24 2023
##########################

#################################### MET analysis #########################
library(sommer)
pheno <- readRDS("pheno")
head(pheno)
blocks <- length(unique(pheno$rep))
loc <- length(unique(pheno$N))

# Fitting genotype by environment models - with a common variance (diagonal model)
fitMET <- mmer(fixed = SDM ~ 1 + N,
                       random = ~ gid + rep + gid:N,
                       rcov = ~ vsr(units),
                       data = pheno)

summary(fitMET)
# Broad-sense heritability
vpredict(fitMET, h2 ~ V1 / ( V1 + V3/loc + V4/(loc*blocks) ) ) # trials level
vpredict(fitMET, h2 ~ V1 / ( V1 + V3 + V4 ) ) # plot level


# Fitting genotype by environment models - unstructured model (US)
fitMET.US <- mmer(fixed = SDM ~ 1 + N,
               random = ~gid + rep + vsr(usr(N), gid),
               rcov = ~ vsr(units),
               data = pheno)
summary(fitMET.US)

# Broad-sense heritability
vpredict(fitMET.US, h2 ~ V1 / ( V1 + (V4)/loc + V6/(loc*blocks) ) ) # trials level
vpredict(fitMET.US, h2 ~ V1 / ( V1 + (V4) + V6 ) ) # plot level


# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
fitMET.US.H <- mmer(fixed = SDM ~ 1 + N,
                  random = ~gid + rep + vsr(usr(N), gid),
                  rcov = ~ vsr(dsr(N), units),
                  data = pheno)

summary(fitMET.US.H)
# Broad-sense heritability
vpredict(fitMET.US.H, h2 ~ V1 / ( V1 + (V4)/loc + (V6+V7)/(loc*blocks) ) ) # trials level
vpredict(fitMET.US.H, h2 ~ V1 / ( V1 + (V4) + (V6+V7) ) ) # plot level

# comparing the models
anova(fitMET, fitMET.US)
anova(fitMET.US, fitMET.US.H)

# predicting BLUPs per environment
BLUPS <- data.frame(
  MET = fitMET$U$gid$SDM,
  US = fitMET.US$U$gid$SDM,
  USH = fitMET.US.H$U$gid$SDM)

head(BLUPS)
cor(BLUPS)

# If you want to include ERM and GRM into your model
#GxE <- kronecker(ERM, Ga)
#pheno$GN <- paste(pheno$gid, pheno$N)

# Fitting genotype by environment models - with a kernels (GRM and ERM)
#fitMET <- mmer(fixed = SDM ~ 1 + rep,
#               random = ~ vsr(N, ERM) + vsr(gid, Ga) + vsr(GN, GxE),
#               rcov = ~ vsr(units),
#               data = pheno)


####################### Stability and adaptability - Finlay-Wilkinson #########################
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)

library(statgenGxE)

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = pheno, genotype = "gid", trial = "N")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "SDM")
summary(dropsFW)

# let's take a look at the output
names(dropsFW)
dropsFW$estimates
dropsFW$envEffs
dropsFW$fittedGeno

## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")

############################# GGE-Biplot Analysis ######################################
library(metan)

model.gge <- gge(pheno, N, gid, SDM, svp = "symmetrical")

model.gge$SDM

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where

######## the end ##########