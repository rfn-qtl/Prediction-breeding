########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 10 - GS using mixed models
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Dec 8 2022
#######################################

# loading data
Y <- readRDS("pheno2step")
head(Y)
str(Y)
Y <- droplevels.data.frame(Y[Y$type == "sc",])

# loading the kernels
Ga <- readRDS("Ga")
dim(Ga)
Gd <- readRDS("Gd")
dim(Gd)

# selecting just the hybrids
Ga <- Ga[match(Y$gid, colnames(Ga)), match(Y$gid, colnames(Ga))]
dim(Ga)
all(Y$gid == colnames(Ga))
Gd <- Gd[match(Y$gid, colnames(Gd)), match(Y$gid, colnames(Gd))]
dim(Gd)
all(Y$gid == colnames(Gd))


##################### GS using A model ####################
library(sommer)

# Cross-validation settings 
nfold <- 5 # training and validation sizes
replicates <- 5 # repeat this procedure 5 times

output_A <- data.frame() # create an empty file to storage KPI
GEBVs <- data.frame() # create an empty file to GEBV

for (r in 1:replicates) {
  cat("Processing the replicate", r, "\n")
  ## test sets
  sets <- sample(rep(1:nfold, length.out = n <- dim(Y)[1])) # different orders of the same vector
  
  for (k in 1:nfold) {
    cat("Processing the fold", k, "\n")
    
    # creating the training and testing sets
    itest = which(sets == k)
    itrain = which(sets != k)
    YGS <- Y
    YGS[itest, "dblups"] <- NA # it creates the validation set
    
    # Running the model
    sol <-  mmer(fixed = dblups ~ 1,
                  random = ~vsr(gid, Gu = Ga),
                  rcov = ~ vsr(units),
                  verbose = F,
                  data = YGS)
  
    # predicting and reorganizing the BLUPS
    BLUPS <- predict.mmer(sol, classify = "gid")$pvals
    predicted.blups <- BLUPS[BLUPS$gid %in% YGS$gid[itest], ]
    observed.blues <- Y[Y$gid %in% YGS$gid[itest],]
    observed.blues <- observed.blues[match(predicted.blups$gid, observed.blues$gid),]
      
    output_A <- rbind(output_A, data.frame(
      method = "A",
      rep = r,
      fold = k,
      h2m = round(as.numeric(vpredict(sol, h2 ~ V1 / (V1 + V2))[1]), 2),
      PA = round(cor(predicted.blups[, 3], observed.blues[, 6]), 2),
      BIC = sol$BIC
      )
    )
    
    # retrieving only the genetic deviations
    blups <- data.frame(gid = names(sol$U$`u:gid`$dblups),
                                 GEBV = sol$U$`u:gid`$dblups)
    
    GEBVs <- rbind(GEBVs, blups)
    
  }
}


# estimate the KPI mean
head(output_A)
(resultA <- apply(output_A[,4:6], 2, mean)) 

# and looking at the distribution
library(ggplot2)
ggplot(output_A, aes(x = method, y = PA)) + 
  geom_boxplot(notch = TRUE, outlier.colour="red", outlier.shape=16, outlier.size=2)+
  geom_jitter(position=position_jitter(0.2))

# and the GEBV
dim(GEBVs)
head(GEBVs)

# lets organize the breeding values and pick the mean 
GEBVs <- as.matrix(tapply(GEBVs$GEBV, GEBVs$gid, mean))
head(GEBVs)
GEBVs <- GEBVs[(match(colnames(Ga), rownames(GEBVs))),]
GEBVs <- as.matrix(GEBVs)
colnames(GEBVs) <- "GEBV"
head(GEBVs)
dim(GEBVs)

###################
## G-GLBUP vs RR-BLUP
###################

# function to estimate the variance
vpq <- function(x){
  p <- sum(x)/(length(x)*2)
  pq <- p * (1 - p)
  return(pq)
}

# load the marker dataset
M <- readRDS("M")
M <- M[colnames(Ga),]
dim(M)

# estimathe the marker contribution to the Vg
vg <- 2*sum(apply(M, 2, vpq))

# and finally the marker effects
a <- t(M) %*% solve(Ga) %*% GEBVs / vg

all(rownames(a) == colnames(M))
colnames(a) <- "marker effect"
head(a)
dim(a)

# correlation between GBLUP and RR-BLUP
rrblup.gbv <- M %*% a
colnames(rrblup.gbv) <- "RR-BLUP" 
head(rrblup.gbv)
(pho <- round(cor(GEBVs, rrblup.gbv),2))

joint <- data.frame(gid = Y$gid, GBLUP = GEBVs, RRBLUP = rrblup.gbv)
head(joint)

library(plotly)
(p <- plot_ly(joint, x = ~GEBV, y = ~RR.BLUP, type = 'scatter', mode = 'markers',  text = ~paste('gid: ', gid)))


##################### GS using A + D model + BLUES and weights ####################

# Cross-validation settings 
nfold <- 5 # training and validation sizes
replicates <- 5 # repeat this procedure 5 times

# add a column to model the dominance
Y$gidD <- Y$gid

output_AD <- data.frame() # create an empty file to storage KPI
GEGVs <- data.frame() # create an empty file to GEBV

for (r in 1:replicates) {
  cat("Processing the replicate", r, "\n")
  ## test sets
  sets <- sample(rep(1:nfold, length.out = n <- dim(Y)[1])) # different orders of the same vector
  
  for (k in 1:nfold) {
    cat("Processing the fold", k, "\n")
    
    # creating the training and testing sets
    itest = which(sets == k)
    itrain = which(sets != k)
    YGS <- Y
    YGS[itest, "blues"] <- NA # it creates the validation set
    
    # Running the model
    sol <-  mmer(fixed = blues ~ 1,
                 random = ~vsr(gid, Gu = Ga) + vsr(gidD, Gu = Gd),
                 rcov = ~ vsr(units),
                 verbose = F,
                 weights = w,
                 tolParInv = 1,
                 data = YGS)
    
    # retrieving only the genetic deviations
    blups <- data.frame(gid = names(sol$U$`u:gid`$blues),
                        A = sol$U$`u:gid`$blues,
                        D = sol$U$`u:gidD`$blues,
                        GV = sol$U$`u:gid`$blues + sol$U$`u:gidD`$blues)
    
    predicted.blups <- blups[blups$gid %in% YGS$gid[itest], ]
    observed.blues <- Y[Y$gid %in% YGS$gid[itest],]
    observed.blues <- observed.blues[match(predicted.blups$gid, observed.blues$gid),]
    
    output_AD <- rbind(output_AD, data.frame(
      method = "AD",
      rep = r,
      fold = k,
      h2m = round(as.numeric(vpredict(sol, h2 ~ V1 / (V1 + V2 + V3))[1]), 2),
      h2mAD = round(as.numeric(vpredict(sol, h2 ~ (V1 +V2) / (V1 + V2 + V3))[1]), 2),
      PA = round(cor(predicted.blups[, 3], observed.blues[, 3]), 2),
      BIC = sol$BIC
    )
    )
    
    GEGVs <- rbind(GEGVs, blups)
    
  }
}

# estimate the KPI mean
head(output_AD)
(resultAD <- apply(na.omit(output_AD[,4:7]), 2, mean)) 

# and looking at the distribution
library(ggplot2)
ggplot(output_AD, aes(x = method, y = PA)) + 
  geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 16, outlier.size = 2)+
  geom_jitter(position=position_jitter(0.2))

# and the GEGVs
dim(GEGVs)
head(GEGVs)
# lets organize the breeding values and pick the mean 
GEGVs <- as.matrix(tapply(GEGVs$GV, GEGVs$gid, mean))
head(GEGVs)


################################# GS using OTS ################################

OTS <- readRDS("OTS")
names(OTS)
output_OTS <- data.frame()

for (i in 1:length(OTS)) {
  cat("Processing OTS", i, "\n")
  
  itrain <- OTS[[i]]
  itest <- Y$gid[!(Y$gid) %in% itrain]
  YGS <- Y
  YGS[Y$gid %in% itest, c("blues")] <- NA
  
  # Running the model
  sol <-  mmer(fixed = blues ~ 1,
               random = ~vsr(gid, Gu = Ga),
               rcov = ~ vsr(units),
               verbose = F,
               weights = w,
               data = YGS)
  
  # predicting and reorganizing the BLUPS
  BLUPS <- predict.mmer(sol, classify = "gid")$pvals
  predicted.blups <- BLUPS[BLUPS$gid %in% YGS$gid[itest], ]
  observed.blues <- Y[Y$gid %in% YGS$gid[itest],]
  observed.blues <- observed.blues[match(predicted.blups$gid, observed.blues$gid),]
  
  output_OTS <- rbind(output_OTS, data.frame(
    method = "OTS_A",
    rep = names(OTS)[i],
    fold = names(OTS)[i],
    h2m = round(as.numeric(vpredict(sol, h2 ~ V1 / (V1 + V2))[1]), 2),
    PA = round(cor(predicted.blups[, 3], observed.blues[, 3]), 2),
    BIC = sol$BIC))
  
}

library(knitr)
kable(output_OTS)
library(plotly)
(p2 <- plot_ly(output_OTS, x = ~OTS, y = ~ PA, type = "bar"))

########### the end ################