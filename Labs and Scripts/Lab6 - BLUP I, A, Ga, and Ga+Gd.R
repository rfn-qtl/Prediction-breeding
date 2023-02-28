########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 6 - BLUP I, A, Ga, and Ga+Gd
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Feb 28 2023
#######################################

#loading phenotypes
pheno <- readRDS("pheno")
head(pheno)
dim(pheno)
str(pheno)

# loading A matrix
A <- readRDS("A")
dim(A)

# Loading Ga
Ga <- readRDS("Ga")
all(pheno$gid %in% rownames(Ga))
dim(Ga)

# Loading Ga and Gd
Gd <- readRDS("Gd")
all(pheno$gid %in% rownames(Gd))
dim(Ga)

# GxE Kernells
(NL <- diag(1, 2)) # two env.: low and ideal N
colnames(NL) <- rownames(NL) <- c("ideal", "low")
NL

# a short example of kronecker
kronecker(NL, Ga[1:3, 1:3])

A.ge <- kronecker(NL, A)
dim(A.ge)
A.ge[1:5, 1:5]

Ga.ge <- kronecker(NL, Ga)
dim(Ga.ge)
Ga.ge[1:5,1:5]

Gd.ge <- kronecker(NL, Gd)
dim(Gd.ge)
Gd.ge[1:5,1:5]

# giving names for columns and rows
GxN <- paste0(colnames(Ga), rep(c("ideal", "low"), each = length(colnames(Ga))))
colnames(A.ge) <- rownames(A.ge) <- colnames(Ga.ge) <- rownames(Ga.ge) <- colnames(Gd.ge) <- rownames(Gd.ge) <- GxN
dim(Ga.ge)
Ga.ge[1:5,1:5]

# let's add two new col in the data.frame
pheno$GN <- as.factor(paste0(pheno$gid, pheno$N)) # helps to model the GxE
pheno$GNd <- as.factor(paste0(pheno$gid, pheno$N)) # helps to model the GxE + D
pheno$gidD <- as.factor(pheno$gid) # helps to model the Gd

##############################################
# comparing the models with different kernels
##############################################

library(sommer)

########################################## model I (just the identity matrix) ##########################
fitI <- mmer(fixed = SDM ~ 1 + N,
             random = ~gid + rep + gid:N,
             rcov = ~ vsr(units),
             data = pheno)
# varcomp
unlist(fitI$sigma)

########################################## model A (pedigree matrix) ##########################
fitA <- mmer(fixed = SDM ~ 1 + N,
             random = ~vsr(gid, Gu = A) + rep + vsr(GN, Gu = A.ge),
             rcov = ~ vsr(units),
             data = pheno)
# varcomp
unlist(fitA$sigma)

########################################## model Ga (Ga matrix) ##########################
fitGa <- mmer(fixed = SDM ~ 1 + N,
             random = ~vsr(gid, Gu = Ga) + rep + vsr(GN, Gu = Ga.ge),
             rcov = ~ vsr(units),
             data = pheno)
# varcomp
unlist(fitGa$sigma)
# BLUPS
gebvs <- predict.mmer(fitGa, classify = "gid")$pvals
head(gebvs, 10)

########################################## model Ga + Gd (Ga and Gd matrices) ##########################
fitGad <- mmer(fixed = SDM ~ 1 + N,
              random = ~vsr(gid, Gu = Ga) + rep + vsr(gidD, Gu = Gd) + vsr(GN, Gu = Ga.ge) + vsr(GNd, Gu = Gd.ge),
              rcov = ~ vsr(units),
              data = pheno)
# varcomp
unlist(fitGad$sigma)

########################################## heritabilities #################################################
(h2 <- round(data.frame(I = as.numeric(vpredict(fitI, h2 ~ V1 / (V1 + V3 + V4))[1]),
                    A = as.numeric(vpredict(fitA, h2 ~ V1 / (V1 + V3 + V4))[1]),
                    Ga = as.numeric(vpredict(fitGa, h2 ~ V1 / (V1 + V3 + V4))[1]),
                    Gad = as.numeric(vpredict(fitGad, h2 ~ V1 / (V1 + V3 + V4 + V5 + V6))[1])), 2))

###########################################    BLUPS      ##################################################
BLUPS <- data.frame(gid = names(fitA$U$`u:gid`[[1]]),
                    I = fitI$U$gid[[1]],
                    A = fitA$U$`u:gid`[[1]],
                    Ga = fitGa$U$`u:gid`[[1]],
                    Gad = fitGad$U$`u:gid`[[1]] + fitGad$U$`u:gidD`[[1]])

head(BLUPS)

# correlation between BV surrogates 
require(PerformanceAnalytics)
# correlarion betwwen BLUPS ontaines from differente methods
chart.Correlation(BLUPS[,2:5], histogram = TRUE, pch = 1, method = "pearson")


############################### selecting the best genotypes by I or Ga #############################
# first, define a quantile
(trshldI <- quantile(BLUPS$I, probs = 0.90))
(trshldGa <- quantile(BLUPS$Ga, probs = 0.90))

# Then, identify the best ones
(selectedI <- BLUPS[BLUPS$I >= trshldI,])
(selectedGa <- BLUPS[BLUPS$Ga >= trshldGa,])

# and estimate the response to selection
(RSI <- round(mean(selectedI$I),4))
(RSGa <- round(mean(selectedGa$Ga),4))

###################### Inbreeding and Ne in the next generation ##################################
# Fi per individual - considering just the parents
(FiI <- round(diag(Ga), 4)[match(selectedI$gid, colnames(Ga))])
(FiGa <- round(diag(Ga), 4)[match(selectedGa$gid, colnames(Ga))])

# Effective size per individual and population
(Ne.iI <- 1 / (2*FiI)); (Ne.iGa <- 1 / (2*FiGa))
(Ne.popI <- sum(Ne.iI)); (Ne.popGa <- sum(Ne.iGa))

# Population inbreeding rate
(Fst.popI <- 1/(2*Ne.popI)); (Fst.popGa <- 1/(2*Ne.popGa))

########## the end ####################