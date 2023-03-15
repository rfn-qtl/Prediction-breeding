########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 7 - Adjusting phenotypes for two-step analysis and OTS
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Dec 5 2022
#######################################

########################################### 
# adjusting phenotypes for two-set analysis 
##########################################

#loading phenotypes
pheno <- readRDS("pheno")
head(pheno)
dim(pheno)
str(pheno)

#modeling
library(sommer)

# using only the classical experimental design, and gid as random
fitR <- mmer(fixed = SDM ~ 1 + N,
             random = ~gid + rep + gid:N,
             rcov = ~ vsr(units),
             data = pheno)

# predicting the BLUPs - main effect
BLUPs <- predict.mmer(object = fitR, classify = "gid")$pvals

# reliability
r <- 1 - BLUPs$standard.error^2 / as.numeric(fitR$sigma$gid)

# deregressed BLUPS
BLUPs$dBLUPS <- BLUPs$predicted.value / r
head(BLUPs)

# using only the classical experimental design, and gid as fixed
fitF <- mmer(fixed = SDM ~ 1 + N + gid + gid:N,
             random = ~rep,
             rcov = ~ vsr(units),
             data = pheno)

# predicting the BLUEs - main effect
BLUEs <- predict.mmer(object = fitF, classify = "gid")$pvals
head(BLUEs)

# weights for ID's - adjust residual for further analysis
w <- diag(fitR$PevU$gid$SDM)

pheno2step <- data.frame(gid = BLUPs$gid,
                         trait = fitF$terms$response[[1]][1],
                         blues = BLUEs$predicted.value,
                         w = w,
                         blups = BLUPs$predicted.value,
                         dblups = BLUPs$dBLUPS)

aux <- pheno[!duplicated(pheno$gid), 2:5]
pheno2step <- merge(pheno2step, aux)
pheno2step <- pheno2step[order(pheno2step$gid),]
dim(pheno2step)
head(pheno2step)
cor(pheno2step[,c(3,5:6)])
saveRDS(pheno2step, "pheno2step")


###########################################
# OTS - reducing the number of individuals
###########################################

# selecting just the hybrids
phenoGS <- pheno2step[as.numeric(pheno2step$gid) >14,]

# loading the marker file
Ga <- readRDS("Ga")
dim(Ga)
rownames(Ga)

# and selecting only the markers concerning SC
Ga <- Ga[match(phenoGS$gid, rownames(Ga)), match(phenoGS$gid, rownames(Ga))]
dim(Ga)
all(phenoGS$gid == rownames(Ga))
all(phenoGS$gid == colnames(Ga))


# PCA and svd analysis
svdGa <- svd(Ga, nu = ncol(Ga), nv = nrow(Ga))
plot(cumsum((svdGa$d[1:ncol(Ga)])^2/sum(svdGa$d^2)), col = "red")
pcsGa <- Ga %*% svdGa$v
rownames(pcsGa) <- rownames(Ga)
dim(pcsGa)
# number od PCs that explain more than 90% of the variance
sum(cumsum((svdGa$d[1:ncol(Ga)])^2/sum(svdGa$d^2)) < 0.90)

# based on the Miztal results, 98% of Va in enough to represent the whole dataset
n_of_pc <- function(K, proportion){
  svdK <- svd(K, nu = ncol(K), nv = nrow(K))
  cs <- sum(cumsum((svdK$d[1:ncol(K)])^2/sum(svdK$d^2)) < proportion)
  return(cs)
}

#however, lets try different proportions
proportions <- c(.90, .95, .98)
size <- c()

for(i in 1:length(proportions)){
  size <- c(size, n_of_pc(Ga, proportions[i]))
}

# different sizes of TS
size

# optimizing training sets - OTS
#install.packages("STPGA")
library(STPGA)

OTS <- list()

system.time(
for(i in 1:length(size)) {
  OTS[[i]] <- GenAlgForSubsetSelectionNoTest(P = Ga, ntoselect = size[i], npop = 5, nelite = 5, mutprob = .8, niterations = 100, plotiters = TRUE, lambda = 1e-5, errorstat = "PEVMEAN", mc.cores = 8)[[1]]
}
)

names(OTS) <- paste("OTS", size, sep = "")
OTS

saveRDS(OTS, "OTS")

########## the end ####################