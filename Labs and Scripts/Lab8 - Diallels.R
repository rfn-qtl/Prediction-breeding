########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 8 - Diallels
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: March 21 2023
#######################################

# loading data
pheno <- readRDS("pheno")
head(pheno)
str(pheno)

library(sommer)

#################### Full diallel designs (use of the overlay) - parents non-related ####################

modFD <- mmec(SDM ~ N + rep,
              random = ~vsc(isc(overlay(female, male))) + gid, 
              nIters = 3,
              data = pheno, 
              verbose = FALSE)

with(pheno, overlay(female, male, sparse = FALSE))

(suma <- summary(modFD)$varcomp)
Vgca <- sum(suma[1,1])
Vsca <- suma[2,1]
Ve <- suma[3,1]
Va <- 4*Vgca
Vd <- 4*Vsca
Vg <- Va + Vd
(H2.FD <- Vg / (Vg + (Ve)))
(h2.FD <- Va / (Vg + (Ve)))

# CGA and SCA
(CAE.FD <- data.frame(gid = rownames(modFD$u),
                     effect = modFD$u,
                     type = c(rep("CGA", length(unique(c(pheno$female, pheno$male)))), 
                              rep("SCA", length(unique(pheno$gid)))))
)


############################### Half diallel or NCII - parents non-related ###########################
# only hybrids
phenoSC <-  droplevels.data.frame(pheno[pheno$type == "sc",])

modHD <- mmec(SDM ~ N + rep,
              random= ~ female + male + gid,
              rcov = ~units, 
              nIters = 3,
              data = phenoSC, 
              verbose = FALSE)

(suma <- summary(modHD)$varcomp)
Vgca <- sum(suma[1:2,1])
Vsca <- suma[3,1]
Ve <- suma[4,1]
Va <- 4*Vgca
Vd <- 4*Vsca
Vg <- Va + Vd
(H2.HD <- Vg / (Vg + (Ve)))
(h2.HD <- Va / (Vg + (Ve)))

# CGA and SCA
(CAE.HD <- data.frame(gid = rownames(modHD$u)[modHD$u != 0],
                     effect = modHD$u[modHD$u != 0],
                     type = c(rep("CGA", length(unique(c(pheno$female, pheno$male)))), 
                              rep("SCA", length(unique(pheno$gid)) - length(unique(c(pheno$female, pheno$male)))))))

# Comparing the models - LRT, AIC, and BIC
lrt <- anova(modFD, modHD)


############################### Half diallel or NCII - parents related ###########################

#loading kinships
Ga <- readRDS("Ga")
dim(Ga)
# Ga for females
Ga[1:7, 1:7]
Ga.f <- Ga[1:7, 1:7]
# Ga for males
Ga[8:14, 8:14]
Ga.m <- Ga[8:14, 8:14]

# dominance
Gd <- readRDS("Gd")
Gd <- Gd[15:ncol(Gd), 15:ncol(Gd)]
dim(Gd)

# running the model - pay attnetion tha now is via mmer
modHDG <- mmer(SDM ~ N + rep,
              random= ~ vsr(female, Gu = Ga.f) + vsr(male, Gu = Ga.m) + vsr(gid, Gu = Gd),
              rcov = ~vsr(units), 
              nIters = 3,
              data = phenoSC, 
              verbose = FALSE)

(suma <- summary(modHDG)$varcomp)
Vgca <- sum(suma[1:2,1])
Vsca <- suma[3,1]
Ve <- suma[4,1]
Va <- 4*Vgca
Vd <- 4*Vsca
Vg <- Va + Vd
(H2.HDG <- Vg / (Vg + (Ve)))
(h2.HDG <- Va / (Vg + (Ve)))

# CGA and SCA
(CAE.HDG <- data.frame(gid = names(c(modHDG$U$`u:female`$SDM, modHDG$U$`u:male`$SDM, modHDG$U$`u:gid`$SDM)),
                     effect = c(modHDG$U$`u:female`$SDM, modHDG$U$`u:male`$SDM, modHDG$U$`u:gid`$SDM),
                     type = c(rep("CGA", length(unique(c(pheno$female, pheno$male)))), 
                              rep("SCA", length(unique(pheno$gid)) - length(unique(c(pheno$female, pheno$male)))))))

# Comparing the models - BIC
modHD$BIC; modHDG$BIC 

# correlation between the two models, with and without kinship
cor(CAE.HD[, 2], CAE.HDG[, 2])
plot(CAE.HD[, 2], CAE.HDG[, 2], ylab = "GEBV", xlab = "PBV")
abline(lm(CAE.HDG[, 2] ~ CAE.HD[, 2]), col = "red")


############# Building heterotic pools (the dataset does not help, but let's how it would be...) ###############

# For that we need to run a full diallel using G matrices
Ga <- readRDS("Ga")
Gd.mf <- kronecker(Ga[1:14, 1:14], Ga[1:14, 1:14])
dim(Gd.mf)
Gd.mf[1:5, 1:5]

# colnames for Gd.mf
out <- c()
aux <- expand.grid(1:14, 1:14)
for(i in 1:nrow(aux)){
if(aux[i,1] == aux[i,2]){out <- c(out, aux[i,1])}
if(aux[i,1] != aux[i,2]){out <- c(out, paste0(aux[i,1], aux[i,2]))}
}

aux2 <- cbind(out, aux)
colnames(aux2) <- c("gid", "female", "male") 
head(aux2)

colnames(Gd.mf) <- rownames(Gd.mf) <- out
Gd.mf[1:5, 1:5]
pheno$gid %in% colnames(Gd.mf)

modFDG <- mmer(SDM ~ N + rep,
               random= ~ vsr(female, Gu = Ga) + vsr(male, Gu = Ga) + vsr(gid, Gu = Gd.mf),
               rcov = ~vsr(units), 
               nIters = 3,
               data = pheno, 
               verbose = FALSE)

# merge datasets
SCA <- data.frame(gid = names(modFDG$U$`u:gid`$SDM), SCA = modFDG$U$`u:gid`$SDM)
SCA <- merge(SCA, aux2)
head(SCA)
SCA <- SCA[,2:4]

library(reshape2)
d <- reshape(SCA, idvar = "female", timevar = "male", direction = "wide")
dim(d)
d <- d[,-1]
colnames(d) <- gsub("SCA.", "", colnames(d), fixed = T)
d <- d[, match(colnames(d), sort(as.numeric(colnames(d))) )]
rownames(d) <- colnames(d)
diag(d) <- 0
d[1:5, 1:5]

# clustering in heterotic pools
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

k <- kmeans(scale(d), centers = 2, nstart = 25)
fviz_cluster(k, data = scale(d), main = "Heterotic Pools")
k$size
sum(k$size)
(heterotic.pools <- data.frame(sort(k$cluster)))


# Selecting testers
GCA <- (modFDG$U$`u:female`$SDM + modFDG$U$`u:male`$SDM) / 2

testers <- data.frame()
for(i in 1:2){
HP <- GCA[rownames(heterotic.pools)[heterotic.pools$sort.k.cluster. == i]]
parent <- names(HP[(HP == max(HP))])
testers <- rbind(testers,
                 data.frame(
  HP = i,
  Tester = parent
))}

testers

############################ the end ##########################