########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 3 - Pedigree
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: November 30 2022
#######################################

# install.packages("pedigreemm")
library(pedigreemm) #load pedigreemm package

# read in the pedigree data
ped  <- read.table("ped.txt", header = TRUE)
head(ped, 15)
tail(ped)
str(ped)

#editPed function orders the pedigre from oldest to newest	
args(editPed)
ped2 <- editPed(ped$male, ped$female, ped$gid)
head(ped2) # the last col is generation
tail(ped2)

#constructs the pedigree object
args(pedigree)
ped3 <- pedigree(ped2$sire, ped2$dam, ped2$label)
head(ped3)
tail(ped3)

# creates the A matrix (relationship matrix, which means 2 times de kinship matrix)
A <- getA(ped3)
A <- as.matrix(A) 
dim(A)
A[1:14, 1:14] # the first seven parents
A[47:55, 47:55] # the last hybrids

# let`s save or A matrix`
saveRDS(A, "A")

####################
# graphs analysis
####################

# svd decomposition - by individuals
svdG <- svd(A, nu = ncol(A), nv = nrow(A))
plot(cumsum((svdG$d[1:ncol(A)])^2/sum(svdG$d^2)), ylab = "proportion accumulated", xlab = "number of individuals", col = "red")

# obtainig the eigenvectors and eigenvalues
pcsG <- A %*% svdG$v
rownames(pcsG) <- colnames(pcsG) <- rownames(A)
dim(pcsG)
pcsG[1:14,1:5]

# PCA graphs
# proportion explained by the first componentes
axispcs <- paste((round(svdG$d[1:ncol(A)]^2/sum(svdG$d^2)*100))[1:3], "%", sep = "")
axispcs

# 3D graph 
library(scatterplot3d)
scatterplot3d(pcsG[,1], pcsG[,2], pcsG[,3], xlab = axispcs[1], ylab = axispcs[2], zlab = axispcs[3], axis = TRUE, color = "red", highlight.3d = FALSE, box = TRUE, angle = 50)

# 2D graph
par(mfrow = c(1,2))
plot(x = pcsG[,1], y = pcsG[,2], xlab = axispcs[1], ylab = axispcs[2], col = "red", main = "PC 1 vs PC 2")
plot(x = pcsG[,1], y = pcsG[,3], xlab = axispcs[1], ylab = axispcs[3], col = "blue", main = "PC 1 vs PC 3")
dev.off()

##############
# heatmaps
##############
#install.packages("superheat")
library(superheat)
args(superheat)

superheat(A, pretty.order.rows = T, pretty.order.cols = T, col.dendrogram = T, clustering.method = "kmeans", 
          dist.method = "euclidean",  bottom.label.text.size = 2, left.label.text.size = 2, legend.text.size = 5)

# saving the graph for papers
png("heatmapA.png", width = 8, height = 8, res = 400, units = "in")
superheat(A, pretty.order.rows = T, pretty.order.cols = T, col.dendrogram = T, clustering.method = "kmeans", 
          dist.method = "euclidean",  bottom.label.text.size = 2, left.label.text.size = 2, legend.text.size = 5)
dev.off()

# running my first BLUP-A with pedigreem
pheno <- readRDS("pheno")
head(pheno)

?pedigreemm()

fit <- pedigreemm(SDM ~ N + rep + row + col + (1|gid) + (1|gid:N), 
                  pedigree = list(gid = ped3), 
                  data = pheno)

str(summary(fit))
summary(fit)$logLik
summary(fit)$coefficients
anova(fit)
summary(fit)$varcor
summary(fit)$AICtab

# heritability
(data.frame(summary(fit)$varcor)[2,4] / (data.frame(summary(fit)$varcor)[2,4] + 
                                          data.frame(summary(fit)$varcor)[3,4] +
                                          data.frame(summary(fit)$varcor)[1,4]))

fixef(fit) # fixed effects
(Random.fit <- ranef(fit)) # random effects
hist(resid(fit), xlab = "Residuals from fit", col = "red", main = NULL)
(Inbreeding <- inbreeding(ped3))

##### the end #########