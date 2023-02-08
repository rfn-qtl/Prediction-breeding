########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 4 - Genomic kinship
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: November 30 2022
#######################################

# read marker data
M  <- readRDS("M")
dim(M)
head(M[,1:6])
tail(M[,1:6])

#source("https://bioconductor.org/biocLite.R")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("impute")
#devtools::install_github(repo = 'italo-granato/snpReady', ref = 'dev')
library(snpReady)
library(ASRgenomics)
# creates the K matrix - genomic relatioship matrix 
args(G.matrix)

G <- snpReady::G.matrix(M, method = "VanRaden", format = "wide", plot = F)

############################### additive #######################
Ga <- G$Ga
dim(Ga)

# let`s load our A matrix and check if they are the same population
A <- readRDS("A")
G2A <- match.G2A(A = A, G = Ga, clean = TRUE, ord = F, mism = TRUE,
                 RMdiff = TRUE)
Ga <- G2A$Gclean
dim(Ga)
Ga[1:4, 1:4]
Ga[52:55, 52:55]

# rescale the matrix to match with the expected values
library(scales)
rescale.G <- function(x){
  out <- scales::rescale(c(x), to = c(0, 2)) 
  out <- matrix(out, nrow = nrow(x), byrow = T)
  rownames(out) <- colnames(out) <- colnames(x)
return(out)
  }

Ga <- rescale.G(Ga)
dim(Ga)
Ga[1:4, 1:4]
Ga[52:55, 52:55]

# let's check the Ga matrix
check_Ga <- kinship.diagnostics(K = Ga, diagonal.thr.small = 0.8,
                               diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_Ga$plot.diag
check_Ga$plot.offdiag

########################### dominance kinship #######################
Gd <- G$Gd[colnames(Ga), colnames(Ga)]
dim(Gd)
Gd[1:4, 1:4]
Ga[52:55, 52:55]

Gd <- rescale.G(Gd)
dim(Gd)
Gd[1:4, 1:4]
Ga[52:55, 52:55]

# let's check the Gd matrix
check_Gd <- kinship.diagnostics(K = Gd, diagonal.thr.small = 0.8,
                                diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_Gd$plot.diag
check_Gd$plot.offdiag

#################### # graphs analysis ##############################

# svd decomposition for Ga and Gd
Ga_pca <- kinship.pca(K = Ga, ncp = 14, label = T, ellipses = T)
Ga_pca$plot.pca
Ga_pca$plot.scree

Gd_pca <- kinship.pca(K = Gd, ncp = 14, label = T, ellipses = T)
Gd_pca$plot.pca
Gd_pca$plot.scree

############### heatmaps ##############

kinship.heatmap(K = Ga, dendrogram = TRUE, row.label = T,
                col.label = T)

kinship.heatmap(K = Gd, dendrogram = TRUE, row.label = T,
                col.label = T)

#saving the final kinship versions
saveRDS(Ga, "Ga")
saveRDS(Gd, "Gd")

############################ Inbreeding and Ne ########################
# Inbreeding per individual - considering just the parents
(Fi <- round(diag(Ga)-1, 2)[1:14])

# Effective size per individual
(Ne.i <- 1/(2*Fi))

# Ne of population
(Ne.pop <- sum(Ne.i))

# Population endogamy rate
(Fi.pop <- 1/(2*Ne.pop))

############# the end ##############