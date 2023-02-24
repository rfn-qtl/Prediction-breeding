##########################
# AGRO 7075	- Prediction-based Breeding
# Lab 2 - Population genomics
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Feb 16 2023
##########################

# loading the files
genotypes <- readRDS("M")
dim(genotypes)
genotypes[1:9, 1:5]  

# considering just the parents, first 14 rows
genotypes <- genotypes[1:14,]
dim(genotypes)

# These parents represent two sources of germplasm N and P.Thus, lets create a factor for that
metadata <- data.frame(row.names = rownames(genotypes), gid = rownames(genotypes))
metadata$Group <- as.factor(rep(c("P", "N"), each = 7))
metadata

# PCA analysis
#BiocManager::install("PCAtools")
library(PCAtools)
p <- pca(t(genotypes), metadata = metadata, removeVar = 0.0)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
PCAtools::biplot(p, colby = 'Group', colkey = c(N = 'forestgreen', P = 'Red'), shape = 'Group', shapekey = c(N = 10, P = 21), legendPosition = 'right')

################################ 
# hierarchical clusters
###############################
hc <- hclust(dist(genotypes))
# very simple dendrogram
plot(hc)

# dendrogram objects
hcd <- as.dendrogram(hc)
# alternative way to get a dendrogram
plot(hcd)

# jaccard
#install.packages("vegan")
library("vegan")
?vegdist
jac <- vegdist(genotypes, method = "jaccard", binary = FALSE, diag = TRUE, upper = TRUE, na.rm = FALSE)
plot(hclust(jac))

#another way to identify clusters is by
wardc <- hclust(dist(genotypes), method = "ward.D2")
library(factoextra)
fviz_dend(wardc, k = 4)

# clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

# or we can use this function to determine the optimal number of clusters
set.seed(123)
genotypes2 <- scale(genotypes)
fviz_nbclust(genotypes2, kmeans, method = "silhouette")
#fviz_nbclust(genotypes2, kmeans, method = "gap_stat")
#fviz_nbclust(genotypes2, kmeans, method = "wss")

# K-Means Clustering
k <- kmeans(genotypes2, centers = 4, nstart = 25)
fviz_cluster(k, data = genotypes2)
k$size
sum(k$size)
(clusters <- data.frame(sort(k$cluster)))

# DAPC
library(adegenet)
grp <- find.clusters(genotypes2, n.pca = 14, n.clust = length(k$size))
dapc1 <- dapc(genotypes2, grp$grp, n.pca = 14, n.da = 2)
scatter(dapc1, scree.da = T, posi.da="bottomleft", posi.pca = 'topright', scree.pca = T, cell=1.5, cex=2, bg="white", cstar=0)
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:4), lab="", ncol=2, xlab="individuals")
assignplot(dapc1)
dev.off()

##########################
# signatures of selection 
##########################
# FST method
M <- t(genotypes)
dim(M)
head(M[,1:14])

# calculate frequencies
X <- matrix(NA, nrow(M), 2)
colnames(X) <- c("P","N")

# P
X[,1] <- apply(M[, 1:7], 1, function(x) sum(x)/(length(x)*2))

# N
X[,2] <- apply(M[, 8:14], 1, function(x) sum(x)/(length(x)*2))


# FST
meansX <- rowMeans(X) # average allele frequency across populations
alleleVar <- meansX * (1 - meansX) # p*q variance
meanDevX <- X - meansX # deviation of each population from mean
FST <- meanDevX^2 / alleleVar # deviation squared divided by var
rownames(FST) <- colnames(genotypes)
head(FST)

# creting the graphs
# install.packages("lokern")
library(lokern)
smoothP = lokerns(FST[,1], n.out = 200)
smoothN = lokerns(FST[,2], n.out = 200)

par(mfrow=c(1,2))
plot(FST[,1], type = "l",xlab = "SNP",ylab = "Fst",col = "gray", main = "P")
lines(smoothP$est, type = "l", col = "red",lwd = 1)
plot(FST[,2], type = "l",xlab = "SNP",ylab = "Fst",col = "gray", main = "N")
lines(smoothN$est,type = "l", col = "red",lwd = 1)


#######################
# opening the black box
#######################

# Quality control
#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")
#devtools::install_github(repo = 'italo-granato/snpReady', ref = 'dev')
library(snpReady)

args(popgen)

# datasets 
genotypes[1:14,1:6]
metadata

#estimating pop gen parameters, for individuals, markers and subgroups
population <- popgen(genotypes, subgroups = metadata$Group, plot = F)

# whole population
population$whole$Markers[1:14,]
population$whole$Genotypes
population$whole$Population
population$whole$Variability

# within groups (P example)
population$bygroup$P$Markers[1:14,]
population$bygroup$P$Genotypes
population$bygroup$P$Population
population$bygroup$P$Variability
population$bygroup$P$exclusive
population$bygroup$P$fixed


# Fst
population$bygroup$F.stats
population$bygroup$F.stats$Genotypes
population$bygroup$F.stats$Markers[1:14,]

# Also, there are many graphs in your folder...


######### Finaly, the LD decay #####################
hapmap <- readRDS("hapmap")
hapmap.1 <- hapmap[hapmap$chrom == 1,]
# just a sample to avoid a long time waiting
sns.sampled <- sort(sample(1:nrow(hapmap.1), 1000))
hapmap.1 <- hapmap.1[sns.sampled,]
M.1 <- readRDS("M")[1:14, hapmap.1$rs]
dim(M.1); dim(hapmap.1)

library(SNPRelate)
map1 <- data.frame(snp.id = as.integer(1:dim(hapmap.1)[1]), 
                   snp.rs.id = hapmap.1$rs, 
                   chr = as.integer(hapmap.1$chrom), 
                   pos = as.integer(hapmap.1$pos))
head(map1)

# creating the GDS file
snpgdsCreateGeno(gds.fn = "./toy.gds",                              # gds filename
                 genmat = M.1,                                        # markers matrix
                 sample.id = rownames(M.1),                           # individual names
                 snp.id = map1$snp.id,                               # snp id (deve ser único)
                 snp.rs.id = map1$snp.rs.id,                         # snp name (pode não ser único)
                 snp.chromosome = map1$chr,                          # cromossomo
                 snp.position = map1$pos,                            # posição dentro do cromossomo (bp)
                 #snp.allele = map$allele,                           # alelos (ref / alt)
                 snpfirstdim = FALSE)                                # argumento para matriz n x s


# Loading gds
genofile <- snpgdsOpen(filename = "./toy.gds")
genofile

# estimating the pairwised LD
LDs <- abs(as.matrix(snpgdsLDMat(genofile, method = "corr", slide = 0, num.thread = 8)$LD))
dim(LDs)
rownames(LDs) <- colnames(LDs) <- hapmap.1$rs
LDs[1:5, 1:5]
saveRDS(LDs, "LDs")
# close GDS
snpgdsClose(genofile)

#install.packages("reshape2")
library(reshape2)
lds <- melt(LDs)
head(lds)
dim(lds)

naive_dist <- function(A, B){
  result = matrix(ncol = length(B), nrow = length(A))
  for (i in 1:length(A))
    for (j in 1:length(B))
      result[i,j] = abs(A[i] - B[j])
  result
}

# estimating the distance between markers
distance <- naive_dist(as.matrix(hapmap.1$pos), as.matrix(hapmap.1$pos)) 
rownames(distance) <- colnames(distance) <- hapmap.1$rs
distance2 <-  melt(distance)
head(distance2)
dim(distance2)
pairLD <- data.frame(marker1 = lds$Var1, marker2 = lds$Var2, r2 = lds$value, dist = distance2$value)
head(pairLD)
# eliminating the markers with itself
pairLD <- pairLD[pairLD$dist != 0,]
# then, subset only markers with significant LD 
dt.cor <- function(x, n){
  dt.calc = dt(sqrt(x) / sqrt((1-x)) * sqrt((n-2)), n) 
  return(dt.calc)
  }
pairLD$p <- dt.cor(pairLD$r2, nrow(M.1))
head(pairLD)
pairLD <- pairLD[pairLD$p < 0.001,]
dim(pairLD)

# average LD and dist
apply(na.omit(pairLD[,3:4]), 2, mean)

# and finally the famous graph
library(ggplot2)
ld.plot <- ggplot(data = pairLD, aes(y = r2, x = dist)) +  
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ x^2) +
  labs(x = 'Distance (bp)', y = expression(r2~(dist^{2}))) +
  xlim(0, max(na.omit(pairLD$dist))) +
  ylim(0, 1) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("LD_decay.pdf", plot = ld.plot, dpi = 300, width = 55, height = 30, units = "cm", bg = "white")

####### The end ###########