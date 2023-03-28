########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 9 - GWAS
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Dec 5 2022
#######################################

# loading phenotypes
Y <- readRDS("pheno2step")
head(Y)
dim(Y)

# loading the marker file
M <- readRDS("M")
dim(M)

# selecting just the hybrids
sc.names <- Y$gid[Y$type == "sc"]
Y <- Y[match(sc.names, Y$gid),]
M <- M[match(sc.names, rownames(M)),]
dim(M)
head(M[,1:6])
tail(M[,1:6])
all(Y$gid == rownames(M))

# loading the hapmap
hapmap <- readRDS("hapmap")
head(hapmap)
colnames(hapmap) <- c("SNP", "Chrm", "Pos")
dim(hapmap)

#double check
all(colnames(M) == hapmap$SNP)

# add a new collum (gid)
M2 <- data.frame(Taxa = rownames(M), M)
head(M2[,1:6])
dim(M2)

######################################## PCA based on Ga ##################################
Ga <- readRDS("Ga")
Ga <- Ga[rownames(M), rownames(M)]
Ga[1:5, 1:5]
dim(Ga)

library(PCAtools)
metadata <- data.frame(row.names = colnames(Ga))
p <- pca(Ga, metadata = metadata, removeVar = 0.0)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p)
PCA <- p$loadings
dim(PCA)
head(PCA[,1:3])

################### FarmCPUpp ##############
library(bigmemory)
library(biganalytics)
library(compiler)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

# preparing final files
head(Y)
myY <- Y[,c(1, 6)] # using dBLUPs because there is no space to add weights in this package
head(myY)

myGD <- M2
head(myGD)

myGM <- hapmap
head(hapmap)

myCV <- PCA

################################################### GWAS models ####################################
args(FarmCPU)

model.K <- FarmCPU(Y = myY, GD = myGD, GM = myGM, cutOff = 0.05, p.threshold = 0.05/nrow(myGM), MAF.calculate = T)
head(model.K$GWAS)

model.PC1 <- FarmCPU(Y = myY, GD = myGD, GM = myGM, CV = myCV[,1], cutOff = 0.05, p.threshold = 0.05/nrow(myGM), MAF.calculate = T)

model.PC2 <- FarmCPU(Y = myY, GD = myGD, GM = myGM, CV = myCV[,1:2], cutOff = 0.05, p.threshold = 0.05/nrow(myGM), MAF.calculate = T)

############################   identifying the real p.threshold value for our data   ########################
args(FarmCPU.P.Threshold)
p.adj <- FarmCPU.P.Threshold(Y = myY, GD = myGD, GM = myGM, trait = "blues", theRep = 10)
p.adj #the real p.threshold value for our data
(p.Bonf <- 0.05/ncol(M)) # p by Bonferroni

# Thus, it is the newest threshold
(thr.p.adj <- -log10(p.adj))
-log10(p.Bonf) #threshold by Bonferroni

############################    running again but using the newest threshold   ###########################
model.p.adj <- FarmCPU(Y = myY, GD = myGD, GM = myGM, p.threshold = p.adj, threshold.output = p.adj, MAF.calculate = TRUE, cutOff = c((p.adj*ncol(M)), 0.05))

# combining all outputs via different models
gwas.all <- data.frame(model.K$GWAS[,1:4], model.PC1$GWAS[, 4], model.PC2$GWAS[, 4], model.p.adj$GWAS[, 4])
colnames(gwas.all)[4:7] <- c("K", "1PC", "2PC", "p.adj")
head(gwas.all)

# correlation amog marker p-values via different models
cor(gwas.all[,4:7], use = "pairwise.complete.obs")


##############################              circular plot            ###########################
library(CMplot)

CMplot(gwas.all, type = "p", plot.type = c("c", "q"),
       cir.chr.h=1.5, cir.legend = F,
       signal.line=1,
       width=10, height=10,
       chr.labels = paste("chrmo", " ", unique(hapmap$Chrm), sep=""),
       threshold=c(p.adj, 0.05/ncol(M)),threshold.col = c('black','red'),
       multracks = TRUE, chr.den.col = NULL, file.output = TRUE,
       signal.col = NULL, outward = TRUE, amplify = T, ylab = colnames(gwas.all[4:7]))


######################## significance, MAF, and allele substitution effect ####################
head(model.p.adj$GWAS)
(qtl <- model.p.adj$GWAS[which(model.p.adj$GWAS$P.value < p.adj),]) # identifying the markers

# boxplots = phenotypes ~ number of copies
par(mfrow = c(nrow(qtl), 1))
for(i in 1:length(qtl$SNP)){
  boxplot(Y$blues ~ M[, qtl$SNP[i]], main = qtl$SNP[i], xlab = "Number of copies", ylab = "SDM", col = "red")
  }

# number of copies per genotype class
library(knitr)
for(i in 1:length(qtl$SNP)) {
  print(kable(qtl$SNP[i]))
  print(kable(table(M[,qtl$SNP[i]])))
}


# heritability and R2
qtl$Va <- round(2 * qtl$maf * (1 - qtl$maf) * qtl$effect^2, 3)
qtl$ha <- round(qtl$Va / var(Y$blues), 2)
qtl$Ac <- round(sqrt(qtl$ha), 2)
qtl$R2 <- qtl$ha^2
kable(qtl)

# when you have more than one significant marker, we need to run a multiple linear regression
# for that, use all the significant marker together
fit.lm.data <- data.frame(blues = Y[,3], M[,qtl$SNP])
colnames(fit.lm.data)[2:ncol(fit.lm.data)] <- qtl$SNP
head(fit.lm.data)
str(fit.lm.data)
fit.lm.data <- droplevels.data.frame(fit.lm.data)
str(fit.lm.data)

# finally, run the multiple regression
fi.lm <- lm(blues ~ ., fit.lm.data)
summary(fi.lm)
fi.lm$coefficients
anova(fi.lm)

# estimating breeding values based on significant markers
bv <- matrix(M[, qtl$SNP]) %*% qtl$effect
cor(bv, Y$blues)

######################################## LD.plot for annotation ##################################
library(SNPRelate)

# Chrm by Chrm
head(hapmap)

hapmap <- hapmap[hapmap$Chrm == qtl$Chrm,]
M <- M[, hapmap$SNP]
dim(M); dim(hapmap)

map1 <- data.frame(snp.id = as.integer(1:dim(hapmap)[1]), 
                   snp.rs.id = hapmap$SNP, 
                   chr = as.integer(hapmap$Chrm), 
                   pos = as.integer(hapmap$Pos))
head(map1)

# creating the GDS file
snpgdsCreateGeno(gds.fn = "./toy.gds",                              # gds filename
                 genmat = M,                                        # markers matrix
                 sample.id = rownames(M),                           # individual names
                 snp.id = map1$snp.id,                               # snp id (deve ser único)
                 snp.rs.id = map1$snp.rs.id,                         # snp name (pode não ser único)
                 snp.chromosome = map1$chr,                          # cromossomo
                 snp.position = map1$pos,                            # posição dentro do cromossomo (bp)
                 #snp.allele = map$allele,                           # alelos (ref / alt)
                 snpfirstdim = FALSE)                                # argumento para matriz n x s

# Loading gds
genofile <- snpgdsOpen(filename = "./toy.gds")
genofile

# estimating the pairwised LD, considering the M used in the GWAS
LDs <- abs(as.matrix(snpgdsLDMat(genofile, method = "corr", slide = 0, num.thread = 8)$LD))
dim(LDs)
rownames(LDs) <- colnames(LDs) <- hapmap$SNP
LDs[1:5, 1:5]
# close GDS
snpgdsClose(genofile)

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
distance <- naive_dist(as.matrix(hapmap$Pos), as.matrix(hapmap$Pos)) 
rownames(distance) <- colnames(distance) <- hapmap$SNP
distance2 <-  melt(distance)
head(distance2)
dim(distance2)
pairLD <- data.frame(marker1 = lds$Var1, marker2 = lds$Var2, r2 = lds$value, dist = distance2$value)
head(pairLD)
pairLD <- pairLD[pairLD$dist != 0,]
dim(pairLD)

# estimate the average LD decay
(averageLD <- round(mean(na.omit(pairLD$r2)), 2))
(decay <- mean(na.omit(pairLD$dist[round(pairLD$r2, 2) == averageLD])) /1E6)

# plot SNP neighbors 
library(gaston)
LD.plot(LDs[c((which(hapmap$SNP == qtl$SNP)-5):(which(hapmap$SNP == qtl$SNP)+5)),c((which(hapmap$SNP == qtl$SNP)-5):(which(hapmap$SNP == qtl$SNP)+5))], 
        hapmap$Pos[c((which(hapmap$SNP == qtl$SNP)-5):(which(hapmap$SNP == qtl$SNP)+5))],
        graphical.par = list(cex = 1.3, bg = "gray"), 
        polygon.par = list(border = NA))

################### the end ############