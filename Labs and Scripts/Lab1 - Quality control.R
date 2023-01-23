##########################
# AGRO 7075	- Prediction-based Breeding
# Lab 1 - Quality control of marker and phenotypic data
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Jan 13 2023
##########################

############################################################
######################### Phenotypes #######################
############################################################

# loading the phenotyping file
pheno <- read.table("pheno.txt", header = T, na.strings = "NA")
head(pheno) # show the first six rows
tail(pheno) # show the last six rows
str(pheno) # present the flies' structure

# adjusting to factors
pheno[,1:9] <- lapply(pheno[,1:9], factor)
str(pheno)

# X and Y coordinate as a variate
pheno$X <- as.numeric(pheno$row)    
pheno$Y <- as.numeric(pheno$col)
str(pheno)

#saving the newest file
saveRDS(pheno, "pheno") # this kind of file is much better for R
# loading again...
pheno <- readRDS("pheno")
head(pheno)

# phenotypic correlation between traits
round(cor(pheno[ ,10:12], use = "pairwise"),2)
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
chart.Correlation(as.matrix(na.omit(pheno[,10:12])), histogram = TRUE, pch = 1)

# identifying outliers
boxplot(pheno$SDM, col = "red")
#install.packages("lme4")
library(lme4)
# outlier detection and elimination
fit <- lm(SDM ~ type + row + col + N + gid, data = pheno)
library(car)
(outlier <- names(outlierTest(fit)$p))
pheno[outlier, "SDM"] <- NA

# check the experimental design and spatial distribution
library(desplot)
d1 <- desplot(pheno, SDM ~ X*Y, out1 = rep, out2 = N,
              out2.gpar=list(col = "green", lwd = 1, lty = 1))
print(d1)

# testing for normality
# First lets check using patterns
shapiro.test(rnorm(length(pheno$SDM))) # normal distribution
shapiro.test(runif(length(pheno$SDM))) # uniform distribution
# then, 
shapiro.test(pheno$SDM)

#install.packages("bestNormalize")
require(bestNormalize)
SDMadj <- bestNormalize(pheno$SDM, standardize = FALSE, allow_orderNorm = TRUE, out_of_sample = FALSE)
SDMadj$chosen_transform
shapiro.test(SDMadj$x.t)
pheno$SDMadj <- SDMadj$x.t
head(pheno)

# What about the residuals?
# Quartile‐Quartile (Q‐Q) normality plot for residuals
fit <- lm(SDM ~ type + row + col + N + gid, data = pheno)
fit2 <- lm(SDMadj ~ type + row + col + N + gid, data = pheno)

par(mfrow = c(2,2)) # organize the plot window in 1 row and 2 col
qqnorm(resid(fit))
qqline(resid(fit), col = "red")
qqnorm(resid(fit2))
qqline(resid(fit2), col = "blue")
hist(pheno$SDM, col = "red", main = "SDM", xlab = "SDM")
hist(pheno$SDMadj, col = "blue", main = "Adjusted SDM", xlab = "Adjusted SDM")
dev.off()

# saving the newest pheno file
str(pheno)
head(pheno)
saveRDS(pheno, "pheno")

#####################################
############# Markers ###############
#####################################

#loading file
geno <- readRDS("geno")
head(geno)
dim(geno)

# let's check the SNPs and genotyping errors
(reads <- unique(unlist(apply(geno[,5:19], 2, unique))))

# create the hapmap file
hapmap <- geno[,c(1, 3, 4)]
head(hapmap)
tail(hapmap)
dim(hapmap)
str(hapmap)
hapmap$chrom <- as.factor(hapmap$chrom)
hapmap$pos <- as.numeric(hapmap$pos)
str(hapmap)

# marker file
M <- t(geno[,6:19])
M[1:5, 1:5]
dim(M)

# check if they are equal
all(colnames(M) == rownames(hapmap)) 

##################################### Quality control ######################

#source("https://bioconductor.org/biocLite.R")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("impute")
#devtools::install_github(repo = 'italo-granato/snpReady', ref = 'dev')
library(snpReady)
args(raw.data)

QC <- raw.data(data = as.matrix(M), 
               frame = "wide", 
               hapmap = hapmap, 
               sweep.sample = 1, 
               maf = 0.15, 
               call.rate = 0.90, 
               base = TRUE, 
               imput = TRUE, 
               imput.type = "wright", 
               outfile = "012", 
               plot = F)

# The report of the quality control approaches
QC$report

# get the newest dataset of markers
M <- QC$M.clean
dim(M)
M[1:5, 1:5]

# and the hapmap
hapmap <- QC$Hapmap
dim(hapmap)
head(hapmap)

# Then, identify and remove remove markers in heterozigosity
non.het.markers <- apply(M, 2, function(x){all(x != 1)})
sum(non.het.markers)
M <- M[, non.het.markers]
dim(M)

# then, correcting the map object again
hapmap <- hapmap[non.het.markers, ] # removing SNPs pruned by LD
head(hapmap)
dim(hapmap)

# verifying if all the markers into M matrix are in the same order into map
identical(as.character(hapmap$rs), colnames(M))

##############################
# pruning markers based on LD
##############################
# SNPRelate package
# to install the package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")
library(SNPRelate)

# hapmap
map2 <- data.frame(snp.id = as.integer(1:dim(hapmap)[1]), snp.rs.id = hapmap$rs, chr = as.integer(hapmap$chrom), pos = as.integer(hapmap$pos))
head(map2)

# creating the GDS file
snpgdsCreateGeno(gds.fn = "./toy.gds",                              # gds filename
                 genmat = M,                                        # markers matrix
                 sample.id = rownames(M),                           # individual names
                 snp.id = map2$snp.id,                               # snp id 
                 snp.rs.id = map2$snp.rs.id,                         # snp name 
                 snp.chromosome = map2$chr,                          # cromossomo
                 snp.position = map2$pos,                            # positionon chrom (bp)
                 #snp.allele = map$allele,                           # alleles (ref / alt)
                 snpfirstdim = FALSE)                                # argument for matrix n x s


# Loading gds
genofile <- snpgdsOpen(filename = "./toy.gds")
genofile

# prune markers by MAF 0.10, CR 0.90 e LD 0.99 (r2)
snps_pruned <- snpgdsLDpruning(gdsobj = genofile,         # gds file
                               remove.monosnp = TRUE,     # only bialellic
                               maf = 0.15,                # MAF
                               missing.rate = 0.90,       # CR
                               method = "corr",           # r2 (method of LD)
                               slide.max.bp = 100000,     # slide window 
                               ld.threshold = 0.98,       # LD threshold
                               start.pos = "first")       # start at

# get SNP ids
snps_pruned <- unlist(unname(snps_pruned))
snps_pruned
length(snps_pruned) # good quality snps 

# removing SNPs prunned by LD from M matrix
M <- M[, as.numeric(snps_pruned)]
M[1:5, 1:5]
dim(M)

# correcting the map object
head(hapmap)
hapmap <- hapmap[as.numeric(snps_pruned), ] # removing SNPs pruned by LD
head(hapmap)
dim(hapmap)

# verifying if all the markers into M matrix are in the same order into map
identical(as.character(hapmap$rs), colnames(M))

# close GDS
snpgdsClose(genofile)


######################## creating the marker set to hybrids as well ####################
phenoSC <- pheno[pheno$type == "sc",]
dim(phenoSC)
head(phenoSC)
sc.grid <- expand.grid(sort(unique(phenoSC$female)), sort(unique(phenoSC$male)))

M.female <- M[match(sc.grid[, 1], rownames(M)),]
M.male <- M[match(sc.grid[, 2], rownames(M)),]
dim(M.female); dim(M.male) ; dim(M)
M.SC <- (M.female + M.male) / 2
dim(M.SC)
rownames(M.SC) <- apply(sc.grid, 1, function(x){paste0(x[1], x[2])})
M.SC[1:5, 1:5]
M[1:5, 1:5]

# Finally combinane both M together
M <- rbind(M, M.SC)
dim(M)
all(unique(phenoSC$gid) %in% rownames(M))

# saving files
saveRDS(hapmap, "hapmap")
saveRDS(M, "M")

################# the end #############