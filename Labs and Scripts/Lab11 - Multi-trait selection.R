########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 11 - Multi-trait selection
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: April 12 2023
#######################################
# loading packages
require(foreach)
require(doParallel)
require(doMC)
library(car)
library(ggplot2)
library(sommer)

# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1) # type the number of cores you want to use
getDoParWorkers()

#loading phenotypes
pheno <- readRDS("pheno")
head(pheno)
dim(pheno)
str(pheno)

# reorganizing the file
pheno <- pheno[,c(1:9, 13:14, 10:12)]
head(pheno)

# index for traits
traits <- colnames(pheno)[12:14]
  
# estimate de phenotypic correlation
pheno.cor <-round(cor(pheno[,traits], use = "pairwise.complete.obs", method = "pearson"), 2)
corrplot::corrplot(pheno.cor, method = 'number', type = "lower", diag = F, col = c("red", "orange", "green", "blue"))

# reorganize de data
pheno.melted <- reshape2::melt(pheno, measure.vars = traits)
head(pheno.melted)  

# load the GRM
Ga <- readRDS("Ga")

# running all GS single-traits in parallel 
results.st <- foreach(i = 1:length(traits), 
                        .packages = c("sommer", "car"), 
                        .combine = "rbind",
                        .export = c("mmer", "predict.mmer", "outlierTest"),
                        .multicombine = TRUE, 
                        .errorhandling = "remove",
                        .verbose = TRUE    
  ) %dopar% {

    # subset the data  
    sample <- droplevels.data.frame(pheno.melted[pheno.melted$variable == traits[i],])
    
    # outlier detection and elimination
    fit <- lm(value ~ rep + gid + N + gid:N, data = sample)
    outlier <- names(outlierTest(fit)$p)
    sample[outlier, "value"] <- NA

    # MME using only the classical experimental design, and gid as random
    fitR <- mmer(fixed = value ~ N + rep,
                random = ~ vsr(gid, Gu = Ga),
                rcov = ~ vsr(units),
                data = sample)

    # predicting the BLUPs - main effect
    BLUPs <- predict(object = fitR, D = "gid")$pvals

    output <- data.frame(BLUPs,
                        trait = unique(sample$variable))

}

head(results.st)
tail(results.st)

# GEBVs per trait
GEBVs <- tapply(results.st$predicted.value, list(results.st$gid, results.st$trait), mean)
head(GEBVs)

######################## single step (ss) multi-trait (MT) GS ################

# MME using only the classical experimental design, and gid as random and no GxE
fitMT <- mmer(cbind(SDM,   SRA,   NAE) ~ N + rep,
              random = ~vsr(gid, Gu = Ga,  Gtc = unsm(3)),
              rcov = ~ vsr(units), 
              tolParInv = 1,
              data = pheno)

# G varcomp via MT-GBLUP
fitMT$sigma$`u:gid`

# correlation among traits via MT-GBLUP
cov2cor(fitMT$sigma$`u:gid`)

# retrieving BLUPs - main effect
BLUPs.MT <- data.frame(fitMT$U$`u:gid`)

# look what happens with the correlations between GEBVs
cor(BLUPs.MT)

# So, the only thing is just weight the GEBV and obtain the SI
# define the economic weights per trait
ecoW <- c(1, 1, 1) # in this case two times more for grain yield
# then, the vector of SI per genotype
MT <- as.matrix(BLUPs.MT) %*% ecoW  
head(MT)


######################################### Selection indices ##################################
# phenotypic covariance between traits
P <- cov(pheno[, traits], use = "pairwise", method = "pearson")

# genetic covariance between traits
G <- cov(GEBVs)
cor(GEBVs)

# Smith-Hazel
# define the economic weights per trait
ecoW <- c(1, 1, 1) # in this case two times more for grain yield
# then, the selection weights per trait
(b.sh <- solve(as.matrix(P)) %*% as.matrix(G) %*% as.matrix(ecoW))
# Finally, the vector of SI per genotype
SH <- GEBVs %*% b.sh  
head(SH)

# Pasek-Baker
# define the desired genetic gains per traits in genetic standard deviations
desired <- c(1, 1, 1)
# G correlation matrix x desired genetic gains in standard deviations
(b.pb <- solve(cov(scale(GEBVs))) %*% as.matrix(desired))
PB <- as.matrix(scale(GEBVs)) %*% b.pb
head(PB)

# correlation between GEBVs and selection indices
GEBVs <- as.data.frame(GEBVs)
GEBVs$SH <- SH
GEBVs$PB <- PB
GEBVs$MT <- MT
cor(GEBVs)

# let's see which method provides the best KPIs
apply(cor(GEBVs[, 4:6]), 2, mean)
apply(cor(GEBVs[, 4:6]), 2, sd)

# defining a threshold to select the individuals to advance
is <- 15 #15%
quantile(GEBVs$SH, ((100 - is) / 100))
GEBVs$Seleted <- GEBVs$SH > quantile(GEBVs$SH, ((100 - is) / 100))
# the number of pre-selected parents
sum(GEBVs$Seleted)
head(GEBVs)
# saving the file
write.csv(GEBVs, "GEBVS.csv")

# 3D plot showing the selected materials and parents
library(ggplot2)
library(ggpubr)
a <- ggplot(GEBVs, aes(x = SDM, 
                       y = SRA,
                       color = Seleted)) + geom_point()

b <- ggplot(GEBVs, aes(x = SDM, 
                       y = NAE,
                       color = Seleted)) + geom_point()

c <- ggplot(GEBVs, aes(x = NAE, 
                       y = SRA,
                       color = Seleted)) + geom_point()

g <-
  ggarrange(a, b, c,
            ncol = 1,
            nrow = 3,
            common.legend = F,
            labels = c('a', 'b', 'c'),
            font.label = c(size = 8),# label.x = c('cycle'),
            legend = 'bottom')

# creating a bottom label
g <- annotate_figure(g, bottom = text_grob(label = 'Traits',face = 'bold', size = 6))

# saving graph
ggsave(filename = "selected_lines.pdf",
       plot = g,
       device = 'pdf',
       width = 250,
       height = 400,
       units = 'mm',
       dpi = 300)

#################### the end ####################################