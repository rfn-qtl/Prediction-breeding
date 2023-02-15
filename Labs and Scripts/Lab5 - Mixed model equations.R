########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 5 - Mixed model equations
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Dec 1 2022
#######################################

#loading phenotypes
pheno <- readRDS("pheno")
head(pheno)
dim(pheno)
str(pheno)

#modeling
library(sommer)

# using only the classical experimental design
fit1 <- mmer(fixed = SDM ~ 1 + N,
              random = ~gid + rep + gid:N,
              rcov = ~ vsr(units),
              data = pheno)

summary(fit1)$varcomp
blocks <- length(unique(pheno$rep))
loc <- length(unique(pheno$N))

# Broad-sense heritability - full model
vpredict(fit1, h2 ~ V1 / (V1 + V3/loc + V4/(loc*blocks))) # trials level
vpredict(fit1, h2 ~ V1 / (V1 + V3 + V4)) # plot level

# H2 Cullis
Vg <- fit1$sigma$gid
ng <- length(unique(pheno$gid))
C22_g <- fit1$PevU$gid$SDM
trC22_g <-sum(diag(C22_g))
av2 <- 2/ng * (trC22_g - (sum(C22_g) - trC22_g) / ng-1) # mean var of a difference between genotypic BLUPS
(H2.Cullis <- 1 - abs(av2) / (2 * Vg))

# weights for ID's - adjust residual for further analysis
w <- diag(C22_g)

# predicting the BLUP - main effect
BLUPs <- predict.mmer(object = fit1, classify = "gid")
BLUPs$pvals

# reliability (r2) = 1 - PEV / (Vg + Vg*Fii)
# near to heritability
# If F == 0, then r = 1 - PEV/Vg
# If F == 1, then r = 1 - PEV/2Vg
(r <- mean(1 - BLUPs$pvals$standard.error^2 / as.numeric(Vg)))

# Accuracy = sqrt(r)
(Accuracy.blup <- mean(sqrt(r))) # BLUP
(Accuracy.pheno <- sqrt(vpredict(fit1, h2 ~ V1 / (V1+V3+V4))[1])) # pheno

# predicting the BLUP per enviroment
BLUPS.env <- predict.mmer(object = fit1, classify = c("gid","N"))
BLUPS.env$pvals


################################### confidence interval for BLUPS ######################
DMS <- BLUPs$pvals[, 4]*1.96
blups2 <- data.frame(BLUPs$pvals[,2:3], "DMS" = DMS)

library(ggplot2)
limits <- aes(ymax = blups2$predicted.value + blups2$DMS,
              ymin = blups2$predicted.value - blups2$DMS)
p <- ggplot(data = blups2, aes(x = reorder(factor(gid), -predicted.value), y = predicted.value))
p + geom_jitter(stat = "identity", colour = "red") +
  geom_errorbar(limits, position = position_dodge(0.5),
                width = 0.10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "gid", y = "BLUP")


################################ correlation between phenotypes and BLUPS #################################
means <- tapply(pheno$SDM, pheno$gid, mean)
means <- means[blups2$gid]
all(blups2$gid == names(means))
blups2$pheno <- means

plot(blups2$pheno, blups2$predicted.value, 
     main = "Phenotype x BLUP", 
     xlab = "Phenotypic mean", 
     ylab = "BLUP",
     col = "blue")
abline(lm(blups2$predicted.value ~ blups2$pheno), col = "red")
text(.7, 0.5, cor(blups2$pheno, blups2$predicted.value, use = "pairwise"), col = "black")

############################### selecting the best genotypes ###############################

# first, define a quantile, in this case 90%
(trshld <- quantile(blups2[, 2], probs = 0.90))

# Then identify the best ones
(selected <- blups2[blups2$predicted.value >= trshld,])

# and estimate the response to selection
# BLUPi = (Yi - Ypop) * hi
# RSi = DSi * hi
# hi = n*h / (1+(n-1)h) # n = number of observations of individual i
(RS <- mean(selected$predicted.value))


##################################### Comparing models - LRT, AIC, and BIC ###################

# to test model factor or compare model, we need to run other models eliminating / including factors 

# in the first example, let's include spatial information
fit2 <- mmer(fixed = SDM ~ 1 + N,
             random = ~gid + rep + gid:N +
             vsr(row) + vsr(col) + spl2Da(X,Y),
             rcov = ~ vsr(units),
             data = pheno)

plot(fit2)
summary(fit2)$varcomp

(lrt12 <- anova(fit1, fit2))

# in the next example, let's remove the GxE factor and test if it is significant
fit0 <- mmer(fixed = SDM ~ 1 + N,
             random = ~gid + rep,
             rcov = ~ vsr(units),
             data = pheno)

(lrt10 <- anova(fit1, fit0))

# how it works behind the scenes
2*(lrt10[1,4] - lrt10[2,4]) #LRT
pchisq(2*(lrt10[1,4] - lrt10[2,4]), 1, lower.tail = F) # chisq test

# wald test for fixed effects, in this case, N level
fit1$Beta$Estimate
wald.test(b = fit1$Beta$Estimate, Sigma = fit1$VarBeta, Terms = 2) # N

####### the end ##########