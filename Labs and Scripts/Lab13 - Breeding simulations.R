########################################
# AGRO 7075	- Prediction-based Breeding
# Lab 13 - Creating and comparing breeding schemes
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: May 8 2023
#######################################
#install.packages("R6")
#install.packages("AlphaSimR")
(start.time <- Sys.time())
require(AlphaSimR)

#################### Parameters ##############
nQTL <- 360
add <- 0.22  
segSites <- 1644
chip.size.1 <- 540
replicates <- 18 # usually the number must be near to 100
cycles <- 20 # it is important define your breeding horizon

# Heritabilities by breeding stage
#        F2,   F3,   F4,   F5,   F6,   F7
h2 <- c(0.03, 0.15, 0.40, 0.60, 0.70, 0.80)
H2 <- c(0.06, 0.20, 0.45, 0.63, 0.72, 0.81)

# crossing block
n.parents <- 40
n.crosses <- 160
progenie.size <- 100

###########################
# crop history of evolution
###########################
set.seed(29121983)

history <- runMacs(nInd = 1000, 
                   nChr = 12,
                   segSites = segSites / 12, 
                   inbred = TRUE,
                   species = "GENERIC", 
                   split = NULL, 
                   ploidy = 2L)

# then, add the simulation parameters to the simulated population
SP = SimParam$new(history)

# Lets add a trait, considering gamma distribution for QTL effect;  A + D effects
SP$restrSegSites(nQTL/12, chip.size.1/12)
SP$addTraitAD(nQtlPerChr = nQTL / 12, gamma = TRUE, mean = 0, var = 1, meanDD = add, varDD = 0.5)

# add a SNP chip to our population
SP$addSnpChip(chip.size.1/12)

# after the historical population (mimic the crop evolution, it is time to simulate the first years of a breeding program, like as a warm up process)
# setting the heritability of the trait, narrow and broad sense
SP$setVarE(h2 = h2[4], H2 = H2[4])
firstPop <- newPop(rawPop = history, isDH = TRUE, simParam = SP)

#########################################################
# simulating the founders
#########################################################
# select the first parents, full inbreed
set.seed(29121983)


founders <- selectInd(firstPop, nInd = n.parents, trait = 1, use = "pheno", gender = "B",
                 selectTop = TRUE, returnPop = TRUE, candidates = NULL,
                 simParam = SP)


#############################################################
# simulating 3 generations of rice trad breeding == burn-in
#############################################################
cat("Burn-in", "\n")

# intensity of selection
nSelF3 <- 11200
nSelF4 <- 1280
nSelF5 <- 40
nSelF6 <- 24
nSelF7 <- 1

pop.trad <- founders
perform <- data.frame()

for (i in 1:3) {
  cat("Processing the cycle", i, "\n")

  f1 <- randCross(pop = pop.trad, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
  
  # setting the heritability of the trait, narrow and broad sense
  SP$setVarE(h2 = h2[1], H2 = H2[1])
  f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
  f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # setting the heritability of the trait, narrow and broad sense
  SP$setVarE(h2 = h2[2], H2 = H2[2])
  f3 <- self(pop = f2, nProgeny = 1, simParam = SP)
  f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # setting the heritability of the trait, narrow and broad sense
  SP$setVarE(h2 = h2[3], H2 = H2[3])
  f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
  f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # setting the heritability of the trait, narrow and broad sense
  SP$setVarE(h2 = h2[4], H2 = H2[4])
  f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
  f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # setting the heritability of the trait, narrow and broad sense
  SP$setVarE(h2 = h2[5], H2 = H2[5])
  f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
  f6sel <- selectInd(pop = f6, nInd = nSelF7, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
  
  # setting the heritability of the trait, narrow and broad sense
  SP$setVarE(h2 = h2[6], H2 = H2[6])
  variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
  
  # then, select the best one regardless the family
  newparents <- f5
  pop.trad <- newparents
  perform <- rbind(perform, data.frame(population = as.numeric(meanG(pop.trad)), cultivar = as.numeric(max(variety@gv))))
}

perform

#########################
# creating the first and small TS to be updated - a population of 1152 - 3 chips
#########################
cat("TS", "\n")
#cross everyone
SP$setVarE(h2 = h2[4], H2 = H2[4])
TS0 <- randCross(pop = pop.trad, nCrosses = 30, nProgeny = 40, simParam = SP)[1:1152]

# recombining the offspring
TS00 <- newPop(rawPop = TS0, simParam = SP)

# obtaining full inbred lines by DH tech
TS000 <- makeDH(pop = TS00, nDH = 1, useFemale = TRUE, simParam = SP)
markers <- RRBLUP(pop = TS000, traits = 1, use = "pheno", snpChip = 1, simParam = SP)

#Evaluate accuracy
TS000 = setEBV(TS000, markers, simParam=SP)
(Ac.GS <- as.numeric(cor(gv(TS000), ebv(TS000))))
capture.output(c(Ac.GS = Ac.GS), file = "GS.accuracies.txt")

# storage the first results
# results C0
resultsC0 <-  data.frame(
  method = c("Current_Trad", "Drift", "GS.F2_HTP.F3"),
  rep = 1,
  cycle = 0,
  PM = as.numeric(meanG(pop.trad)),
  Va = as.numeric(varA(pop.trad)),
  Variety = perform[3,2],
  Years = c(5, 5, 3)
)

#####################################################
# simulating more cycles of breeding by the current method
#####################################################
cat("Current_Trad", "\n")

# intensity of selection
nSelF3 <- 11200
nSelF4 <- 1600
nSelF5 <- 160
nSelF6 <- 40
nSelF7 <- 1

results.current.trad <- data.frame()

for(k in 1:replicates){
  
  pop.base <- pop.trad
  
  for (i in 1:cycles) {
    
    cat("running cycle", i, "in replicate", k, "\n")

    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[2], H2 = H2[2])
    f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
    f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f6, nInd = nSelF7, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
       
    newparents <- f6
    pop.base <- newparents

    results.current.trad <- rbind(results.current.trad, data.frame(
      method = "Current_Trad",
      rep = k,
      cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Variety = as.numeric(gv(variety)),
      Years = 5
      ))
  }
}


########################################################################
# simulating more cycles of breeding by the current method but by DRIFT
#########################################################################
cat("Drift", "\n")

results.current.drift <- data.frame()

for(k in 1:replicates){
  
  pop.base <- pop.trad
  
  for (i in 1:cycles) {
    
    cat("running cycle", i, "in replicate", k, "\n")
    
    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nSelF3, trait = 1, "rand", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[2], H2 = H2[2])
    f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
    f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "rand", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f4 <- self(pop = f3sel, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "rand", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f5sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "rand", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f6 <- self(pop = f5sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f6, nInd = nSelF7, trait = 1, "rand", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
    
    newparents <- f6
    pop.base <- newparents
    
    results.current.drift <- rbind(results.current.drift, data.frame(
      method = "Drift",
      rep = k,
      cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Variety = as.numeric(gv(variety)),
      Years = 5
    ))
  }
}


###########################################################################
# simulating more cycles of breeding by the GS in F2 and HTP in F3 rows
###########################################################################
cat("GS.F2_HTP.F3", "\n")

# intensity of selection
nTGF2 <- 11200
nSelF3 <- 960
nSelF4 <- 480
nSelF5 <- 40
nSelF6 <- 1

results.gshtp <- data.frame()

for(k in 1:replicates){
  
  pop.base <- pop.trad
  TS <- TS000
  popList <- list(TS)
  snps <- markers
  
  for (i in 1:cycles) {
    
    cat("running cycle", i, "in replicate", k, "\n")
    
    f1 <- randCross(pop = pop.base, nCrosses = n.crosses, nProgeny = 1, simParam = SP)
    
    SP$setVarE(h2 = h2[1], H2 = H2[1])
    f2 <- self(pop = f1, nProgeny = progenie.size, simParam = SP)
    f2sel <- selectInd(pop = f2, nInd = nTGF2, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # GS
    f2sel <- setEBV(pop = f2sel, solution = snps, value = "bv", simParam = SP)
    f2sel <- selectInd(pop = f2sel, nInd = nSelF3, trait = 1, "ebv", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense - in his case HTP
    SP$setVarE(h2 = h2[3], H2 = H2[3])
    f3 <- self(pop = f2sel, nProgeny = 1, simParam = SP)
    f3sel <- selectInd(pop = f3, nInd = nSelF4, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[4], H2 = H2[4])
    f4 <- self(pop = f3, nProgeny = 1, simParam = SP)
    f4sel <- selectInd(pop = f4, nInd = nSelF5, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[5], H2 = H2[5])
    f5 <- self(pop = f4sel, nProgeny = 1, simParam = SP)
    f6sel <- selectInd(pop = f5, nInd = nSelF6, trait = 1, "pheno", gender = "B", selectTop = TRUE, simParam = SP)
    
    # setting the heritability of the trait, narrow and broad sense
    SP$setVarE(h2 = h2[6], H2 = H2[6])
    variety <- self(pop = f6sel, nProgeny = 1, simParam = SP)
    
    # updating the TS and marker effects
    TSi <- f4
    popList <- c(popList, TSi)
    if(length(popList) <= 3) {
      TS <- mergePops(popList)
    }
    if(length(popList) > 3) {
      TS <- mergePops(popList[(length(popList)-2):length(popList)]) 
    }
    snps <- RRBLUP2(pop = TS, traits = 1, use = "pheno", snpChip = 1, simParam = SP, maxIter = 20)
    
    newparents <- f5
    pop.base <- newparents
    
    results.gshtp <- rbind(results.gshtp, data.frame(
      method = "GS.F2_HTP.F3",
      rep = k,
      cycle = i,
      PM = as.numeric(meanG(pop.base)),
      Va = as.numeric(varA(pop.base)),
      Variety = as.numeric(gv(variety)),
      Years = 3
    ))
  }
}


########### saving the outputs
cat("saving the output", "\n")
output <- rbind(resultsC0, results.current.trad, results.current.drift, results.gshtp)
saveRDS(output, paste("output", "LSU_Rice_BP", sep = ""))
######################################################
stopImplicitCluster()
(end.time <- Sys.time())
(start.time - end.time)


################################################ graphs #######################################

library("tidyr")
library('dplyr')
library("ggplot2")
library('ggpubr')

rm(list = ls()) # clear the environment

data <- readRDS("outputLSU_Rice_BP")
unique(data$method)
data$method <- as.factor(data$method)
str(data)

PM <-
  data %>%
  group_by(method, cycle) %>% 
  summarise(PM = mean(PM)) %>% 
  ggplot(data = ., mapping = aes(x = cycle, y = PM, col = method)) + 
  geom_line(size = 0.4) +
  geom_point(size = 0.7) +
  scale_x_continuous(limits = c(0, 20),
                     breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)) +
  labs(x = 'cycles',
       y = 'Population mean') +
  theme(axis.title = element_text(size = 5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 6, face = 'bold', margin = margin()),
        legend.position = 'bottom',
        legend.text = element_text(size = 5, hjust = 0),
        legend.title = element_text(size = 6, face = 'bold'),
        legend.key.height = unit(3, 'mm'))

Variety <-
  data %>%
  group_by(method, cycle) %>% 
  summarise(Variety = mean(Variety)) %>% 
  ggplot(data = ., mapping = aes(x = cycle, y = Variety, col = method)) + 
  geom_line(size = 0.4) +
  geom_point(size = 0.7) +
  scale_x_continuous(limits = c(0, 20),
                     breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)) +
  labs(x = 'cycles',
       y = 'Variety performance') +
  theme(axis.title = element_text(size = 5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 6, face = 'bold', margin = margin()),
        legend.position = 'bottom',
        legend.text = element_text(size = 5, hjust = 0),
        legend.title = element_text(size = 6, face = 'bold'),
        legend.key.height = unit(3, 'mm'))

Va <-
  data %>%
  group_by(method, cycle) %>% 
  summarise(Va = mean(Va)) %>% 
  ggplot(data = ., mapping = aes(x = cycle, y = Va, col = method)) + 
  geom_line(size = 0.4) +
  geom_point(size = 0.7) +
  scale_x_continuous(limits = c(0, 20),
                     breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)) +
  labs(x = 'cycles',
       y = 'Genetic variability') +
  theme(axis.title = element_text(size = 5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 6, face = 'bold', margin = margin()),
        legend.position = 'bottom',
        legend.text = element_text(size = 5, hjust = 0),
        legend.title = element_text(size = 6, face = 'bold'),
        legend.key.height = unit(3, 'mm'))


# grouping graphs
g <-
  ggarrange(PM, Variety, Va,
            ncol = 1,
            nrow = 3,
            common.legend = TRUE,
            labels = c('a', 'b', 'c'),
            font.label = c(size = 8),# label.x = c('cycle'),
            legend = 'right')

# creating a bottom label
g <- annotate_figure(g, bottom = text_grob(label = 'cycles',face = 'bold', size = 6))

g 

##### graphs for RS after 15 years of breeding

fifteen <-
  data %>%
  group_by(method, cycle)

scheme <- unique(fifteen$method)
n.cycles <- c(3, 3, 5)

final <- data.frame()
for(i in 1:length(scheme)){
  final <- rbind(final,
                 fifteen[fifteen$method == scheme[i] & fifteen$cycle == n.cycles[i],]) 
}

# remove the cycle zero intercept
final$RS <-  final$PM - mean(data[data$cycle == 0, "PM"])
final$RS_percyear <-  (final$PM - mean(data[data$cycle == 0, "PM"])) / mean(data[data$cycle == 0, "PM"]) * 100 / 15

# and looking at the distribution and RS
library(ggplot2)
p <- ggplot(final, aes(x = reorder(method, RS), y = RS)) + 
  geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  geom_jitter(position=position_jitter(0.2)) +
  labs(x = "Breeding method", y = "Response to selection after 15 years of breeding") +
labs(x = "", fill = "Breeding Method")

p

############ the end #############