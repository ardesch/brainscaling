# This script runs alternative phylogenetic regression models,
# including OU, EB, and RMA.

# Load libraries
library("nlme")
library("ape")
library("geiger")
library("phytools")

# Read consensus tree from 10ktrees in NEXUS format
tree <- read.nexus("../data/consensusTree_10kTrees_Primates_Version3.nex")

# Read data into a dataframe
mydata <- read.csv("../pgls/volumetricData.csv")

# Match species names in the data to species names in the phylogenetic tree
matching <- data.frame(tree$tip.label)
matching$subject <- c("lagothrix", "aotus", "cebus", "gorilla", "human", "gibbon", "bonobo", "chimpanzee", "orangutan", "pithecia", "lophocebus", "macaque", "colobus", "galago")
comp_data <- merge(mydata, matching)
col_idx <- match("tree.tip.label", names(comp_data))
comp_data <- comp_data[, c(col_idx, (1:ncol(comp_data))[-col_idx])]
colnames(comp_data)[colnames(comp_data) == "tree.tip.label"] <- "names"
missing_data <- name.check(tree, comp_data, comp_data$names) # check that this returns OK

## RMA (phytools) ##

# Cerebral volume ~ cortical surface area
x <- setNames(log(comp_data$SupraTentorialVol), comp_data$names)
y <- setNames(log(comp_data$cort_surf), comp_data$names)
fit.ml <- phyl.RMA(x, y, tree, method="lambda", h0=2/3)

# Gray matter volume ~ white matter volume
x <- setNames(log(comp_data$CortexVol), comp_data$names)
y <- setNames(log(comp_data$CerebralWhiteMatterVol), comp_data$names)
fit.ml <- phyl.RMA(x, y, tree, method="lambda", h0=1)

# White matter volume ~ cortical surface area
x <- setNames(log(comp_data$CerebralWhiteMatterVol), comp_data$names)
y <- setNames(log(comp_data$cort_surf), comp_data$names)
fit.ml <- phyl.RMA(x, y, tree, method="lambda", h0=2/3)

# Cortical surface area ~ cc cross-sectional area
x <- setNames(log(comp_data$cort_surf), comp_data$names)
y <- setNames(log(comp_data$total_cc_area), comp_data$names)
fit.ml <- phyl.RMA(x, y, tree, method="lambda", h0=1)

## Brownian motion, Ornstein-Uhlenbeck, and early-burst models ##

mydata <- read.csv("../pgls/volumetricData.csv")
matching <- data.frame(tree$tip.label)
matching$subject <- c("lagothrix", "aotus", "cebus", "gorilla", "human", "gibbon", "bonobo", "chimpanzee", "orangutan", "pithecia", "lophocebus", "macaque", "colobus", "galago")
comp_data <- merge(mydata, matching)
col_idx <- match("tree.tip.label", names(comp_data))
comp_data <- comp_data[, c(col_idx, (1:ncol(comp_data))[-col_idx])]
colnames(comp_data)[colnames(comp_data) == "tree.tip.label"] <- "names"
missing_data <- name.check(tree, comp_data, comp_data$names) # check that this returns OK
row.names(comp_data) <- comp_data$names
comp_data = comp_data[, -2]
comp_data = comp_data[, -1]
comp_data = log(comp_data) # log transform the whole data set

## Compare trait AICs ##
bmFit <- fitContinuous(phy=tree, dat=comp_data, model="BM") 
ouFit <- fitContinuous(phy=tree, dat=comp_data, model="OU") 
ebFit <- fitContinuous(phy=tree, dat=comp_data, model="EB") 
