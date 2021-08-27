# This script runs PGLS on network density and cerebral volume.

# Load libraries
library("ape")
library("geiger")
library("caper")
source("tidyPgls.R")

# Read consensus tree from 10ktrees in NEXUS format
tree <- read.nexus("../data/consensusTree_10kTrees_Primates_Version3.nex")

# Unnormalized data (for figures)
# Read data into a dataframe
mydata <- read.csv("../pgls/density.csv")

# Match species names in the data to species names in the phylogenetic tree
matching <- data.frame(tree$tip.label)
matching$subject <- c("lagothrix", "aotus", "cebus", "gorilla", "human", "gibbon", "bonobo", "chimpanzee", "orangutan", "pithecia", "lophocebus", "macaque", "colobus", "galago")
comp_data <- merge(mydata, matching)
col_idx <- match("tree.tip.label", names(comp_data))
comp_data <- comp_data[, c(col_idx, (1:ncol(comp_data))[-col_idx])]
colnames(comp_data)[colnames(comp_data) == "tree.tip.label"] <- "names"
missing_data <- name.check(tree, comp_data, comp_data$names) # check that this returns OK

# Create comparative data frame
comp.data <- comparative.data(tree, comp_data, names.col = "names", vcv.dim=2, warn.dropped=TRUE)

# Analyses with estimated lambda
mylambda <- "ML"
df1 <- tidyPgls(pgls(densityLH~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # LH density vs supratentorial volume
df2 <- tidyPgls(pgls(densityRH~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # RH density vs supratentorial volume
results = rbind(df1, df2)

write.csv(results, "../pgls/pglsDensity.csv")

# Normalized data (for regression coefficients)
# Read data into a dataframe
mydata <- read.csv("../pgls/densityNormalized.csv")

# Match species names in the data to species names in the phylogenetic tree
matching <- data.frame(tree$tip.label)
matching$subject <- c("lagothrix", "aotus", "cebus", "gorilla", "human", "gibbon", "bonobo", "chimpanzee", "orangutan", "pithecia", "lophocebus", "macaque", "colobus", "galago")
comp_data <- merge(mydata, matching)
col_idx <- match("tree.tip.label", names(comp_data))
comp_data <- comp_data[, c(col_idx, (1:ncol(comp_data))[-col_idx])]
colnames(comp_data)[colnames(comp_data) == "tree.tip.label"] <- "names"
missing_data <- name.check(tree, comp_data, comp_data$names) # check that this returns OK

# Create comparative data frame
comp.data <- comparative.data(tree, comp_data, names.col = "names", vcv.dim=2, warn.dropped=TRUE)

# Analyses with estimated lambda
mylambda <- "ML"
df1 <- tidyPgls(pgls(densityLH~SupraTentorialVol, data=comp.data, lambda=mylambda)) # LH density vs supratentorial volume
df2 <- tidyPgls(pgls(densityRH~SupraTentorialVol, data=comp.data, lambda=mylambda)) # RH density vs supratentorial volume
results = rbind(df1, df2)

write.csv(results, "../pgls/pglsDensityNormalized.csv")