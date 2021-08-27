# This script runs PGLS on connectivity asymmetry and cerebral volume.

# Load libraries
library("ape")
library("geiger")
library("caper")
source("tidyPgls.R")

# Read consensus tree from 10ktrees in NEXUS format
tree <- read.nexus("../data/consensusTree_10kTrees_Primates_Version3.nex")

# Normalized data (for easier interpretation of coefficients)

# Read data into a dataframe
mydata <- read.csv("../pgls/connectivityAsymmetryNormalized.csv")

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
df1 <- tidyPgls(pgls(conn_asymmetry~SupraTentorialVol, data=comp.data, lambda=mylambda)) # connectivity asymmetry vs supratentorial volume
df2 <- tidyPgls(pgls(conn_asymmetry~cort_surf, data=comp.data, lambda=mylambda)) # connectivity asymmetry vs cortical surface area
results = rbind(df1, df2)

write.csv(results, "../pgls/pglsConnectivityAsymmetryNormalized.csv")

# Unnormalized data (for figures)

# Read data into a dataframe
mydata <- read.csv("../pgls/connectivityAsymmetry.csv")

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
df1 <- tidyPgls(pgls(conn_asymmetry~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # connectivity asymmetry vs supratentorial volume
df2 <- tidyPgls(pgls(conn_asymmetry~log(cort_surf), data=comp.data, lambda=mylambda)) # connectivity asymmetry vs cortical surface area
results = rbind(df1, df2)

write.csv(results, "../pgls/pglsConnectivityAsymmetry.csv")


