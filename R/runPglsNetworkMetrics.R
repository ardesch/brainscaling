# This script runs PGLS on characteristic path length/clustering coefficient
# and cerebral volume.

# Load libraries
library("ape")
library("geiger")
library("caper")
source("tidyPgls.R")

# Read consensus tree from 10ktrees in NEXUS format
tree <- read.nexus("../data/consensusTree_10kTrees_Primates_Version3.nex")

## Normalized, lh
# Read data into a dataframe
mydata <- read.csv("../pgls/networkMetricsNormalized_lh.csv")

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
df1 <- tidyPgls(pgls(Cnorm~SupraTentorialVol, data=comp.data, lambda=mylambda)) # C lh
df2 <- tidyPgls(pgls(Lnorm~SupraTentorialVol, data=comp.data, lambda=mylambda)) # L lh

## Normalized, rh
# Read data into a dataframe
mydata <- read.csv("../pgls/networkMetricsNormalized_rh.csv")

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
df3 <- tidyPgls(pgls(Cnorm~(SupraTentorialVol), data=comp.data, lambda=mylambda)) # C rh
df4 <- tidyPgls(pgls(Lnorm~(SupraTentorialVol), data=comp.data, lambda=mylambda)) # L rh

results = rbind(df1, df2, df3, df4)

write.csv(results, "../pgls/pglsNetworkMetricsNormalized.csv")

## Unnormalized, lh
# Read data into a dataframe
mydata <- read.csv("../pgls/networkMetrics_lh.csv")

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
df1 <- tidyPgls(pgls(Cnorm~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # C lh
df2 <- tidyPgls(pgls(Lnorm~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # L lh

## Normalized, rh
# Read data into a dataframe
mydata <- read.csv("../pgls/networkMetrics_rh.csv")

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
df3 <- tidyPgls(pgls(Cnorm~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # C rh
df4 <- tidyPgls(pgls(Lnorm~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # L rh

results = rbind(df1, df2, df3, df4)

write.csv(results, "../pgls/pglsNetworkMetrics.csv")