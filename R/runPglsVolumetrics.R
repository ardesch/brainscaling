# This script runs PGLS on volumetrics such as white matter volume, cortical gray matter volume, cerebral volume etc.

# Load libraries
library("ape")
library("geiger")
library("caper")
source("tidyPgls.R")

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

# Create comparative data frame
comp.data <- comparative.data(tree, comp_data, names.col = "names", vcv.dim=2, warn.dropped=TRUE)

# Analyses with estimated lambda
mylambda <- "ML"
df1 <- tidyPgls(pgls(log(cort_surf)~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # cortical surface area vs cerebral volume
df2 <- tidyPgls(pgls(log(cort_surf)~log(CortexVol), data=comp.data, lambda=mylambda)) # cortical surface area vs cortical gray matter volume
df3 <- tidyPgls(pgls(log(cort_thickness)~log(CortexVol), data=comp.data, lambda=mylambda)) # cortical thickness vs cortical gray matter volume
df4 <- tidyPgls(pgls(log(CerebralWhiteMatterVol)~log(CortexVol), data=comp.data, lambda=mylambda)) # white matter vs cortical gray matter volume
df5 <- tidyPgls(pgls(log(cort_surf)~log(CerebralWhiteMatterVol), data=comp.data, lambda=mylambda)) # cortical surface area vs white matter
df6 <- tidyPgls(pgls(log(total_cc_area)~log(SupraTentorialVol), data=comp.data, lambda=mylambda)) # corpus callosum cross-sectional area vs cerebral volume
df7 <- tidyPgls(pgls(log(total_cc_area)~log(cort_surf), data=comp.data, lambda=mylambda)) # corpus callosum cross-sectional area vs cortical surface area
df8 <- tidyPgls(pgls(log(total_cc_area)~log(CerebralWhiteMatterVol), data=comp.data, lambda=mylambda)) # corpus callosum cross-sectional area vs white matter volume
df9 <- tidyPgls(pgls(log(total_cc_area)~log(CortexVol), data=comp.data, lambda=mylambda)) # corpus callosum cross-sectional area vs cortical gray matter volume

results <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)
write.csv(results, "../pgls/pglsVolumetrics.csv")
