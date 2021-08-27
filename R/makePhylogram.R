# This script produces a phylogram as used in Figure 1.

# Load libraries
library("ggtree")
library("treeio")

# Read consensus tree from 10ktrees in NEXUS format
tree = read.nexus("../data/consensusTree_10kTrees_Primates_Version3.nex")

# Plot image
image <- ggtree(tree) + theme_tree2()
plot(image)