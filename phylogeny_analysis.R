# Phylogenetic analysis of pinniped COX1 genes
# Author: Mona Rosanally

# Load required packages
if(!require(msa)) install.packages("msa", dependencies=TRUE)
if(!require(phangorn)) install.packages("phangorn", dependencies=TRUE)
if(!require(Biostrings)) install.packages("Biostrings", dependencies=TRUE)

library(msa)
library(phangorn)
library(Biostrings)

# Set working directory (adjust as needed)
setwd("~/Documents/Phylogenetic_Analysis_Pinniped_COX1_Genes")

# ----------------------------
# 1. Load sequence data
# ----------------------------
Pinniped_init <- readDNAStringSet("Pinnipeds_initdata.txt")
Pinniped_new <- readDNAStringSet("Pinnipeds_newdata.txt")

# ----------------------------
# 2. Multiple Sequence Alignment
# ----------------------------
msa_init <- msa(Pinniped_init)
msa_new <- msa(Pinniped_new)

# Save MSA output to text
sink("msa_new_output.txt")
print(msa_new, show="complete")
sink()

# Convert MSA to phangorn format
phyDat_new <- as.phyDat(msa_new)

# ----------------------------
# 3. Neighbor-Joining Tree
# ----------------------------
dist_new <- dist.ml(phyDat_new, model="F81")
NJ_tree <- NJ(dist_new)

# Save high-resolution NJ tree
png("NJ_tree.png", width=2000, height=1500, res=300)
plot.phylo(NJ_tree, cex=1, main="Neighbor-Joining Tree")
dev.off()

pdf("NJ_tree.pdf", width=10, height=7)
plot.phylo(NJ_tree, cex=1, main="Neighbor-Joining Tree")
dev.off()

# ----------------------------
# 4. Model Testing for ML Tree
# ----------------------------
mt <- modelTest(phyDat_new, model=c("JC","F81","K80","HKY","SYM","GTR"))
best_model <- as.character(mt$Model[which.min(mt$AICc)])
best_model <- gsub("\\+.*","",best_model)  # remove +G+I if present

# Initialize ML model
fitStart <- eval(get(best_model, attr(mt,"env")), attr(mt,"env"))

# ----------------------------
# 5. Maximum Likelihood Tree
# ----------------------------
fit <- optim.pml(fitStart, rearrangement="stochastic", optGamma=TRUE, optInv=TRUE, model=best_model)

# Bootstrap analysis
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)

# Save high-resolution ML tree with bootstrap
png("ML_tree.png", width=2000, height=1500, res=300)
plotBS(midpoint(fit$tree), bs, p=50, type="p", cex=0.8, main="Maximum Likelihood Tree")
dev.off()

pdf("ML_tree.pdf", width=10, height=7)
plotBS(midpoint(fit$tree), bs, p=50, type="p", cex=0.8, main="Maximum Likelihood Tree")
dev.off()