#!/u/mlicht/miniconda3/bin/Rscript
library(bio3d)

#intercept variables from shell script
shell_var <- commandArgs(trailingOnly = TRUE)
vars <- strsplit(shell_var, "=", fixed = TRUE)

HOME_DIRECTORY <- vars[[1]][2]
STRUCTURE_PATH <- vars[[2]][2]

# Define input directory for the structure
INPUT <- paste(HOME_DIRECTORY, STRUCTURE_PATH, sep = "/")

pdb <- read.pdb(INPUT)

#Determine interacting chains
chain_A <- atom.select(pdb, chain = "A") # Select chain A

not_chain_A <- atom.select(pdb, chain = "A", inverse = TRUE) #select everything else

# Get the interacting chains
interacting_chains <- binding.site(a = pdb, b.inds = chain_A, a.inds = not_chain_A, cutoff = 5)

chains <- unique(as.vector(interacting_chains$chain))

cat(chains, sep = " ")
