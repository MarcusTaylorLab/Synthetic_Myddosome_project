#!/usr/local/bin/Rscript
library(bio3d)

#intercept variables from shell script
shell_var <- commandArgs(trailingOnly = TRUE)
vars <- strsplit(shell_var, "=", fixed = TRUE)

HOME_DIRECTORY <- vars[[1]][2]
STRUCTURE_PATH <- vars[[2]][2]

STRUCTURE_TITLE <- strsplit(
  strsplit(STRUCTURE_PATH, "/", fixed = TRUE)[[1]][2], ".", fixed = TRUE)[[1]][1]

# Define directories for the structure
INPUT <- paste(HOME_DIRECTORY, STRUCTURE_PATH, sep = "/")
OUTPUT <- paste(HOME_DIRECTORY, "output_files", STRUCTURE_TITLE, sep = "/")

# read pdb-file
pdb <- read.pdb(INPUT)

chain_interactions <- data.frame()

for (CHAIN in unique(as.vector(pdb$atom$chain))){
  
  #Determine interacting chains
  chain_select <- atom.select(pdb, chain = CHAIN) # Select chain
  not_chain_select <- atom.select(pdb, chain = CHAIN, inverse = TRUE) #select everything else
  
  # Get the interacting chains
  interacting_res <- binding.site(a = pdb, b.inds = chain_select, a.inds = not_chain_select, cutoff = 5)
  
  chain_interactions_temp <- data.frame(
    description = STRUCTURE_TITLE,
    chain = CHAIN,
    total_interactions = length(unique(as.vector(interacting_res$chain)))
  )
  
  chain_interactions_temp$chain_interactions = paste(unique(as.vector(interacting_res$chain)), collapse = " ")
  
  chain_interactions <- rbind(chain_interactions, chain_interactions_temp)
  
  rm(chain_interactions_temp, interacting_res, not_chain_select, chain_select)
}

min_chains <- chain_interactions %>% 
  filter(
    total_interactions == min(total_interactions)
  )

# Chain with minimal interaction partners
min_chains$chain[1]

# Interaction partners of selected chain
min_chains$chain_interactions[[1]]

setwd(OUTPUT)

#provide filename
filename <- paste("interfacing_chains_", STRUCTURE_TITLE, ".csv", sep = "")

write.csv(chain_interactions, file = filename)
