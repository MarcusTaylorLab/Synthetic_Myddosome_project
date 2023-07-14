pdb <- read.pdb("/Users/u_lichtenstein/hpc_data/sym/BDLD_57_x10_model_2_model_AG.pdb")

int <- atom.select(pdb, chain = LETTERS[seq(from = 1, to = 24, by = 2)]) 
ind <- atom.select(pdb, chain = LETTERS[seq(from = 2, to = 24, by = 2)]) 

# trim PDB to selection
pdb2 <- trim.pdb(pdb, int)
pdb3 <- trim.pdb(pdb, ind)

# Define the old chain identifiers
old_chains <- rev(LETTERS[seq(from = 1, to = 24, by = 2)])

# Define the new chain identifiers
new_chains <- letters[seq(from = 1, to = 12, by = 1)]
new_chains2 <- LETTERS[seq(from = 1, to = 12, by = 1)]

# Iterate through each chain
for (i in seq_along(old_chains)) {
  # Rename the chain
  pdb2$atom$chain[pdb2$atom$chain == old_chains[i]] <- new_chains[i]
  print(old_chains[i])
}

for (i in seq_along(new_chains)) {
  # Rename the chain
  pdb2$atom$chain[pdb2$atom$chain == new_chains[i]] <- new_chains2[i]
  print(new_chains[i])
}

# Define the old chain identifiers
old_chains <- LETTERS[seq(from = 2, to = 24, by = 2)]

# Define the new chain identifiers
new_chains <- letters[seq(from = 13, to = 24, by = 1)]  
new_chains2 <- LETTERS[seq(from = 13, to = 24, by = 1)]  

# Iterate through each chain
for (i in seq_along(old_chains)) {
  # Rename the chain
  pdb3$atom$chain[pdb3$atom$chain == old_chains[i]] <- new_chains[i]
}

for (i in seq_along(new_chains)) {
  # Rename the chain
  pdb3$atom$chain[pdb3$atom$chain == new_chains[i]] <- new_chains2[i]
  print(new_chains[i])
}

renamed_pdb <- cat.pdb(pdb2, pdb3, rechain = FALSE)

setwd("/Users/u_lichtenstein/hpc_data/sym/")

# write the new pdb object to file
write.pdb(renamed_pdb, file="BDLD_57_x10_model_2_model_AG_selected_chains.pdb")

