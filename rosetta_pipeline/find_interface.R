library(pacman)
p_load(bio3d)

INPUT_PATH <- paste()

pdb <- read.pdb("/Users/u_lichtenstein/hpc_data/from_rcsb/6I3N.pdb")

#Determine interacing chains
chain_A <- atom.select(pdb, chain = "A") # Select chain A

not_chain_A <- atom.select(pdb, chain = "A", inverse = TRUE) #select everything else

# Get the interacting chains
interacting_chains <- binding.site(a = pdb, b.inds = chain_A, a.inds = not_chain_A, cutoff = 5)

chains <- unique(as.vector(interacting_chains$chain))

print(chains)
