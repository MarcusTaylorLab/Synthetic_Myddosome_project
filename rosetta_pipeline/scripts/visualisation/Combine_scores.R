library(pacman)
p_load(tidyr, data.table, ggplot2, biomaRt)

filenames <- list.files("/Users/u_lichtenstein/hpc_data/merge_2", 
                        pattern="*.sc", 
                        full.names=TRUE)

score_list <- lapply(filenames, read.csv, sep = "", header = TRUE)

extend_fx <- function(x){
  
for (col in (paste("dG_", LETTERS[2:26], sep = ""))){
  
  #get the names of all interface energy columns. dG_AB and so on
  NAMES <- colnames(x)[-c(1:4)]
  
  if (!col %in%  NAMES) {
    #If the column doesn't exist add it with default value
    x[[col]] <- NA
  }
 }
  return(x)
}

score_all <- lapply(score_list, extend_fx)
score_all <- rbindlist(score_all, fill = TRUE)

# Get the number of observations in all structures and save as variable
observations <- colnames(score_all)[length(score_all)]

# Calculate number of interfaces per row
score_all$n_rows <- rowSums(!is.na(score_all[, 6:length(score_all)]))

# Remove columns with only NA values
#score_all <- score_all[, !apply(score_all, 2, function(x) all(is.na(x)))]

# Get actual names of proteins
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# get list of all available info
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

Symbol_fx <- function (x) {
  getBM(attributes=c('uniprot_gn_symbol','description'), 
        filters = 'pdb', 
        values = x, 
        mart = ensembl)
}

score_all <-
  score_all %>% 
  filter(n_rows != 0,
        !is.na(description)
         ) %>% 
  group_by(description) %>% 
  mutate(
    name = strsplit(description, split = "_")[[1]][1],
    name_struc = paste(name, ": ", Symbol_fx(name)[1,1]),
    interface_min = min(c_across(dG_A:observations), na.rm = TRUE),
    #mean_dG = mean(c_across(dG_AB:dG_AZ), na.rm = TRUE),
    sum_sec_interface = sum(c_across(dG_A:observations), na.rm = TRUE) - interface_min,
    type = fcase(
      substr(name, start = 1, stop = 3) == "6E9", "DHF",
      substr(name, start = 1, stop = 3) == "HA1", "HA13",
      substr(name, start = 1, stop = 3) != "6E9" & substr(name, start = 1, stop = 3) != "HA1", "DD"
    )
)

ggplot(
  data = score_all,
  aes(
    x = interface_min,
    y = sum_sec_interface,
    color = name_struc,
    shape = type
  )
)+
  geom_point(
    size = 3
  )

ggsave(
  filename = "dG_interface.pdf",
  scale = 0.5,
  units = "mm",
  height = 400,
  width = 600
)
  
