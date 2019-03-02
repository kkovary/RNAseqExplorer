### Additional functions ###

geneSearch <- function(searchList, geneList){
  # The goal of this function is to take in a list of genes
  # of interest and search a gene synonym database for
  # a match with a cannonical name as well as the searched name.
  
  geneList %>% filter(tolower(GN_Syn) %in% tolower(searchList))
}

nameMatch <- function(x){
  # The goal of this function is to take the matches from a gene
  # search list and a synonym list, match them with the RNA-seq
  # data, and output a dataframe with the searched for names
  # along with the cannonical name.
  
  
}

# Searches the GN_Syn column for all matching GNs
# and returns the filtered df

