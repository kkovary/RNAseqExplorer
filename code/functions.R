### Additional functions ###

geneSearch <- function(searchList, snyList, rnaSeqData){
  # The goal of this function is to take in a list of genes
  # of interest and search a gene synonym database for
  # a match with a cannonical name as well as the searched name.
  
  df = snyList %>% filter(tolower(GN_Syn) %in% tolower(searchList)) %>%
    mutate(match = GeneName == GN_Syn, plotName = ifelse(match, GeneName, 
                                                           paste0(GN_Syn,' (',GeneName,')')))
  
  df = df %>% inner_join(rnaSeqData, ., by = 'GeneName') %>%
    select(plotName, NC_0_1:PPARG_6_3) %>% dplyr::rename(GeneName = plotName) %>% gather("Sample", "TPM", 2:ncol(.)) %>% 
    separate(Sample, into = c("siRNA", "Day", "Replicate"), sep = "\\_") %>%
    mutate(siRNA = paste0('si', siRNA), Day = as.numeric(Day))
  
  return(df)
}

convertUniprot <- function(csvFilepath){
  # This function reads in a csv file of gene/protein names
  # and organizes them to be used with this app to be read
  # in as 'geneSyns'.
  require(readxl)
  
  read_xlsx(csvFilepath) %>% select(Entry, `Gene names  (primary )`, `Gene names`) %>% 
    rename(GeneName = `Gene names  (primary )`, GN_Syn = `Gene names`) %>%
    mutate(GN_Syn = strsplit(GN_Syn, ' ')) %>% 
    unnest() %>% filter(!duplicated(tolower(GN_Syn))) %>% 
    write_csv(path = paste0(getwd(),'/Data/GeneNames.csv'))
}



