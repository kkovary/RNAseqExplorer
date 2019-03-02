### Additional functions ###

geneSearch <- function(searchList){
  # The goal of this function is to take in a list of genes
  # of interest and search a gene synonym database for
  # a match with a cannonical name as well as the searched name.
  
  df = geneSyns %>% filter(tolower(GN_Syn) %in% tolower(searchList)) %>%
    mutate(match = GeneName == GN_Syn, plotName = ifelse(match, GeneName, 
                                                           paste0(GN_Syn,' (',GeneName,')')))
  
  df = df %>% inner_join(normalized_data_genelevel_tpm, ., by = 'GeneName') %>%
    select(plotName, NC_0_1:PPARG_6_3) %>% gather("Sample", "TPM", 2:ncol(.)) %>% 
    separate(Sample, into = c("siRNA", "Day", "Replicate"), sep = "\\_") %>%
    mutate(siRNA = paste0('si', siRNA), Day = as.numeric(Day)) %>%
    filter(siRNA %in% input$siRNA, Day %in% input$time)
  
  return(df)
}

