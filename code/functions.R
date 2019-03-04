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

volPlotDataFun <- function(numerator, denominator, data, pvalueCut, foldChangeCut){
  # This function takes in the RNAseq data, and calculates the fold change
  # and pvalue from two conditions (numerator and denominator) and labels
  # which genes are significant or not based on the user defined pvalue and
  # fold change cutoffs.
  
  data %>% filter(Condition %in% c(numerator, denominator)) %>%
    group_by(GeneName) %>% 
    summarise(foldChange = mean(TPM[Condition == numerator], na.rm = T) / mean(TPM[Condition == denominator], na.rm = T),
              pvalue = t.test(TPM[Condition == numerator], TPM[Condition == denominator], na.rm = T)$p.value) %>%
    mutate(Hit = ifelse(pvalue <= pvalueCut & abs(log2(foldChange)) >= log2(foldChangeCut), TRUE, FALSE))
  
}

convertUniprot <- function(csvFilepath){
  # This function reads in a csv file of gene/protein names
  # and organizes them to be used with this app to be read
  # in as 'geneSyns'.
  require(readxl)
  
  read_xlsx(csvFilepath) %>% select(Entry, `Gene names  (primary )`, `Gene names`) %>% 
    rename(GeneName = `Gene names  (primary )`, GN_Syn = `Gene names`) %>%
    mutate(GeneName = strsplit(GeneName, '; ')) %>% unnest() %>% mutate(GN_Syn = strsplit(GN_Syn, ' ')) %>% 
    unnest() %>% filter(!duplicated(tolower(GN_Syn))) %>% 
    
    write_csv(path = paste0(getwd(),'/Data/GeneNames.csv'))
}

