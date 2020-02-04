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


# Dynamics categorizer
cat_heatmap <- function(data = norm, cat = cat_quant_long, d0_siRNA_ = NULL, 
                        dm_siRNA_ = NULL, dm_NC_ = NULL, dm_PPARG_ = NULL, 
                        ins_siRNA_ = NULL, ins_NC_ = NULL, ins_PPARG_ = NULL,
                        output = 'heatmap'){
  if(!is.null(d0_siRNA_)){
    cat <- filter(cat, d0_siRNA %in% d0_siRNA_)
  }
  if(!is.null(dm_siRNA_)){
    cat <- filter(cat, dm_siRNA %in% dm_siRNA_)
  }
  if(!is.null(dm_NC_)){
    cat <- filter(cat, dm_NC %in% dm_NC_)
  }
  if(!is.null(dm_PPARG_)){
    cat <- filter(cat, dm_PPARG %in% dm_PPARG_)
  }
  if(!is.null(ins_siRNA_)){
    cat <- filter(cat, ins_siRNA %in% ins_siRNA_)
  }
  if(!is.null(ins_NC_)){
    cat <- filter(cat, ins_NC %in% ins_NC_)
  }
  if(!is.null(ins_PPARG_)){
    cat <- filter(cat, ins_PPARG %in% ins_PPARG_)
  }
  data = data %>% filter(GeneName %in% cat$GeneName) %>%
    group_by(GeneName, siRNA, day) %>% 
    summarise(meanTPM = mean(tpm, na.rm = T)) %>%
    group_by(GeneName) %>% mutate(normTPM = meanTPM / meanTPM[day == 0 & siRNA == 'NC'])
  
  # font_size <- 10 - 0.04*nrow(cat)
  # if(font_size < 0){
  #   font_size <- 0.1
  # }
  
  brbg_hcl <- colorspace::diverging_hcl(11,
                                        h = c(180, 50), c = 80, l = c(20, 95), power = c(0.7, 1.3))
  
  switch(output,
         'heatmap' = if(nrow(data) !=0){
           ggplot(data, aes(x=as.factor(day), y=GeneName, fill=log2(normTPM))) + 
             geom_tile() + 
             #scale_fill_viridis(name="TPM") + 
             scale_fill_gradient2(low = '#2166ac',mid = '#f7f7f7', high = '#b2182b', limits = c(-4,4), oob = squish) +
             #scale_fill_continuous(palette = 'Blue-Red') +
             coord_equal() + 
             facet_wrap(~siRNA, ncol=2) + theme_classic() +
             xlab('Day') + theme(axis.text.y = element_text(size = 10),
                                 axis.text.x = element_text(angle = 90, size = 6))
         } else{
           warning('no genes in selected range')
         },
         'cat' = cat,
         'data' = data)
  
  
}
