library(tidyverse)
library(viridis)
library(scales)
library(gtools)
# Categorize expression dynamics

# Based on conditions
d0 <- c('+','-','=')

dm <- expand.grid(d0,d0) %>%
  unite(col = 'dmi', Var1:Var2,sep = '') %>% unlist() %>% as.vector()

ins <- dm

expand.grid(d0,dm,ins)

# if you want to include filtering for multiple conditions (i.e. + & =)
#a <- sapply(1:7, function(x) nrow(permutations(3, x, repeats.allowed = T))) 

a <- sapply(1:7, function(x) nrow(permutations(3, x, repeats.allowed = T)))
b <- sapply(1:7, function(x) nrow(combinations(7, x, repeats.allowed = F)))
sum(a*b)

sum(sapply(1:7, function(x) nrow(permutations(3, x, repeats.allowed = T))))

# Data import
norm <- read_csv('~/Documents/GitHub/ShinyApps/RNAseqExplorer/data/normalized_data_genelevel_tpm.csv')
norm <- norm %>% pivot_longer(NC_0_1:PPARG_6_3, names_to = 'sample', values_to = 'tpm') %>%
  separate(sample, into = c('siRNA','day','rep'), sep = '_', remove = F) %>%
  mutate(day = as.numeric(day))

keep <- norm %>% group_by(GeneName) %>% summarise(low = max(tpm) < 2) %>%
  filter(!low)
norm <- norm %>% filter(GeneName %in% keep$GeneName)

# Identify categories in data
dm <- c(0.5,1,2)
ins <- c(3,4,5,6)
cat_quant <- norm %>% group_by(GeneName) %>% 
  summarise(
    #D0 Section
    d0_siRNA_fc = mean(tpm[siRNA == 'NC' & day == 0], na.rm = T) / mean(tpm[siRNA == 'PPARG' & day == 0], na.rm = T),
    d0_siRNA_pval = t.test(tpm[siRNA == 'NC' & day == 0], tpm[siRNA == 'PPARG' & day == 0])$p.value,
    #DM Section
    dm_siRNA_fc = mean(tpm[siRNA == 'NC' & day %in% dm], na.rm = T) / mean(tpm[siRNA == 'PPARG' & day %in% dm], na.rm = T),
    dm_siRNA_pval = t.test(tpm[siRNA == 'NC' & day %in% dm], tpm[siRNA == 'PPARG' & day %in% dm])$p.value,
    dm_NC_fc = mean(tpm[siRNA == 'NC' & day %in% dm], na.rm = T) / mean(tpm[siRNA == 'NC' & day == 0], na.rm = T),
    dm_NC_pval = t.test(tpm[siRNA == 'NC' & day %in% dm], tpm[siRNA == 'NC' & day == 0])$p.value,
    dm_PPARG_fc = mean(tpm[siRNA == 'PPARG' & day %in% dm], na.rm = T) / mean(tpm[siRNA == 'PPARG' & day == 0], na.rm = T),
    dm_PPARG_pval = t.test(tpm[siRNA == 'PPARG' & day %in% dm], tpm[siRNA == 'PPARG' & day == 0])$p.value,
    #Ins Section
    ins_siRNA_fc = mean(tpm[siRNA == 'NC' & day %in% ins], na.rm = T) / mean(tpm[siRNA == 'PPARG' & day %in% ins], na.rm = T),
    ins_siRNA_pval = t.test(tpm[siRNA == 'NC' & day %in% ins], tpm[siRNA == 'PPARG' & day %in% ins])$p.value,
    ins_NC_fc = mean(tpm[siRNA == 'NC' & day %in% ins], na.rm = T) / mean(tpm[siRNA == 'NC' & day == dm], na.rm = T),
    ins_NC_pval = t.test(tpm[siRNA == 'NC' & day %in% ins], tpm[siRNA == 'NC' & day == dm])$p.value,
    ins_PPARG_fc = mean(tpm[siRNA == 'PPARG' & day %in% ins], na.rm = T) / mean(tpm[siRNA == 'PPARG' & day == dm], na.rm = T),
    ins_PPARG_pval = t.test(tpm[siRNA == 'PPARG' & day %in% ins], tpm[siRNA == 'PPARG' & day == dm])$p.value
  ) %>% pivot_longer(cols = d0_siRNA_fc:ins_PPARG_pval, names_to = 'stat', values_to = 'value') %>%
  separate(stat, into = c('time','comparison','stat')) %>%
  unite('condition', time:comparison) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(p.adj = p.adjust(pval, 'BH'),
         exp_cat = ifelse(
           fc >= 1.5 & p.adj <= 0.05,
           '+',
           ifelse(
             fc <= 2/3 & p.adj <= 0.05,
             '-',
             '='
           )
         )
  )

cat_quant_long <- cat_quant %>% select(-fc, -pval, -p.adj) %>%
  pivot_wider(names_from = condition, values_from = exp_cat) %>% drop_na()

cat_quant_long_sum <- cat_quant %>% select(-fc, -pval, -p.adj) %>%
  pivot_wider(names_from = condition, values_from = exp_cat) %>% drop_na() %>%
  unite(col = 'category', d0_siRNA:ins_PPARG) %>% select(category) %>%
  unlist %>% table()

filter(cat_quant_long, 
       d0_siRNA == '=',
       dm_siRNA == '=',
       dm_NC == '+',
       ins_NC == '-',
       ins_siRNA == '='
)

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
  
  font_size <- 10 - 0.04*nrow(cat)
  if(font_size < 0){
    font_size <- 0.1
  }
  
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
             xlab('Day') + theme(axis.text.y = element_text(size = font_size),
                                 axis.text.x = element_text(angle = 90))
         } else{
           warning('no genes in selected range')
         },
         'cat' = cat,
         'data' = data)
  
  
}

cat_heatmap( 
  d0_siRNA_ = NULL, 
  dm_siRNA_ = '+', 
  dm_NC_ = '-',
  dm_PPARG_ = NULL,
  ins_siRNA_ = '+',
  ins_NC_ = '+',
  ins_PPARG_ = NULL
)


norm_func <- function(x, metric){
  # Function that will normalize all values by the first time point
  zero <- filter(x, day == '0')$tpm
  x %>% filter(day != '0') %>% group_by(day) %>% 
    summarise(fc = mean(tpm, na.rm = T) / mean(zero, na.rm = T),
              pval = t.test(tpm, zero, na.rm = T)$p.value)
  # switch(metric,
  #        'fc' = return(c(NA,x$fc)),
  #        'pval' = return(c(NA,x$pval))
  #        )
}

test <- norm %>% group_by(GeneName, siRNA) %>% nest() %>% 
  rowwise() %>% mutate(fc = list(norm_func(data, 'fc'))) %>%
  ungroup() %>% dplyr::select(-data) %>% unnest(cols = fc)

test <- test %>% mutate(p.adj = p.adjust(pval, 'BH'))

test[is.na(test$fc),'fc'] <- 0
test[is.na(test$p.adj),'p.adj'] <- 1
test[is.na(test$pval),'pval'] <- 1

test$dir <- 'unchanged'
test[test$fc >= 2 & test$p.adj < 0.05, 'dir'] <- 'up'
test[test$fc <= 0.5 & test$p.adj < 0.05, 'dir'] <- 'down'

test2 <- test %>% dplyr::select(GeneName, siRNA, day, dir) %>%
  group_by(GeneName, siRNA) %>% nest() %>% rowwise() %>%
  filter(sum(is.na(data$dir)) == 0) %>% ungroup() %>%
  unnest(cols = data) %>%
  pivot_wider(names_from = day, values_from = dir) %>%
  unite(pattern, `0.5`:`6`, sep = '_')

length(unique(test2$pattern))
table(test2$pattern)
filter(test2, pattern == 'up_up_up_up_up_up_up')


