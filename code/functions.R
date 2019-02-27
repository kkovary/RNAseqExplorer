### Additional functions ###

geneSearch <- function(geneList){
  geneList
}



#geneList = c('Acp2', 'Nr1c1', 'COAB') %>% tolower()

#filter_all(test, tolower(GN_Syn) %in% geneList)

#test = gn %>% mutate(GN_Syn = str_split(GN_Syn, ' '))


#test %>% map_lgl(tolower(GN_Syn), ~ geneList)

#test %>% filter(map_lgl(tolower(GN_Syn), ~ which(. %in% geneList)))

test %>% map(GN_Syn, ~ !is.na(.))

any_not_na <- function(x) {
  !all(map_lgl(x, is.na))
}


dat_cleaned <- dat %>%
  rownames_to_column("ID") %>%
  group_by(ID) %>%
  nest() %>%
  filter(map_lgl(data, any_not_na)) %>%
  unnest() %>%
  select(-ID)




