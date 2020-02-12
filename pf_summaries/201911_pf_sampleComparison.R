# pfv6 sample summary plots
# author: jessica ribado
# 2019-11-04

###############################################################################
# install and load neccessary packages
###############################################################################
for(p in c('dplyr', 'ggplot2', 'ggrepel', 'countrycode', 
           'scatterpie', 'maps', 'ggpubr',
           'yarrr', "viridis", "RColorBrewer")){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}
theme_set(theme_bw())

# set colors for years
year_colors <- setNames(c(brewer.pal(7, "Greys")[-7], rev(viridis_pal(option="C")(10))[-1]), 
                        seq(2001,2015))



###############################################################################
# load and edit data
###############################################################################
# load Pf3k metadata 
pf3k_metadata <- read.delim(gzfile("P:/pf3k_v5_data/pf3k_release_5_metadata_20170804.txt.gz"), 
                            stringsAsFactors = F)
pf3k_rmCross <- pf3k_metadata %>% filter(!is.na(collection_year))


# load pfv6 metadata
pf6_dir <- "C:/Users/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Pfv6 Malaria"
pfv6_metadata <- read.delim(paste(pf6_dir, "Pf_6_samples.txt", sep="/"), 
                            stringsAsFactors = F)
names(pfv6_metadata) <- tolower(names(pfv6_metadata))
pfv6_rmCont <- pfv6_metadata %>% filter(qc.pass == "True")


# find unique samples in each study to create a new variable 
pf3k_only <- pf3k_rmCross$sample[!pf3k_rmCross$sample %in% pfv6_rmCont$sample]
pfv6_only <- pfv6_rmCont$sample[!pfv6_rmCont$sample %in% pf3k_rmCross$sample]

# combine datasets
pf3k_merge <- dplyr::select(pf3k_rmCross, sample, country, site, collection_year) %>%
  dplyr::mutate(study_inc = ifelse(sample %in% pf3k_only, "pf3k", "both")) %>%
  dplyr::rename("year"=collection_year)
pfv6_merge <-  dplyr::select(pfv6_rmCont, sample, country, site, year) %>%
  dplyr::mutate(study_inc = ifelse(sample %in% pfv6_only, "pfv6", "both"))
pf_combine <- rbind.data.frame(pf3k_merge, pfv6_merge, stringsAsFactors = F) %>% unique()
pf_combine$country <- ifelse(pf_combine$country == "Viet Nam", "Vietnam", pf_combine$country) 
continents <- cbind.data.frame(
                 country = unique(pf_combine$country), 
                 continent = countrycode(unique(pf_combine$country), 'country.name', 'continent'),
                 stringsAsFactors=F) 
pf_combine <- inner_join(pf_combine, continents)


###############################################################################
# visualization parameters
###############################################################################
# set order for countries based on number of samples
country_descOrd <- pf_combine %>%
  group_by(country) %>% summarise(sample_total = n()) %>% 
  arrange(desc(sample_total))%>% .[["country"]] %>% as.character()

# count sampels per country and year
pf_counts <- pf_combine %>%
  group_by(study_inc, continent, country, year) %>%
  summarise(year_freq = n())


# barplot by country
pf_country <- pf_counts %>% ungroup() %>%
  mutate(label = ifelse(year_freq > 10 , year_freq, NA),
         study_inc = factor(study_inc, levels = c("pf3k", "both", "pfv6")),
         continent = factor(continent, levels = c("Africa", "Asia", "Americas", "Oceania")))%>%
  ggplot(aes(x=factor(country, levels=rev(country_descOrd)), y=year_freq, 
             label=label)) +
  geom_col(aes(y = year_freq, fill = as.character(year)), 
           position = position_stack(reverse = TRUE),
           width = .5) +
  geom_text(size = 2, color="white", position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(name = "Collection Year", values=year_colors) +
  labs(x="Country", y="Number of samples") +
  facet_grid(continent~study_inc, scales="free", space="free") +
  theme_bw() + coord_flip() 


#ggsave("pfv6_countryBar.png", plot = pfv6_country, path = fig_dir, 
       width = 6, height = 7, units = c("in"), dpi = 300)





  