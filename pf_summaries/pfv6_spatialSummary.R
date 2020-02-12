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


theme_j <- function () {
  theme_bw(base_size=16) %+replace%
    theme(
      # font sizes and color
      panel.background  = element_blank(),
      plot.background   = element_rect(fill="transparent", colour=NA),
      plot.title        = element_text(size = rel(.85)),
      strip.background  = element_rect(fill="transparent", colour=NA),
      strip.text        = element_text(face="bold", size=rel(1)),
      axis.title        = element_text(size=rel(0.8)),
      axis.text         = element_text(size=rel(0.6), color="grey30"),
      # legend
      legend.title         = element_text(size=rel(0.8)),
      legend.text          = element_text(size=rel(0.6)),
      legend.background    = element_rect(fill="transparent", colour=NA),
      legend.key           = element_rect(fill="transparent", colour=NA),
      legend.justification = "top"
    )
}

theme_set(theme_j())


###############################################################################
# load and edit data
###############################################################################
pf6_dir <- "C:/Users/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Pfv6 Malaria"
pfv6_metadata <- read.delim(paste(pf6_dir, "Pf_6_MalariaGenMaterials/Pf_6_samples.txt", sep="/"), 
                            stringsAsFactors = F)
names(pfv6_metadata) <- tolower(names(pfv6_metadata))
pfv6_metadata$country <- ifelse(pfv6_metadata$country == "Viet Nam", "Vietnam", pfv6_metadata$country) 
# remove lab crosses from the dataset
pfv6_rmCont <- pfv6_metadata %>% filter(qc.pass == "True")


###############################################################################
# visualization parameters
###############################################################################
# attach continent to country names
pfv6_countries  <- unique(pfv6_rmCont$country)
pfv6_continents <- cbind.data.frame(country = pfv6_countries, 
                                    continent = countrycode(pfv6_countries, 'country.name', 'continent'),
                                    stringsAsFactors=F) 
pfv6_rmCross <- inner_join(pfv6_rmCont, pfv6_continents)

# set order for countries based on number of samples
country_descOrd <- pfv6_rmCross %>%
  group_by(country) %>% summarise(sample_total = n()) %>% 
  arrange(desc(sample_total))%>% .[["country"]] %>% as.character()

# set colors for years
year_colors <- setNames(c(brewer.pal(7, "Greys")[-7], rev(viridis_pal(option="C")(10))[-1]), 
                        seq(2001,2015))


###############################################################################
# plots
###############################################################################
#set output dir for figures
fig_dir = "C:/Users/jribado/OneDrive - IDMOD/Malaria/pfv6_summary"

# ceate summary data frame with counts per year and region
pfv6_counts <- pfv6_rmCross %>%
  group_by(continent, country, year) %>%
  summarise(year_freq = n())


# barplot by country
pfv6_country <- pfv6_counts %>%
  mutate(label = ifelse(year_freq > 10 , year_freq, NA)) %>%
  ggplot(aes(x=factor(country, levels=rev(country_descOrd)), y=year_freq, 
             label=label)) +
  geom_col(aes(y = year_freq, fill = as.character(year)), 
           position = position_stack(reverse = TRUE),
           width = .5) +
  geom_text(size = 2, color="white", position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(name = "Collection Year", values=year_colors) +
  labs(x="Country", y="Number of Pfv6 Samples") +
  theme_bw() + coord_flip() +
  theme(legend.justification = c(1, 0), 
        legend.position = c(1, 0),
        legend.box.margin=margin(c(30,20,20,20)))
ggsave("pfv6_countryBar.png", plot = pfv6_country, path = fig_dir, 
              width = 6, height = 7, units = c("in"), dpi = 300)


pfv6_year <-pfv6_counts %>%
  mutate(label = ifelse(year == min(year), as.character(country), NA_character_)) %>%
  ggplot(aes(x=year, y=year_freq,
             color=country, group=country)) +
  geom_point(pch=19, size=4.5, alpha=0.75) + 
  geom_line() +
  geom_text_repel(aes(label=label), box.padding = 0.25,
                  hjust = 0.5, vjust = -1) +
  #scale_color_manual(values=rev(country_colors)) +
  facet_wrap(~continent, ncol=2) +
  labs(x="Collection Year", y="Number of Samples") +
  theme_bw() + guides(color=FALSE)
ggsave("pfv6_yearLine.png", plot = pfv6_year, path = fig_dir, 
        width = 10, height = 8, units = c("in"), dpi = 300)


###############################################################################
# spatial graphs
###############################################################################
# expand year data to be in wide format
pfv6_siteCounts <- pfv6_rmCross %>%
  group_by(continent, country, site, year, lat, long) %>%
  summarise(year_freq = n()) %>% ungroup() %>%
  tidyr::spread(year, year_freq) %>%
  mutate_if(is.integer, ~replace(., is.na(.), 0))
# merge cooordinates with frequency data
pfv6_coordCounts <- pfv6_siteCounts %>%
  mutate(lat_fixed = ifelse(continent == "Africa", 0-lat, lat),
         nSite = rowSums(.[6:20]),
         nTotal = nSite/sum(nSite))


# map countries to continents for easy subsetting. some issues arise with 1-to-1 mapping, but none that affect the countries of interest.  
worldmap   <- map_data("world")
continents <- cbind.data.frame(
  region = unique(worldmap$region), 
  continent = countrycode(unique(worldmap$region), 'country.name', 'continent'),
  stringsAsFactors=F) 
worldmap <- inner_join(worldmap, continents)


# plotting function
scatterpie_plot <- function(continent_sub, xlim_min, xlim_max, ylim_min, ylim_max, prop_size){
  continent_map <- ggplot() + 
    geom_map(data = filter(worldmap, continent == continent_sub), 
             map  = filter(worldmap, continent == continent_sub), 
             aes(x=long, y=lat, map_id=region), col = "white", fill = "gray90") +
    xlim(c(xlim_min, xlim_max)) + ylim(c(ylim_min, ylim_max))
  continent_prop <- continent_map + 
    geom_scatterpie(aes(x=long, y=lat_fixed, group = site), 
                    data = filter(pfv6_coordCounts, continent == continent_sub), 
                    cols = colnames(pfv6_coordCounts)[6:20]) +
    scale_fill_manual(values=year_colors) + theme_void() 
  continent_propAll <- continent_map + 
    geom_scatterpie(aes(x=long, y=lat_fixed, group = site, r = nTotal*prop_size), 
                    data = filter(pfv6_coordCounts, continent == continent_sub), 
                    cols = colnames(pfv6_coordCounts)[6:20]) +
    scale_fill_manual(values=year_colors) + theme_void()
  
    return(list(a=continent_prop + theme(plot.margin = unit(c(1,1,1,1), "lines")),
                b=continent_propAll + theme(plot.margin = unit(c(1,1,1,1), "lines"))))
}

# determine the map limits per continent
coord_lims <- dplyr::group_by(pfv6_coordCounts, continent) %>%
  dplyr::summarise(xlim_min=min(long),
                   xlim_max=max(long),
                   ylim_min=min(lat),
                   ylim_max=max(lat))
  

# plot 
africa_plots  <- scatterpie_plot("Africa", -25, 55, 25, -20, 30)
asia_plots    <- scatterpie_plot("Asia", 90, 145, -8, 30, 30)
oceania_plots <- scatterpie_plot("Oceania", 140, 155, -12, 0, 100)
america_plots <- scatterpie_plot("Americas", -80, -72, -5, 7, 200)

# combined
scatterpi_comb <- ggarrange(asia_plots$a, asia_plots$b,
                            africa_plots$a, africa_plots$b,
                            oceania_plots$a, oceania_plots$b,
                            america_plots$a, america_plots$b, 
  ncol=2, nrow=4, common.legend = T, legend="left",
  labels=c("Asia", "", "Africa", "", "Oceania", "", "South America", ""))
ggsave("pfv6_scatterpieYear.png", plot = scatterpi_comb, path = fig_dir, 
       width = 8, height = 12, units = c("in"), dpi = 300)
