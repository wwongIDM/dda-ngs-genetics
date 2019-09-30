# pf3k sample summary plots
# author: jessica ribado
# 2019-09-23

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


###############################################################################
# load and edit data
###############################################################################
pf3k_metadata <- read.delim(gzfile("P:/pf3k_v5_data/pf3k_release_5_metadata_20170804.txt.gz"), 
                            stringsAsFactors = F)
# fix Vietnam spelling
pf3k_metadata$country <- ifelse(pf3k_metadata$country == "Viet Nam", "Vietnam", pf3k_metadata$country) 
# remove lab crosses from the dataset
pf3k_rmCross <- pf3k_metadata %>% filter(!is.na(collection_year))



###############################################################################
# visualization parameters
###############################################################################
# attach continent to country names
pf3k_countries  <- unique(pf3k_rmCross$country)
pf3k_continents <- cbind.data.frame(country = pf3k_countries, 
                                    continent = countrycode(pf3k_countries, 'country.name', 'continent'),
                                    stringsAsFactors=F) 
pf3k_rmCross <- inner_join(pf3k_rmCross, pf3k_continents)

# set order for countries based on number of samples
country_descOrd <- pf3k_rmCross %>%
  group_by(country) %>% summarise(sample_total = n()) %>% 
  arrange(desc(sample_total))%>% .[["country"]] %>% as.character()
  
# set colors for countries
country_colors <- setNames(c(piratepal(palette = "nemo"), piratepal(palette = "up"), piratepal(palette = "espresso")),
                           unique(pf3k_rmCross$country))
# set colors for years
year_colors <- c(brewer.pal(5, "Greys")[-5], rev(viridis_pal(option="C")(9))[-1])


###############################################################################
# plots
###############################################################################
#set output dir for figures
fig_dir = "C:/Users/jribado/OneDrive - IDMOD/pf3k_summary"

# ceate summary data frame with counts per year and region
pf3k_counts <- pf3k_rmCross %>%
  group_by(continent, country, collection_year) %>%
  summarise(year_freq = n())


# barplot by country
pf3k_country <- pf3k_counts %>%
  mutate(label = ifelse(year_freq > 10 , year_freq, NA)) %>%
  ggplot(aes(x=factor(country, levels=rev(country_descOrd)), y=year_freq, 
             label=label)) +
    geom_col(aes(y = year_freq, fill = as.character(collection_year)), 
           position = position_stack(reverse = TRUE),
           width = .5) +
    geom_text(size = 3, color="white", position = position_stack(vjust = 0.5)) + 
    scale_fill_manual(name = "Collection Year", values=year_colors) +
    labs(x="Country", y="Number of Pf3K Samples") +
    theme_bw() + coord_flip() +
    theme(legend.justification = c(1, 0), 
          legend.position = c(1, 0),
          legend.box.margin=margin(c(30,20,20,20)))
# ggsave("pf3k_countryBar.png", plot = pf3k_country, path = fig_dir, 
#       width = 6, height = 5, units = c("in"), dpi = 300)


# line graph by year
pf3k_year <-pf3k_counts %>%
  mutate(label = ifelse(collection_year == min(collection_year), as.character(country), NA_character_)) %>%
  ggplot(aes(x=collection_year, y=year_freq,
             color=country, group=country)) +
    geom_point(pch=19, size=4.5, alpha=0.75) + 
    geom_line() +
    geom_text_repel(aes(label=label), box.padding = 0.25,
                    hjust = 0.5, vjust = -1) +
    scale_color_manual(values=rev(country_colors)) +
    facet_grid(~continent) +
    labs(x="Collection Year", y="Number of Samples") +
    theme_bw() + guides(color=FALSE)
# ggsave("pf3k_yearLine.png", plot = pf3k_year, path = fig_dir, 
#        width = 6, height = 4, units = c("in"), dpi = 300)


###############################################################################
# spatial graphs
###############################################################################
# this was unexpectedly more difficult than expected to do programatically  without providing Google maps a CC to use their API, or finding older packages with new API restrictions that no longer work with the specified functions. One option:
# http://datapages.com/gis-map-publishing-program/gis-open-files/global-framework/global-heat-flow-database/shapefiles-list
# africa_shp_path = "C:/Users/jribado/OneDrive - IDMOD/pf3k_summary/Africa_Shapefile/Africa.shp"
# africa_shp <- rgdal::readOGR(dsn = africa_shp_path)
# africa <- fortify(as.data.frame(africa_shp))

# same issue with getting site (lat,lon) data with API restrictions or depreciation of different open source tools. Manually curate the list of coordinates for Google since there are only 26 - may be needed downstream for genomic data plotting too. Note: coordinate for Kolle, Mali are not available; usin Bamako coordinates as a proxy. 
sites_df <- dplyr::select(pf3k_rmCross, country, site) %>% unique()
# write.table(sites_df, paste(fig_dir, "pf3k_sites_coords.txt", sep="/"), sep="\t", quote=F, row.names = F)
coords_df <- read.delim(paste(fig_dir, "pf3k_sites_coords.txt", sep="/"), sep="\t", stringsAsFactors = F)


# expand year data to be in wide format
pf3k_siteCounts <- pf3k_rmCross %>%
  group_by(continent, country, site, collection_year) %>%
  summarise(year_freq = n()) %>% ungroup() %>%
  tidyr::spread(collection_year, year_freq) %>%
  mutate_if(is.integer, ~replace(., is.na(.), 0))
# merge cooordinates with frequency data
pf3k_coordCounts  <- inner_join(coords_df, pf3k_siteCounts) %>%
  mutate(nSite = rowSums(.[6:16]),
         nTotal = nSite/sum(nSite))


# map countries to continents for easy subsetting. some issues arise with 1-to-1 mapping, but none that affect the countries of interest.  
worldmap   <- map_data("world")
continents <- cbind.data.frame(
  region = unique(worldmap$region), 
  continent = countrycode(unique(worldmap$region), 'country.name', 'continent'),
  stringsAsFactors=F) 
worldmap <- inner_join(worldmap, continents)


# map the year percentages per site both as total in site and as overall poroprtion of total samples
# africa
africa <- ggplot() + 
  geom_map(data = filter(worldmap, continent == "Africa"), 
           map  = filter(worldmap, continent == "Africa"), 
           aes(x=long, y=lat, map_id=region), col = "white", fill = "gray90") +
  xlim(c(-20,40)) + ylim(c(20,-20))
africa_prop    <- africa + 
  geom_scatterpie(aes(x=longitude, y=flip, group = site), 
                  data = filter(pf3k_coordCounts, continent =="Africa"), 
                  cols = colnames(pf3k_coordCounts)[6:16]) +
  scale_fill_manual(values=year_colors) + theme_void() 
africa_propAll <- africa + 
  geom_scatterpie(aes(x=longitude, y=flip, group = site, r = nTotal*25), 
                  data = filter(pf3k_coordCounts, continent =="Africa"), 
                  cols = colnames(pf3k_coordCounts)[6:16]) +
  scale_fill_manual(values=year_colors) + theme_void()

# asia
asia <- ggplot() + 
  geom_map(data = filter(worldmap, continent == "Asia"), 
           map  = filter(worldmap, continent == "Asia"), 
           aes(x=long, y=lat, map_id=region), col = "white", fill = "gray90") +
  xlim(c(90,110)) + ylim(c(10,25))
asia_prop    <- asia + 
  geom_scatterpie(aes(x=longitude, y=latitude, group = site, r = .5), 
                  data = filter(pf3k_coordCounts, continent =="Asia"), 
                  cols = colnames(pf3k_coordCounts)[6:16]) +
  scale_fill_manual(values=year_colors) + theme_void()
asia_propAll <- asia + 
  geom_scatterpie(aes(x=longitude, y=latitude, group = site, r = nTotal*15), 
                  data = filter(pf3k_coordCounts, continent =="Asia"), 
                  cols = colnames(pf3k_coordCounts)[6:16]) +
  scale_fill_manual(values=year_colors) + theme_void()


# combine plots
scatterpi_comb <- ggarrange(
  africa_prop + theme(plot.margin = unit(c(1,1,1,1), "lines")) + 
    labs(fill = "Collection \nYear") , 
  africa_propAll + theme(plot.margin = unit(c(1,1,1,1), "lines")), 
  asia_prop + theme(plot.margin = unit(c(1,1,1,1), "lines")), 
  asia_propAll + theme(plot.margin = unit(c(1,1,1,1), "lines")), 
  ncol=2, nrow=2, common.legend = TRUE, legend="left")
ggsave("pf3k_scatterpieYear.png", plot = scatterpi_comb, path = fig_dir, 
               width = 8, height = 6, units = c("in"), dpi = 300)

                    