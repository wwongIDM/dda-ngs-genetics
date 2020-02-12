# pfv6 fws plot

################################################################################
for(p in c('data.table', 'dplyr', 'ggplot2')){
  if(!p %in% installed.packages()[,1]){
    install.packages(p)
    library(p, character.only = TRUE)
  }
  library(p, character.only = TRUE)
}

################################################################################
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


pf6_dir <- "C:/Users/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Pfv6 Malaria/Pf_6_MalariaGenMaterials"
pf6_meta <- fread(paste(pf6_dir, "Pf_6_samples.txt", sep="/"))
pf6_fws  <- fread(paste(pf6_dir, "Pf_6_fws.txt", sep="/"))
pf6_comb <- inner_join(pf6_meta, pf6_fws)
pf6_comb <- dplyr::mutate(pf6_comb, Population = factor(Population, 
                                                        levels = c("WAF", "EAF", "CAF", "WSEA", "ESEA", "SAS", "SAM", "OCE")))

fws_plot <- pf6_comb %>%
  ggplot(aes(x=Population, y=Fws)) +
  ggbeeswarm::geom_quasirandom(alpha=0.5, pch=19) +
  geom_hline(yintercept=0.95, colour="#990000", linetype="dashed") +
  labs(x="Pfv6 Population", y="Within-host diversity\n(FWS)")
ggsave("pfv6_fwsScatter.png", plot=fws_plot, path="C:/Users/jribado/Dropbox (IDM)/Data, Dynamics, and Analytics Folder/Projects/Pfv6 Malaria/plots", width=5, height=5, units = c("in"), dpi = 300)
