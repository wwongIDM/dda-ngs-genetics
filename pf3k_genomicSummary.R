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

#set output dir for figures
fig_dir = "C:/Users/jribado/OneDrive - IDMOD/pf3k_summary"

# set colors for countries
country_colors <- setNames(c(piratepal(palette = "nemo"), 
                             piratepal(palette = "up"), 
                             piratepal(palette = "espresso")),
                           unique(pf3k_rmCross$country))
# set colors for years
year_colors <- c(brewer.pal(5, "Greys")[-5], rev(viridis_pal(option="C")(9))[-1])

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
# highlight genomic data quality metrics from metadata file
###############################################################################
# genome wide average coverage
coverage <- pf3k_rmCross %>%
  ggplot(aes(x=bases_mapped/1e9, y=mean_coverage)) +
  scale_color_manual(values=country_colors) + 
  geom_hline(yintercept=50, linetype="dashed", color = "grey60", size=1) +
  geom_hline(yintercept=10, linetype="dashed", color = "grey15", size=1) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(low="#004C3FFF", high="#DFF2F1FF") +
  labs(x="Bases mapping to reference genome (Millions)",
       y="Mean genome coverage") +
  guides(fill=guide_legend(title="Sample count")) +
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.box.margin=margin(c(10,0,0,10)))
ggsave("pf3k_genomeCoverage.png", plot = coverage, path = fig_dir, 
       width = 5, height = 4, units = c("in"), dpi = 300)

###############################################################################
# highlight (high level) genetic diversity with vcf tools output
###############################################################################
# read in and edit files
vcf_dir <- "P:/pf3k_v5_data/vcf_stats"
vcf_files <- dir(vcf_dir, all.files = FALSE, full.names = T, pattern="vcfSubset_") 

# load vcf outputs
setwd(vcf_dir)
vcf_list <- lapply(vcf_files, data.table::fread)
names(vcf_list) <- gsub(".*Subset_|\\.txt.*", "", vcf_files)
# read in headers
vcf_headers <- read.delim(paste(vcf_dir, "bcftools_headers.txt", sep="/"), 
                          header=F, stringsAsFactors = F) %>%
  mutate_if(is.character, ~gsub(".*\\#\ |.*\\]", "", .)) %>%
  dplyr::arrange(V1)
# for loops are not ideal, but is th easiest way to fix the column names 
for (i in 1:length(gsub(".*Subset_|\\.txt.*", "", vcf_files))){
  vcf_list[[i]] <- tidyr::separate(vcf_list[[i]], V1, c("chr", "stat"), sep=":") %>%
    dplyr::mutate(chr = gsub(".*Pf3D7_|\\_v3.*", "", chr)) 
  colnames(vcf_list[[i]]) <- c("chr", as.character(vcf_headers[i,])[1:ncol(vcf_list[[i]])-1])
}

# summary stats
summ_stats <- vcf_list$SN %>% 
  dplyr::filter(grepl("number of SNPs|number of others", key)) %>%  
  dplyr::filter(!grepl("API|stats", chr)) %>%
  ggplot (aes(x=chr, y=value/1e5, fill=key)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Genetic category",
                    values=RColorBrewer::brewer.pal(3, "Dark2"),
                    labels=c("Structural variants, \ncomplex rearrangements, \netc.",
                             "Single or multi-nucleotide \nvariants, indels, etc."))+
  labs(x="Chromosome", y="Number of sites (10K)")+
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.box.margin=margin(c(10,0,0,10)))
ggsave("pf3k_summStats.png", plot = summ_stats, path = fig_dir, 
       width = 6, height = 4, units = c("in"), dpi = 300)


# snps per allele frequency bin
snp_af <- vcf_list$AF %>%  
  dplyr::filter(!grepl("API|stats", chr)) %>%
  dplyr::mutate(af_bin = cut(`allele frequency`, breaks=(0:20)*.05)) %>%
  dplyr::filter(!is.na(af_bin)) %>%
  dplyr::group_by(chr, af_bin) %>%
  dplyr::summarise(total_snps = sum(`number of SNPs`)) %>%
  ggplot(aes(x=af_bin, y=total_snps, color=chr, group=chr)) +
  # ggbeeswarm::geom_quasirandom(alpha=0.5,method = "smiley") +
  geom_point() +
  geom_line() +
  # facet_grid(chr~.) +
  scale_y_log10() +
  labs(x="Allele Frequency Bin", y="Number of SNPs") + 
  theme(axis.text.x = element_text(angle = 90))
ggsave("pf3k_snpsAFBin.png", plot = snp_af, path = fig_dir, 
       width = 6, height = 4, units = c("in"), dpi = 300)

# insertion length density 
insert_dist <- vcf_list$IDD %>% 
  dplyr::filter(!grepl("API|stats", chr)) %>%
  ggplot(aes(x=`length (deletions negative)`, y=`count`/1e4, fill=chr)) +
  geom_density(stat = "identity") +
  facet_grid(chr~.) +
  labs(x="Insertion length", y="Count (10K)")+
  guides(fill=FALSE)
ggsave("pf3k_insertDist.png", plot = insert_dist, path = fig_dir, 
       width = 6, height = 12, units = c("in"), dpi = 300)


# transitions and transversions
trans_heatmap <- vcf_list$ST %>% 
  dplyr::filter(!grepl("API|stats", chr)) %>%
  tidyr::separate(type, c("original", "change"), sep=">", remove=F) %>%
  dplyr::mutate(change=ifelse(type %in% c("C>T", "T>C", "A>G", "G>A"), "transition", "transversion")) %>%
  ggplot(aes(y=type, x=chr, fill=scale(count))) +
  geom_tile() + 
  scale_fill_gradient2(name="Count Z-score", low="blue", high="red") +
  facet_grid(change~., space="free", scales="free")+ 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Conversion type", x="Chromosome") +
  theme(legend.position = "top")
ggsave("pf3k_basechangeHeatmap.png", plot = trans_heatmap, path = fig_dir, 
       width = 3, height = 4, units = c("in"), dpi = 300)