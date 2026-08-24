library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)
library(data.table)
library(ggrepel)

# folder containing the kraken reports
kraken_dir <- "~/HPC_data/iron_reducing_bacteria/results/krakenuniq_collapsed/"

# list files like N1_kraken_report_standard, N2_kraken_report_standard, etc.
files <- list.files(
  path = kraken_dir,
  pattern = "^N[0-9]+_kraken_report_standard$",
  full.names = TRUE
)

# function to read one kraken standard report
read_kraken_report <- function(filex) {
  # Kraken standard report usually has 6 columns:
  # percent, clade_reads, taxon_reads, rank, taxid, name
  # names often have leading spaces due to hierarchy
  df <- fread(filex)
  colnames(df) = c("percent","reads","taxReads","kmers","dup","cov","taxID","rank","taxName")
  sample_name <- basename(filex) %>%
    str_remove("_kraken_report_standard$")
  
  df %>%
    mutate(
      sample = sample_name,
      name = str_trim(taxName)
    )
}

# read all files
kraken_all0 <- map_dfr(files, read_kraken_report)

# assign groups
kraken_all <- kraken_all0 %>%
  mutate(
    sample_num = as.integer(str_extract(sample, "[0-9]+")),
    group = case_when(
      sample_num %in% 1:6 ~ "spear",
      sample_num %in% 7:11 ~ "soil",
      TRUE ~ NA_character_
    )
  )

target_taxa <- tibble(
  genus = c("Rhodoferax", "Ferribacterium", "Desulfuromonas", "Desulfuromonas",
            "Geobacter", "Geobacter", "Geothrix", "Methylomonas",
            "Shewanella", "Shewanella"),
  species = c("ferrireducens", "limneticum", "acetoxidans", "spp.",
              "sulfurreducens", "metallireducens", "fermentans", "spp.",
              "oneidensis", "algae")
) %>%
  mutate(taxon_label = paste(genus, species))

kraken_targets <- kraken_all %>%
  mutate(
    word_count = str_count(name, "\\S+"),
    genus_from_name = word(name, 1),
    species_from_name = if_else(word_count >= 2, word(name, 2), NA_character_),
    taxon_label = case_when(
      # species-level rows: "Genus species"
      word_count >= 2 ~ paste(genus_from_name, species_from_name),
      # genus-level rows: "Genus" -> convert to "Genus spp."
      word_count == 1 ~ paste(genus_from_name, "spp."),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(taxon_label %in% target_taxa$taxon_label)

unique(kraken_targets$taxon_label)

# --------------------------------------------------
# choose abundance measure
# --------------------------------------------------

# Here I use direct taxon reads normalized by total direct taxon reads in the report.
# This gives a relative abundance-like value within each sample.

plot_df <- kraken_targets %>%
  mutate(
    rel_abundance = percent, # / 100,
    taxon_label = factor(taxon_label, levels = target_taxa$taxon_label)
  ) %>%
  select(sample, group, taxon_label, reads, rel_abundance)

# add zeros for taxa absent in a sample
plot_df_full <- plot_df %>%
  left_join(
    kraken_all %>%
      distinct(sample, group),
    by = "sample"
  ) %>%
  mutate(
    group = coalesce(group.x, group.y)
  ) %>%
  select(sample, group, taxon_label, reads, rel_abundance) %>%
  distinct()

# preserve sample order N1, N2, ...
sample_levels <- paste0("N", 1:11)

plot_df_full <- plot_df_full %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    group = factor(group, levels = c("spear", "soil"))
  )


unique(plot_df_full$taxon_label)


per_sample = ggplot(plot_df_full, aes(x = sample, y = rel_abundance, fill = group)) +
  geom_col() +
  facet_wrap(~ taxon_label, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Relative abundance",
    # title = "Selected KrakenUniq taxa across metal and soil samples"
  ) +
  scale_fill_manual(values = c("darkgreen","lightblue"),breaks = c("spear","soil"), labels=c("Spears","Soil")) + 
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "italic"),
    legend.title=element_blank()
    )

img_dir = "/Users/jrodriguez/Projects/iron_reducing_bacteria/"
# ggsave(paste0(img_dir,"per_sample_barplot_abundance_targets.pdf"), per_sample)

summary_df <- plot_df_full %>%
  group_by(group, taxon_label) %>%
  summarise(
    mean_abundance = mean(rel_abundance, na.rm = TRUE),
    sd_abundance = sd(rel_abundance, na.rm = TRUE),
    n = n(),
    se_abundance = sd_abundance / sqrt(n),
    .groups = "drop"
  )

grouped_abundance = ggplot(summary_df, aes(x = taxon_label, y = mean_abundance, fill = group)) +
  geom_col(position = position_dodge(width = 0.9)) +
  # geom_errorbar(
  #   aes(
  #     ymin = pmax(mean_abundance - se_abundance, 0),
  #     ymax = mean_abundance + se_abundance
  #   ),
  #   position = position_dodge(width = 0.8),
  #   width = 0.2
  # ) +
  theme_bw() +
  labs(x = NULL,y = "Mean relative abundance") +
  scale_fill_manual(values = c("darkgreen","lightblue"),breaks = c("spear","soil"),labels=c("Spears","Soil")) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    text = element_text(size = 20)
  )

# ggsave(paste0(img_dir,"grouped_barplot_abundance_targets.pdf"), grouped_abundance)


#####################
##### TEST DEBUG ZONE
# target_species <- c("Rhodoferax ferrireducens","Ferribacterium limneticum","Desulfuromonas acetoxidans","Geobacter sulfurreducens","Geobacter metallireducens","Geothrix fermentans","Shewanella oneidensis","Shewanella algae")
target_species <- c("Rhodoferax ferrireducens","Geobacter sulfurreducens","Geobacter metallireducens","Methanococcoides methylutens")

kraken_all$code = paste0(kraken_all$name," ",kraken_all$sample)

kraken_breadth0 = kraken_all %>% filter(rank == "species")

lab_df <- kraken_breadth0 %>%
  dplyr::filter(taxName %in% target_species)

kraken_breadth <- kraken_breadth0 %>%
  mutate(
    highlight = ifelse(taxName %in% lab_df$taxName, "target", "other")
  )


breadth_plot = ggplot() +
  geom_point(
    data = kraken_breadth %>% dplyr::filter(!(taxName %in% lab_df$taxName)),
    aes(x = log10(kmers), y = log10(reads)),
    color = "grey90",
    alpha = 0.5,
    size = 1
  ) +
  geom_hline(yintercept = 2, lty = "dashed", col="lightgreen",lwd = 1.5) +
  geom_vline(xintercept = 3, lty = "dashed", col="lightgreen",lwd = 1.5) +
  geom_point(
    data = lab_df,
    aes(x = log10(kmers), y = log10(reads)),
    color = "red",
    size = 3
  ) +
  geom_label_repel(
    data = lab_df,
    aes(x = log10(kmers), y = log10(reads), label = code),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() + 
  theme(
    legend.title = element_blank(),
    text = element_text(size = 20))

breadth_plot

ggsave(paste0(img_dir,"breadth_cov_plot.pdf"), breadth_plot)

