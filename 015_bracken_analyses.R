###
# 0. SCRIPT USED TO ANALYZE THE RESULTS FROM KRAKENUNIQ ####
###

library(tidyverse)
library(phyloseq)
library(ape)
library(data.table)
library(microbiome)
library(RColorBrewer)
library(scales)
library(viridis)
library(microViz)
library(ggpubr)
library(ggrepel)
library(vegan)

# library(pheatmap)

options(scipen = 10)
setwd("~/HPC_data/iron_reducing_bacteria//results/bracken/")

# CAN BE TRIED ALSO WITH THE ABUNDANCE, SO WE CAN FILTER BETTER THE %.
ind = "counts" # "abund"

# 0. SETUP ####
# OTU count table
level = "G" # "Genus" # "Species"
# Read abundances bracken 
rawd = fread(paste0("bracken_abundances_",level,"_ALL")) %>% as.data.frame()
if(ind == "counts"){
  df0 = rawd[,c(1,seq(4,ncol(rawd),2))]
} else {
  df0 = rawd[,c(1,seq(5,ncol(rawd),2))]
}
s1 = paste0("_bracken_report_",level,".abundance_num")
sx0 = gsub(s1, "", colnames(df0))
sx1 = sx0[-1]
df1 = df0[,-1]
head(df0)
colnames(df1) = c(sx1)
df1 = df1 %>% relocate(c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"))
df1


# Metadata read 
mtd = fread("~/HPC_data/iron_reducing_bacteria/metadata.txt") # Use the metadata.txt if we want the original metadata and not layers 4 and 4B merged.
mtd$Origin = factor(mtd$Origin,levels = c("Spear","Soil"))
rownames(mtd) = mtd$ID
mtd$ID = factor(mtd$ID, levels=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"))

# Taxonomy read
taxonomy_tab = as.matrix(df0$name)
colnames(taxonomy_tab) = level
# assign rownames to the tax level chosen
rownames(taxonomy_tab) = taxonomy_tab[,1] # [,tax_lvl_int]
rownames(df1) = rownames(taxonomy_tab)


# CLEAN CROSSCHECK NAMES
inters_names = intersect(colnames(df1), rownames(mtd))
otu_cleaned = df1[, inters_names]
metadata_cleaned = mtd[match(inters_names, rownames(mtd)),]
mx1 = sample_data(metadata_cleaned, errorIfNULL = FALSE)
# mx1$caribou_dna = as.numeric(mx1$caribou_dna)
sample_names(mx1) = metadata_cleaned$ID

# Read into phyloseq
p0 = phyloseq(otu_table(as.matrix(otu_cleaned), taxa_are_rows = TRUE),
              tax_table(taxonomy_tab),mx1)
dim(otu_table(p0))

# View(otu_table(p0))
# FILTER SAMPLES #

# Samples N3 has very low counts
p1 = prune_samples(sample_sums(p0) >= 20000, p0)

# FILTER ON % ABUNDANCE (Using percent abundance!) ####
## p2  = transform_sample_counts(p1, function(x) x / sum(x) )
## ## Mean abundance >X
## p3 = filter_taxa(p2, function(x) mean(x) > 0.001, TRUE)
## dim(otu_table(p3))

# Alternatively, we can also filter using the total counts: (choose one or the other)

# FILTER USING TOTAL COUNTS (At least 10% taxa with more than X counts)
mins = 1000 # For all groupings
p3 = filter_taxa(p1, function(x) sum(x > mins) > (0.1*length(x)),TRUE)
dim(otu_table(p3))
# View(otu_table(p3))

# Do some barplots
p4 = subset_taxa(p3, G == "Desulfuromonas")  # Achromobacter
plot_bar(p4, fill="G")
p = plot_bar(p4, fill="G")
p$data$Sample = factor(p$data$Sample, levels=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11"))

# ALPHA DIVERSITY ####
# Correct by median sequencing depth
total = median(sample_sums(p1))
standf = function(x, t = total) round(t * (x / sum(x))) 
otu_table(p1)
pad = transform_sample_counts(p1, standf)
otu_table(pad)
ad = prune_taxa(taxa_names(pad) != "Homo", pad) # Remove human comtamination

# Plots
# Alpha div per date ####
theme_set(theme_bw())
palpha_origin = plot_richness(ad, x = "Origin", measures="Shannon", color = "Origin") +
  geom_point(size=4) +
  # geom_boxplot(alpha = 0.8, size=1)  +
  scale_color_manual(breaks = c("Spear","Soil"), values = c("darkgreen","lightblue")) +
  # scale_x_discrete(limits=c("",""),labels=c("1978","2021")) +
  stat_compare_means(label.x = 1.35, size = 8) + 
  xlab("Sample") +
  ggrepel::geom_text_repel(
    aes(label = ID),  # Label with Sample column
    size = 8,             # Text size
    box.padding = 0.5,    # Padding around text
    point.padding = 0.5,  # Padding around points
    force = 10,           # Repulsion force
    max.overlaps = 10     # Maximum number of overlapping labels
    ) + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 20)) +
  guides(color="none")

palpha_origin

img_dir = "/Users/jrodriguez/Projects/iron_reducing_bacteria/"
ggsave(filename = paste0(img_dir,"Shannon_alpha_origin_KrakenUniq.pdf"), palpha_origin)

# BETA-DIVERSITY ####
set.seed(1)

if(level == "S"){
  pb1 = prune_taxa(phyloseq::taxa_names(p1) != "Homo sapiens", p1)
} else if (level == "G") {
  pb1 = prune_taxa(phyloseq::taxa_names(p1) != "Homo", p1)
} else {stop("DEFINE LEVEL TO REMOVE HUMAN CONTAMINATION TAXA.")}

# Filter type
FILTER = "COUNTS" # "ABUND"; #  "ABUND"
dist.met = "bray";
meth1 = "NMDS" # "PCoA" ;
# meth1 = "PCoA" ;
min_abu = NULL

if(FILTER == "COUNTS"){
  #' Filter on counts and observation %
  #' "At least x counts in 10% of the samples"
  min_abu = 0.1
  pb2 = filter_taxa(pb1, function(x) sum(x > mins) > (min_abu*length(x)),TRUE)
  pb3 = transform_sample_counts(pb2, function(otu) otu/sum(otu))
  ord.nmds.bray = ordinate(pb3, method = meth1, distance = dist.met)
  } else if (FILTER == "ABUND") {

  #' Filter on average abundance %
  #' "At least x average abundance in samples"
  min_abu = 0.005
  pb.abun = transform_sample_counts(pb1, function(otu) otu/sum(otu))
  pb3 = filter_taxa(pb.abun, function(x) mean(x) > min_abu, TRUE)
  ord.nmds.bray = ordinate(pb3, method = meth1, distance = dist.met)
  # otu_table(pb3)
}
tit1 = paste0(meth1,"; ",level," level.\n",dist.met," | ",FILTER, " filtered: ", min_abu * 100,"%")

## NMDS Plot ####
View(otu_table(pb3))
plot_beta = plot_ordination(pb3,
                ord.nmds.bray,
                color = "Origin") +
  geom_point(size = 5) + 
  # ggtitle(tit1) +
  scale_color_manual(breaks = c("Spear","Soil"), values = c("darkgreen","lightblue")) +
  geom_text_repel(aes(label=ID), size=8) + 
  theme(text = element_text(size = 18),
        legend.title = element_blank())
# otu_table(pb3)

plot_beta

ggsave(paste0(img_dir,"NMDS_beta_origin_KrakenUniq.pdf"), plot_beta)

# Using aitchinson distance

plot_beta_aitch = pb2 %>% tax_transform(rank = "unique", trans = "normalize") %>% dist_calc(dist = "aitchison") %>%
  ord_calc(method = "auto") %>%  
  ord_plot(axes = c(1, 2),
           color="Origin",
           alpha = 0.5,
           size = 5) + 
  geom_text_repel(aes(label=ID), size=8) + 
  scale_color_manual(breaks = c("Spear","Soil"), values = c("darkgreen","lightblue")) +
  theme(text = element_text(size = 18))

plot_beta_aitch

ggsave(paste0(img_dir,"Aitch_beta_origin_KrakenUniq.pdf"), plot_beta_aitch)

# PERMANOVA ####
# Association of B-div with a metadata variable
# Calculate bray curtis distance matrix 
pb.dist = phyloseq::distance(pb3, method = "bray")
# make a data frame from the sample_data
pb.smp = data.frame(sample_data(pb2))
set.seed(1)
panv.date = adonis2(pb.dist ~ Origin, data = pb.smp, permutations = 999, method="bray")
panv.date
message("Permanova p-value for origin: ", panv.date$`Pr(>F)`[1])


# DIFFERENTIAL ABUNDANCE ####
# USING ALL OF THE SOIL SAMPLES
library("DESeq2")

## By Origin ####

das.date = phyloseq_to_deseq2(pb1, ~ Origin)
diagdd.date = DESeq(das.date) # , test="Wald", fitType="parametric")
res.date = results(diagdd.date) # , cooksCutoff = FALSE)

alpha = 0.05
logchange = 1.5

res.plot <- cbind(
  as(res.date, "data.frame"),
  as(tax_table(pb1)[rownames(res.date), ], "matrix")
) %>%
  as.data.frame()

res.plot$taxon_id <- rownames(res.plot)

res.plot <- res.plot %>%
  mutate(
    taxon_label = as.character(G),
    neglog10padj = -log10(padj),
    neglog10padj = ifelse(is.infinite(neglog10padj), NA, neglog10padj),
    pass_both = !is.na(padj) & padj < alpha & abs(log2FoldChange) >= logchange,
    colv = ifelse(pass_both, "coral1", "lightgrey")
  )

dif_abu <- ggplot(res.plot, aes(x = neglog10padj, y = log2FoldChange)) +
  geom_point(aes(color = colv), size = 2, alpha = 0.9) +
  geom_hline(yintercept = logchange, color = "red", lty = "dashed") +
  geom_hline(yintercept = -logchange, color = "red", lty = "dashed") +
  geom_vline(xintercept = -log10(alpha), color = "red", lty = "dashed") +
  scale_color_identity() +
  # scale_x_continuous(limits = c(0, 10), oob = squish) +
  geom_text_repel(data=res.plot[res.plot$colv=="coral1",],aes(label=G), min.segment.length = 0.1, max.overlaps = 14) +
  ggtitle(
    "Differential abundances soil vs. spears",
    subtitle = "-logFC: More abundant in spears"
  ) +
  xlab("-log10(padj)") +
  ylab("log2FoldChange") +
  theme_bw() +
  theme(text = element_text(size = 18)) + 
  coord_cartesian(xlim = c(0, 10))

dif_abu

ggsave(paste0(img_dir,"volcano_plot_spears.pdf"), dif_abu)

sig.date.f = res.plot %>% filter(abs(log2FoldChange) > logchange) %>% arrange((padj))
View(sig.date.f)
fwrite(sig.date.f, "~/Projects/iron_reducing_bacteria/Differential_Abundance_KrakenUniq_list.txt",sep=";",dec=",",quote=F,col.names = T,row.names = T)


# DO THE SAME ONLY WITH SAMPLE N9 FOR SOIL GROUP. ####
py = prune_samples(!(sample_names(pb1) %in% c("N7", "N8","N10","N11")), pb1)
px = prune_samples(sample_sums(py) >= 20000, py)
das.date = phyloseq_to_deseq2(px, ~ Origin)
diagdd.date = DESeq(das.date) # , test="Wald", fitType="parametric")
res.date = results(diagdd.date) # , cooksCutoff = FALSE)

alpha = 0.05
logchange = 1.5

target_taxa <- c("Rhodoferax","Ferribacterium","Desulfuromonas","Geobacter","Geothrix","Methylomonas","Shewanella")

res.plot <- cbind(
  as(res.date, "data.frame"),
  as(tax_table(pb1)[rownames(res.date), ], "matrix")
) %>%
  as.data.frame()

res.plot$taxon_id <- rownames(res.plot)

res.plot <- res.plot %>%
  mutate(
    taxon_label = as.character(G),
    neglog10padj = -log10(padj),
    neglog10padj = ifelse(is.infinite(neglog10padj), NA, neglog10padj),
    pass_both = !is.na(padj) & padj < alpha & abs(log2FoldChange) >= logchange,
    colv = ifelse(pass_both, "coral1", "lightgrey")
  )

# label all significant target taxa
#lab_targets <- res.plot %>% filter(pass_both, taxon_label %in% target_taxa)

# label additional significant non-target taxa
# choose the top N most extreme by significance and fold change
# n_extra <- 100
# lab_extra <- res.plot %>%
#   filter(pass_both, !(taxon_label %in% target_taxa)) %>%
#   mutate(rank_score = neglog10padj * abs(log2FoldChange)) %>%
#   arrange(desc(rank_score)) %>%
#   slice_head(n = n_extra)
# 
# lab_df <- bind_rows(lab_targets, lab_extra) %>%
#   distinct(taxon_label, .keep_all = TRUE)

dif_abu_N9 <- ggplot(res.plot, aes(x = neglog10padj, y = log2FoldChange)) +
  geom_point(aes(color = colv), size = 2, alpha = 0.9) +
  geom_hline(yintercept = logchange, color = "red", lty = "dashed") +
  geom_hline(yintercept = -logchange, color = "red", lty = "dashed") +
  geom_vline(xintercept = -log10(alpha), color = "red", lty = "dashed") +
  scale_color_identity() +
  # scale_x_continuous(limits = c(0, 10), oob = squish) +
  geom_text_repel(data=res.plot[res.plot$colv=="coral1",],aes(label=G),min.segment.length = 0.1,max.overlaps = 20) +
  ggtitle(
    "Differential abundances N9 sample vs. spears",
    subtitle = "-logFC: More abundant in spears"
  ) +
  xlab("-log10(padj)") +
  ylab("log2FoldChange") +
  theme_bw() +
  theme(text = element_text(size = 18)) + 
  coord_cartesian(xlim = c(0, 10))

dif_abu_N9

# res.plot[rownames(res.plot) == "Methylomonas",]
ggsave(paste0(img_dir,"volcano_plot_spears_N9.pdf"), dif_abu_N9)

sig.date.N9 = res.plot %>% filter(abs(log2FoldChange) > logchange) %>% arrange((padj))
View(sig.date.N9)
fwrite(sig.date.N9, "~/Projects/iron_reducing_bacteria/Differential_Abundance_KrakenUniq_list_ONLY_N9.txt",sep=";",dec=",",quote=F,col.names = T,row.names = T)

stop("DONE!")

############
############

# ORIGINAL VOLCANO PLOT WITH ALL OF THEM ####

# library("DESeq2")
# ## By Origin ####
# das.date = phyloseq_to_deseq2(pb1, ~ Origin)
# diagdd.date = DESeq(das.date) # , test="Wald", fitType="parametric")
# res.date = results(diagdd.date) # , cooksCutoff = FALSE)
# 
# # Format in DF
# alpha = 0.001
# sig.date0 = res.date[which(res.date$padj < alpha), ]
# sig.date1 = cbind(as(sig.date0, "data.frame"), as(tax_table(pb1)[rownames(sig.date0), ], "matrix"))
# sig.date1 %>% arrange(padj)
# 
# ### Plot pvalue x log-fold change ####
# logchange = 1.5
# sig.date = sig.date1 %>% mutate(colv = ifelse(abs(log2FoldChange) > logchange,"coral1" ,"lightgrey"))
# dif_abu = ggplot(data = sig.date, aes(x = -log10(padj), y = log2FoldChange, color = colv)) + 
#   geom_point() + 
#   geom_hline(yintercept =  logchange,color="red",lty = "dashed") +
#   geom_hline(yintercept = -logchange,color="red",lty = "dashed") +
#   scale_color_identity() +
#   geom_label_repel(data=sig.date[sig.date$colv=="coral1",],aes(label=G),min.segment.length = 0.1,max.overlaps = 20) +
#   ylim(-6,6) + 
#   ggtitle("Differential abundances soil vs. spears",subtitle="-logFC: More abundant in spears")
# 
# dif_abu
# 
# sig.date.f = sig.date %>% filter(abs(log2FoldChange) > logchange) %>% arrange((padj))
# dim(sig.date.f)
# fwrite(sig.date.f, "~/Projects/iron_reducing_bacteria/Differential_Abundance_KrakenUniq_list.txt",sep=";",dec=",",quote=F,col.names = T,row.names = T)
# ggsave(filename = paste0(img_dir,"/Differential_Abundance_KrakenUniq.pdf"),dif_abu)
# 
# 
# 
# # MY VERSION FOR N9: ####
# # res.date[rownames(res.date) == "Rhodoferax",]
# # Format in DF
# alpha = 0.05
# sig.date0 = res.date[which(res.date$padj < alpha), ]
# sig.date1 = cbind(as(sig.date0, "data.frame"), as(tax_table(pb1)[rownames(sig.date0), ], "matrix"))
# sig.date1 %>% arrange(padj)
# # View(res.date)
# 
# ### Plot pvalue x log-fold change ####
# logchange = 1.5
# sig.date = sig.date1 %>% mutate(colv = ifelse(abs(log2FoldChange) > logchange,"coral1" ,"lightgrey"))
# dif_abu_N9 = ggplot(data = sig.date, aes(x = -log10(padj), y = log2FoldChange, color = colv)) + 
#   geom_point() + 
#   geom_hline(yintercept =  logchange,color="red",lty = "dashed") +
#   geom_hline(yintercept = -logchange,color="red",lty = "dashed") +
#   scale_color_identity() +
#   geom_text_repel(data=sig.date[sig.date$colv=="coral1",],aes(label=G),min.segment.length = 0.1,max.overlaps = 20) +
#   # ylim(-6,6) + 
#   ggtitle("Differential abundances N9 sample vs. spears", subtitle="-logFC: More abundant in spears") + 
#   theme(text = element_text(size = 18)) + 
#   coord_cartesian(xlim = c(0, 10))
# dif_abu_N9
# 
# sig.date.f = sig.date %>% filter(abs(log2FoldChange) > logchange) %>% arrange((padj))
# dim(sig.date.f)


