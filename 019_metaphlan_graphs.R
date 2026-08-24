library(tidyverse)
library(readxl)
library(phyloseq)
library(ape)
library(data.table)
library(microbiome)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(pheatmap)
library(viridis)
# library(rstatix)

options(scipen = 10)

setwd("~/HPC_data/iron_reducing_bacteria/metaphlan/out")

# Function has been modified in the source code, obtained from: https://rdrr.io/github/g-antonello/gautils2/src/R/metaphlan_to_phyloseq.R#sym-metaphlan_to_phyloseq
# It was meant to read profile files directly.
source("~/HPC_data/SCRIPTS/caribou/_OLD/metaphlan_to_phyloseq.R") 

# Metadata read 
mtd = fread("../../metadata.txt")
mtd[] = lapply(mtd, function(x) gsub(",", ".", x))
rownames(mtd) = mtd$ID

# READ SATURATION ####
# Number of reads mapped
nreads = fread("read_counts.txt")
nreads$V1 = gsub("N", "", nreads$V1)
nreads$V1 = as.integer(nreads$V1)
nreads = nreads %>% arrange(V1)

# Number of total reads
tot_reads = fread("reads_sequenced.txt")
tot_reads$V1 = gsub("N", "", tot_reads$V1)
tot_reads$V1 = as.integer(tot_reads$V1)
tot_reads = tot_reads %>% arrange(V1)

reads_all = cbind(nreads,tot_reads)[,c(1,2,4)]
reads_all$tag = mtd$Origin
colnames(reads_all)=c("ID","map","tot","tag")

ggplot(data=reads_all,aes(x=tot,y=map,label=ID)) + 
  geom_point(size = 2.8, aes(color=tag)) +
  theme_minimal() +
  scale_x_continuous(breaks=c(0,50e6,100e6,150e6,200e6,250e6,300e6),labels = c("0","50", "100","150","200","250","300")) +
  xlab("Total reads (in millions)") +
  ylab("Metaphlan mapped reads") +
  # geom_smooth(method="loess",se=FALSE,aes(color=tag)) +
  ggrepel::geom_text_repel()

# READ METAPHLAN TABLE ####
# Read full metaphlan table
mpa_abu = fread("merged_ALL_abundance_table.txt") %>% as.data.frame()
sx1 = gsub("\\.profile", "", colnames(mpa_abu))
colnames(mpa_abu) = sx1
mpa_abu = mpa_abu %>% select("clade_name","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11")
# Move the unclassified row to the top
mpa_abu = mpa_abu %>% filter(clade_name == "UNCLASSIFIED") %>% bind_rows(mpa_abu %>% filter(clade_name != "UNCLASSIFIED"))
# Check same order for samples
# all.equal(mtd$ID, colnames(mpa_abu)[-1])

# Transform to absolute counts.
getCounts <- function(m = mpa_abu, rds = nreads){
  mpa = mpa_abu[,-1]
  samps = colnames(mpa)
  colnames(mpa) = samps
  mpa_counts = as.data.frame(t(apply(mpa,1 , function(x) round((x*nreads$V2)/100))))
  mpa_counts = bind_cols(mpa_abu[,1],mpa_counts)
  colnames(mpa_counts) = c("clade_name",samps)
  return(mpa_counts)
}
mpa_counts = getCounts()

# READ TO PHYLOSEQ ####
# theme_set(theme_classic())
input = mpa_counts
# input = mpa_abu
tax_lvl="Species"
physeq = metaphlan_to_phyloseq(mpa = input,
                               metadata = mtd,
                               tax_lvl = tax_lvl,
                               version = 4)
physeq
## PLOT RAW READS ####
# create ordered sample variable
sample_data(physeq)$ID <- factor(
  sample_names(physeq),
  levels = c("N1","N2","N3","N4","N5","N6",  # metal
             "N7","N8","N9","N10","N11")     # soil
)

plot_bar(physeq, x = "ID", fill = tax_lvl) +
  theme_bw()

## FILTER SAMPLES ####
# remove UNCLASSIFIED
p0 = subset_taxa(physeq, Class != "UNCLASSIFIED")
# remove samples with ~0 counts
p1 = subset_samples(p0, !ID %in% c("N3","N4"))

# NORMALIZE ####
# Remove taxa with average > X
# p2 = filter_taxa(p1, function(x) mean(x) > 1e-4, TRUE) 

# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
p2 = filter_taxa(p1, function(x) sum(x > 1000) > (0.1*length(x)), TRUE) 
dim(otu_table(p2))
otu_table(p2)
plot_bar(p2, x = "ID", fill = tax_lvl) +
  theme_bw()

### TRANSFORM COUNTS ####
# TRANSFORM abundances to the median sequencing depth
# total = median(sample_sums(p2))
# standf = function(x, t=total) round(t * (x / sum(x))) 
# p3 = transform_sample_counts(p2, standf)
# dim(otu_table(p3))
# otu_table(p3)

# ALTERNATIVELY, TRANSFORM TO TOTAL COL SUM = 1
p3 = transform_sample_counts(p2, function(x) x / sum(x) ) 
colSums(otu_table(p3))

# Filter taxa for low coefficient of variation
# p4 = filter_taxa(p3, function(x) sd(x)/mean(x) >  1, TRUE)
# dim(otu_table(p4))

# HEATMAP ABUNDANCE ####
# colx = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
colx = viridis(30)

# metadata
### add more variables... layer?
my_sample_col = as.data.frame(mtd$Origin)
colnames(my_sample_col) = c("Origin")
# my_sample_col$layer = gsub("Layer " , "" , mtd$layer)
# colnames(my_sample_col) = c("Date","Layer")
row.names(my_sample_col) = mtd$ID
ann_colors = list(Origin = c("Spear"="darkgreen", "Soil"="lightblue"))

# x1 = otu_table(physeq)
x1 = otu_table(p3)
htm = pheatmap(x1, scale="row",
         color = colx, 
         annotation_col = my_sample_col,
         annotation_colors = ann_colors[1],
         border_color = NA)
htm
img_dir = "/Users/jrodriguez/Projects/iron_reducing_bacteria/"
ggsave(paste0(img_dir,"heatmap_metaphlan.pdf"), htm)
  
  
  
