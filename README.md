# microbial-corrosion-layers
Scripts used to analyze shotgun metagenome data from corrosion layers related to the paper under review:

**Title:** "_Distinct microbial communities associated with corrosion layer transformation in reburied Archaeological iron in an anoxic peat environment_" 

**Authors**: Frydendhal _et al_., 2026

**Brief explanation of the code provided to analyse DNA data:**

(_All code in Bash or R._)

**Preprocessing data**; 

_filtering 001_fastp_

**Collapse/Merge reads**;

_008_collapse_

_009_merge_collapsed+forward_reads_

**KrakenUniq classification and abundance analyses**;

_010_krakenuniq_

_014_bracken_

_015_bracken_analyses (plots)_

_016_barplot_abundance_target_species (plots)_

**Metaphlan**

_017_metaphlan_

_018_metaphlan_read_counts_

_019_metaphlan_graphs (plots)_
