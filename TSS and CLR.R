####transforming using CLR and TSS in the phyloseq object####


# Load necessary libraries
library(phyloseq)
library(microbiome) # for CLR transformation

# Assuming your phyloseq object is named `ps`

# 1. TSS Transformation
ps_tss <- transform_sample_counts(ps, function(x) x / sum(x))

# 2. CLR Transformation
ps_clr <- microbiome::transform(ps_tss, "clr")

# Check the transformations
ps_tss
ps_clr