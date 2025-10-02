# ============================================================
# Load required packages
# ============================================================
install.packages(c("tidyverse","data.table","readxl","BiocManager"))
BiocManager::install(c("dada2","phyloseq","DECIPHER","Biostrings","phangorn"))

library(dada2)
library(phyloseq)
library(readxl)
library(tidyverse)
library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# ============================================================
# Paths (EDIT if needed)
# ============================================================
# Your FASTQ folder from the screenshot
fastq_dir <- "/Users/delaneyrager/Desktop/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/Individual_Seqs"

# Metadata Excel file (put the correct filename here)
metadata_file <- "/Users/delaneyrager/Desktop/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/KitTestingMetadata.xlsx"

# SILVA database files (make sure these exist in your system)
silva_train   <- "/Users/delaneyrager/Desktop/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/silva_nr99_v138_train_set.fa.gz"
silva_species <- "/Users/delaneyrager/Desktop/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/silva_species_assignment_v138.fa.gz"

# ============================================================
# Step 1. List FASTQs and extract sample names
# ============================================================
fnFs <- sort(list.files(fastq_dir, pattern="_R1_001\\.fastq(\\.gz)?$", full.names=TRUE))
fnRs <- sort(list.files(fastq_dir, pattern="_R2_001\\.fastq(\\.gz)?$", full.names=TRUE))

# Extract sample names (e.g. MB001B_S150 from MB001B_S150_L001_R1_001.fastq.gz)
get_sample <- function(x) sub("_L001_.*", "", basename(x))
sample.names <- get_sample(fnFs)

# ============================================================
# Step 2. Filter and trim
# ============================================================
filt_dir <- file.path(fastq_dir, "filtered"); dir.create(filt_dir, showWarnings=FALSE)

filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240,200), maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)

print(out)

# Drop any missing filtered files
keep <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
sample.names <- sample.names[keep]

# ============================================================
# Step 3. Learn errors and run DADA2
# ============================================================
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

cat("Total reads kept: ",
    sum(seqtab.nochim), " of ", sum(seqtab), " (",
    round(100*sum(seqtab.nochim)/sum(seqtab),1), "%)\n", sep="")

# ============================================================
# Step 4. Assign taxonomy
# ============================================================
tax <- assignTaxonomy(seqtab.nochim, silva_train, multithread=TRUE)
tax <- addSpecies(tax, silva_species)

# ============================================================
# Step 5. Load metadata (Excel)
# ============================================================
meta <- read_excel(metadata_file)
meta <- as.data.frame(meta)

# IMPORTANT: make sure your Excel sheet has a column "SampleID"
# matching the names in sample.names
rownames(meta) <- meta$SampleID

# ============================================================
# Step 6. Build phyloseq object
# ============================================================
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  tax_table(as.matrix(tax)),
  sample_data(meta)
)

ps

# Save for later use
saveRDS(ps, file=file.path(fastq_dir, "phyloseq_object.rds"))
###########################################
# DADA2 + phyloseq pipeline (corrected)
###########################################

# --- Load libraries ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2")
}
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())

# --- Define path ---
path <- "~/Desktop/MiSeq_SOP"  # adjust if needed
list.files(path)

# --- Identify fastq files ---
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# --- Quality profiles ---
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# --- Filter and trim ---
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# --- Learn errors ---
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# --- Denoising ---
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# --- Merge pairs ---
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# --- Make sequence table ---
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# --- Track reads through pipeline ---
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# --- Taxonomy assignment ---
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/MiSeq_SOP/silva_nr_v132_train_set.fa", multithread=TRUE)

# --- Build sample metadata ---
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- ifelse(samdf$Day > 100, "Late", "Early")
rownames(samdf) <- samples.out

# -----------------------------
# Create phyloseq objects
# -----------------------------
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(taxa)
)

# Mock-only phyloseq object (for QC)
ps.mock <- prune_samples(sample_names(ps) == "Mock", ps)

# Biological samples only (for downstream analysis)
ps.analysis <- prune_samples(sample_names(ps) != "Mock", ps)

# Add DNA sequences and rename taxa in analysis object
dna <- Biostrings::DNAStringSet(taxa_names(ps.analysis))
names(dna) <- taxa_names(ps.analysis)
ps.analysis <- merge_phyloseq(ps.analysis, dna)
taxa_names(ps.analysis) <- paste0("ASV", seq(ntaxa(ps.analysis)))

# Inspect both
ps.mock
ps.analysis

# -----------------------------
# Example downstream analysis
# -----------------------------
plot_richness(ps.analysis, x="Day", measures=c("Shannon", "Simpson"), color="When")

ps.prop <- transform_sample_counts(ps.analysis, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps.analysis), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps.prop)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

# Optional: rarefaction
rarefied_ps <- rarefy_even_depth(ps.analysis, rngseed=123, verbose=TRUE)
rarefied_ps
