
#may need to download this from dada2
taxa <- assignTaxonomy(seqtab.nochim, "/Users/myleemccool/desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#dad2 does not like the newest silva verison - will just use the old one for now

#### downloading all taxonomy training sets #### 
#download silva training set (using this code)
download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz", 
              destfile = "silva_nr99_v138.1_train_set.fa.gz")

# install RDP taining set from dada2 website in taxonomy section (this code should work)
download.file("https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz", destfile = "rdp_train_set_18.fa.gz")

#install greengenes training set 
download.file("https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz",
              destfile = "gg_13_8_train_set_97.fa.gz")

#### 
# Step 1: Run taxonomic classification with SILVA
silva_taxa <- assignTaxonomy(seqtab, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
silva_taxa <- addSpecies(silva_taxa, "silva_species_assignment_v138.1.fa.gz")  # Add species-level assignments if available

# Identify sequences that were not classified
unclassified_silva <- is.na(silva_taxa[, 1])  # Checking which sequences didn't get classified at the Domain level

# Step 2: Run taxonomic classification with RDP on unmatched sequences
if (any(unclassified_silva)) {
  rdp_taxa <- assignTaxonomy(seqtab[, unclassified_silva], "rdp_train_set_18.fa.gz", multithread = TRUE)
  # Merge results, keeping SILVA classifications where available
  silva_taxa[unclassified_silva, ] <- rdp_taxa
}

# Identify sequences that were still not classified
unclassified_rdp <- is.na(silva_taxa[, 1])

# Step 3: Run taxonomic classification with Greengenes on unmatched sequences
if (any(unclassified_rdp)) {
  gg_taxa <- assignTaxonomy(seqtab[, unclassified_rdp], "gg_13_8_train_set_97.fa.gz", multithread = TRUE)
  # Merge results, keeping SILVA and RDP classifications where available
  silva_taxa[unclassified_rdp, ] <- gg_taxa
}

# Resulting taxa table
final_taxa <- silva_taxa

# View the first few rows
head(final_taxa)

# Optionally, save the final taxonomy table to a CSV file
write.csv(final_taxa, "combined_taxonomy_results.csv", row.names = TRUE)
#inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
