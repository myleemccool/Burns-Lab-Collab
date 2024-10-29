#### Assigning taxonomy ####
#create a log file
logfile <- "taxonomy_pipline_log.txt"

#the sink() function allows you to incorperate logging
sink(logfile, append = TRUE, split = TRUE) #send output to loog file

#log start time and session info -- maybe not nessicary but could be useful
cat("### Taxonomy Pipeline Log ###\n")
cat("Start time: ", Sys.time(), "\n")
sessionInfo()
cat("\n")

cat("### Step 1: Downloading Taxonomy Databases ###\n")

# Download SILVA training set
cat("Downloading SILVA training set...\n")
download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz", 
              destfile = "silva_nr99_v138.1_train_set.fa.gz")
cat("SILVA training set downloaded successfully.\n")

# Download RDP training set
cat("Downloading RDP training set...\n")
download.file("https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz", 
              destfile = "rdp_train_set_18.fa.gz")
cat("RDP training set downloaded successfully.\n")

# Download Greengenes training set
cat("Downloading Greengenes training set...\n")
download.file("https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz",
              destfile = "gg_13_8_train_set_97.fa.gz")
cat("Greengenes training set downloaded successfully.\n")

cat("\n### Step 2: Running Taxonomic Classification with SILVA ###\n")

# Run taxonomic classification with SILVA
cat("Running assignTaxonomy for SILVA...\n")
silva_taxa <- assignTaxonomy(seqtab, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
cat("SILVA taxonomic classification completed. Total number of taxa classified: ", sum(!is.na(silva_taxa[, 1])), "\n")

# Add species-level assignments
cat("Adding species-level assignments using SILVA...\n")
silva_taxa <- addSpecies(silva_taxa, "silva_species_assignment_v138.1.fa.gz")
cat("Species-level assignment added.\n")

# Identify sequences that were not classified at the domain level
unclassified_silva <- is.na(silva_taxa[, 1])
cat("Number of sequences not classified by SILVA: ", sum(unclassified_silva), "\n")

# Step 3: Running taxonomic classification with RDP
if (any(unclassified_silva)) {
  cat("\n### Step 3: Running Taxonomic Classification with RDP ###\n")
  rdp_taxa <- assignTaxonomy(seqtab[, unclassified_silva], "rdp_train_set_18.fa.gz", multithread = TRUE)
  cat("RDP taxonomic classification completed. Total number of taxa classified: ", sum(!is.na(rdp_taxa[, 1])), "\n")
  
  # Merge results
  silva_taxa[unclassified_silva, ] <- rdp_taxa
  cat("Merged SILVA and RDP classification results.\n")
}

# Identify sequences that were still not classified
unclassified_rdp <- is.na(silva_taxa[, 1])
cat("Number of sequences still not classified after RDP: ", sum(unclassified_rdp), "\n")

# Step 4: Running taxonomic classification with Greengenes
if (any(unclassified_rdp)) {
  cat("\n### Step 4: Running Taxonomic Classification with Greengenes ###\n")
  gg_taxa <- assignTaxonomy(seqtab[, unclassified_rdp], "gg_13_8_train_set_97.fa.gz", multithread = TRUE)
  cat("Greengenes taxonomic classification completed. Total number of taxa classified: ", sum(!is.na(gg_taxa[, 1])), "\n")
  
  # Merge results
  silva_taxa[unclassified_rdp, ] <- gg_taxa
  cat("Merged SILVA, RDP, and Greengenes classification results.\n")
}

# Resulting taxa table
final_taxa <- silva_taxa
cat("\n### Final Taxonomy Table ###\n")
cat("Number of taxa in final table: ", nrow(final_taxa), "\n")

# Save the final taxonomy table
write.csv(final_taxa, "combined_taxonomy_results.csv", row.names = TRUE)
cat("Final taxonomy table saved as 'combined_taxonomy_results.csv'.\n")

# Inspect taxonomic assignments
taxa.print <- final_taxa
rownames(taxa.print) <- NULL
cat("Preview of final taxonomy assignments:\n")
print(head(taxa.print))

# End of log
cat("\nEnd time: ", Sys.time(), "\n")
sink()  # Stop logging

#### Now lets creat a phyloseq object :) ####
#first start off by loading the phyloseq object using the library function
#library() loads function into enviorment 
#packageVersion retreives currently installed version of the package you are lodaing
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
#ggplot is useed to create visulizations and complex plots 
library(ggplot2); packageVersion("ggplot2")
# this is optional but you can adjust the theme: lots of different fun themes to use
theme_set(theme_bw())
#This part is also not super important, it is essentially creating a data table using fuctions in r 
#This extracts the names from the objects 
samples.out <- rownames(seqtab.nochim)
#for this data this fuction splits sample IDs from the subject sample names
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
# creates a gender variable
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
#splits samples out at D and converts to integers stored in day
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
#created new columns labled subject, gender, and day
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
#assigns "early" and "late" to new created column titled when
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
#this is the main code needed to construct a phyloseq object 
#seqtab.nochim is the OTU (operational taxonomic units) table
#taxa_are_rows=FALSE means the OTUs are in columns (samples are in rows).
#sample_data(samdf): Adds samdf as the sample metadata for the phyloseq object - change the samdf according to your own data
#tax_table(taxa): Adds taxa as the taxonomy table, containing taxonomic information (e.g., phylum, class) for the OTUs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
#sample_names(ps) != "Mock" creates a logical vector where TRUE means a sample's name is not "Mock".
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
#After running this code, ps will be a phyloseq object that excludes the sample labeled "Mock."