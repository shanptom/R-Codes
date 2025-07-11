phyloseq_to_excel <- function(ps_object, file_name = "phyloseq_export.xlsx") {
  # Load required packages
  if (!requireNamespace("phyloseq", quietly = TRUE) || !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Required packages 'phyloseq' and 'openxlsx' are not installed.")
  }
  
  # Extract and orient OTU table (taxa as rows)
  otu_mat <- phyloseq::otu_table(ps_object)
  if (!phyloseq::taxa_are_rows(otu_mat)) {
    otu_mat <- t(otu_mat)
  }
  otu_table_df <- as.data.frame(otu_mat)
  
  # Extract taxonomy table
  tax_table_df <- as.data.frame(phyloseq::tax_table(ps_object))
  
  # Extract reference sequences (if available)
  if (!is.null(phyloseq::refseq(ps_object))) {
    ref_seqs <- as.data.frame(phyloseq::refseq(ps_object))
    colnames(ref_seqs) <- "Reference_Sequences"
  } else {
    ref_seqs <- data.frame(Reference_Sequences = "No reference sequences available")
  }
  
  # Create list of data frames
  data_list <- list(
    OTU_Table = otu_table_df,
    Taxonomy_Table = tax_table_df,
    Reference_Sequences = ref_seqs
  )
  
  # Write to Excel
  openxlsx::write.xlsx(data_list, file = file_name, rowNames = TRUE)
}




tables_to_phyloseq <- function(
    count_file,
    taxonomy_file,
    metadata_file,
    refseq_file = NULL
) {
  # Helper function to read CSV or TSV based on file extension
  read_table <- function(file) {
    if (grepl("\\.csv$", file, ignore.case = TRUE)) {
      read.csv(file, row.names = 1, check.names = FALSE)
    } else {
      read.delim(file, row.names = 1, check.names = FALSE)
    }
  }
  
  # Load required packages
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Please install the 'phyloseq' package.")
  if (!is.null(refseq_file) && !requireNamespace("Biostrings", quietly = TRUE)) stop("Please install the 'Biostrings' package to handle reference sequences.")
  
  # Load count table
  count_df <- read_table(count_file)
  otu <- phyloseq::otu_table(as.matrix(count_df), taxa_are_rows = TRUE)
  
  # Load taxonomy table
  tax_df <- read_table(taxonomy_file)
  tax <- phyloseq::tax_table(as.matrix(tax_df))
  
  # Load metadata
  meta_df <- read_table(metadata_file)
  meta <- phyloseq::sample_data(meta_df)
  
  # Create phyloseq object
  ps <- phyloseq::merge_phyloseq(otu, tax, meta)
  
  # Load optional reference sequences
  if (!is.null(refseq_file)) {
    refseq_df <- read_table(refseq_file)
    if (ncol(refseq_df) != 1) stop("Reference file must have only one column of sequences.")
    ref_seqs <- Biostrings::DNAStringSet(refseq_df[[1]])
    names(ref_seqs) <- rownames(refseq_df)
    ps <- phyloseq::merge_phyloseq(ps, phyloseq::refseq(ref_seqs))
  }
  
  return(ps)
}



build_phylogenetic_tree <- function(ps_object, 
                                    prefix = "output",
                                    k = 4, 
                                    inv = 0.2,
                                    optimize_tree = FALSE) {
  # Load required libraries
  if (!requireNamespace("DECIPHER", quietly = TRUE) ||
      !requireNamespace("phangorn", quietly = TRUE) ||
      !requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Please install the required packages: DECIPHER, phangorn, phyloseq")
  }
  
  # Extract reference sequences
  ASVs <- ps_object@refseq
  if (is.null(ASVs)) stop("No reference sequences found in the phyloseq object.")
  
  # Step 1: Align sequences
  alignment <- DECIPHER::AlignSeqs(ASVs, anchor = NA)
  saveRDS(alignment, paste0(prefix, "_alignment.rds"))
  
  # Step 2: Convert alignment to phyDat format
  phang_align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  saveRDS(phang_align, paste0(prefix, "_phang.align.rds"))
  
  # Step 3: Compute distance matrix
  dm <- phangorn::dist.ml(phang_align)
  saveRDS(dm, paste0(prefix, "_dm.rds"))
  
  # Step 4: Construct NJ tree
  treeNJ <- phangorn::NJ(dm)
  saveRDS(treeNJ, paste0(prefix, "_treeNJ.rds"))
  
  # Step 5: Fit initial tree
  fit <- phangorn::pml(treeNJ, data = phang_align)
  fitGTR <- phangorn::update(fit, k = k, inv = inv)
  
  # Optional: Optimize tree
  if (optimize_tree) {
    fitGTR <- phangorn::optim.pml(
      fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
      rearrangement = "stochastic", control = phangorn::pml.control(trace = 0)
    )
  }
  
  saveRDS(fitGTR, paste0(prefix, "_fitGTR.rds"))
  
  # Step 6: Merge tree into phyloseq object
  ps_with_tree <- phyloseq::merge_phyloseq(ps_object, phyloseq::phy_tree(fitGTR$tree))
  saveRDS(ps_with_tree, paste0(prefix, "_tree.rds"))
  
  return(ps_with_tree)
}


