if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.19")
BiocManager::install("phyloseq")

install.packages("devtools")
devtools::install_github("ChiLiubio/microeco")
devtools::install_github('schuyler-smith/phylosmith')
devtools::install_github("gauravsk/ranacapa")
devtools::install_github("zdk123/SpiecEasi")

if(!require("file2meco")) install.packages("file2meco", repos = BiocManager::repositories())



