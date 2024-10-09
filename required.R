# data processing
library(tidyverse)
library(data.table)
library(Matrix)
library(Seurat)
library(Signac)

# plotting
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(jcolors)

# r markdown
library(xaringanExtra)

path <- "/projects/marralab/cayan_prj/PrecisionMed/Objects/"
files <- list.files(path)

wholeCohort <- fread("/projects/marralab/cayan_prj/PrecisionMed/Data/WholeCohort.csv")

atac_barcodes <- fread("/projects/marralab/Single_Cell_Projects/CellRanger/Software/CellRanger_ARC/cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz",
                       header = F) %>%
  mutate(V1 = paste(V1, "1", sep = "-"))
rna_barcodes <- fread("/projects/marralab/Single_Cell_Projects/CellRanger/Software/CellRanger_ARC/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz",
                      header = F) %>%
  mutate(V1 = paste(V1, "1", sep = "-"))

# functions
# converting between barcodes of different modalities
convert_rna_indices <- function(rna_bc){
  rna_indices <- match(rna_bc, rna_barcodes$V1)
  valid_atac_bc <- atac_barcodes$V1[rna_indices]
  return(valid_atac_bc)
}
convert_atac_indices <- function(atac_bc){
  atac_indices <- match(atac_bc, atac_barcodes$V1)
  valid_rna_bc <- rna_barcodes$V1[atac_indices]
  return(valid_rna_bc)
}

plotQC <- function(f){
  qc <- fread(paste0(path, f, "/QCMetrics_byCell.tsv"))
  logLib <- log10(qc$sum)
  logGenes <- log10(qc$detected)
  
  libLine <- median(logLib) - 3 * mad(logLib)
  genesLine <- median(logGenes) - 3 * mad(logGenes)
  qc[is.na(qc)] <- "no assignment"
  
  p1 <- ggplot(qc, aes(x = log10(sum), y = log10(detected), colour = doublet)) +
    geom_point(size = 1, alpha = 0.25) +
    geom_hline(yintercept = genesLine, linetype="dashed", colour = "grey") +
    geom_vline(xintercept = libLine, linetype="dashed", colour = "grey") +
    theme_linedraw() +
    scale_color_d3() +
    labs(x = "Log10(Number of UMIs)",
         y = "Log10(Number of Genes)",
         colour = "scDblFinder\nPrediction")
  
  p2 <- ggplot(qc, aes(x = log10(sum), y = subsets_Mito_percent, colour = discard)) +
    geom_point(size = 1, alpha = 0.25) +
    # geom_hline(yintercept = genesLine, linetype="dashed", colour = "grey") +
    geom_vline(xintercept = libLine, linetype="dashed", colour = "grey") +
    theme_linedraw() +
    scale_color_d3() +
    labs(x = "Log10(Number of UMIs)",
         y = "Percent of Mitochondrial Reads",
         colour = "Discard")
  return(plot_grid(p1, p2))
}


