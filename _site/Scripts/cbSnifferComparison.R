# Purpose: 
# To run cbSniffer on: 
# - all cells with high-quality RNA
# - cells with both high-quality RNA and ATAC
# - scRNA-seq BAM files

source("../required.R")
args_file <- fread("/projects/marralab/cayan_prj/PrecisionMed/Data/WholeCohort.csv")

for(i in seq(1, nrow(args_file))[-10]){
  f <- args_file$folder_name[i]
  
  file_path_atac <- args_file$file_path_atac[i]
  atac_bam_path <- gsub("filtered_peak_bc_matrix", "possorted_bam.bam", file_path_atac)
  
  file_path_rna <- args_file$file_path[i]
  rna_bam_path <- gsub("filtered_feature_bc_matrix", "possorted_genome_bam.bam", file_path_rna)
  
  if(!dir.exists(paste0(path, f, "/cbSnifferOuts/"))){
    dir.create(paste0(path, f, "/cbSnifferOuts/"))
  }
  cbDir <- paste0(path, f, "/cbSnifferOuts/")
  
  # cell barcodes for all high-quality RNA cells
  qcFile <- fread(paste0(path, f, "/QCMetrics_byCell.tsv")) %>%
    dplyr::filter(discard == FALSE)
  
  rnaBarcodes <- qcFile$barcode
  data.frame(barcodes = rnaBarcodes) %>%
    write_tsv(paste0(cbDir, "all_rnaBarcodes.tsv"))
  
  atacBarcodes <- convert_rna_indices(rnaBarcodes)
  data.frame(barcodes = atacBarcodes) %>%
    write_tsv(paste0(cbDir, "all_atacBarcodes.tsv"))
  
  # # cell barcodes for cells with high-quality RNA and ATAC
  # bcFile <- fread(paste0(path, f, "/barcodesClusters_labelled.csv"))
  # hq_atacBarcodes <- convert_rna_indices(bcFile$barcode)
  # data.frame(barcodes = hq_atacBarcodes) %>%
  #   write_tsv(paste0(cbDir, "hq_atacBarcodes.tsv"))
  
  atac_cb_string <- paste("python3 ../../sander_mutationcalling_cb_sniffer.py",
                          atac_bam_path,
                          paste0(path, f, "/self_POGSNVcalls.tsv"),
                          paste0(cbDir, "all_atacBarcodes.tsv"),
                          paste0(cbDir, "allCells_cbATAC"),
                          # "-mq 10",
                          "-bq 30") # !! for separate modalities, need to convert barcodes back to ATAC (or need to store in metadata)
  system(atac_cb_string)
  rna_cb_string <- paste("python3 ../../sander_mutationcalling_cb_sniffer.py",
                         rna_bam_path,
                         paste0(path, f, "/self_POGSNVcalls.tsv"),
                         paste0(cbDir, "all_rnaBarcodes.tsv"),
                         paste0(cbDir, "allCells_cbRNA"),
                         # "-mq 10",
                         "-bq 30")
  system(rna_cb_string)
}


