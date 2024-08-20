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

path <- "/projects/marralab/cayan_prj/PrecisionMed/Objects/"
files <- list.files(path)

wholeCohort <- fread("/projects/marralab/cayan_prj/PrecisionMed/Data/WholeCohort.csv")