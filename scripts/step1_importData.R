# Import libraries----
library(tidyverse)
library(EnsDb.Mmusculus.v79)
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(gridExtra)

# Set  export directories----
figure_directory <- "results/figures"
table_directory <- "results/tables"
versions_folder <- "results/versions"

# Functions----
## Save plot
save_plot <- function(filename, plot, width=6, height=4, dpi=300) {
  # Create a dated subfolder in the versions folder
  date_subfolder <- format(Sys.Date(), "%Y%m%d")  # Current date in YYYYMMDD format
  version_directory <- file.path(versions_folder, date_subfolder)
  if (!dir.exists(version_directory)) dir.create(version_directory)
  
  # Save the plot in the figures and versions folders
  ggsave(
    filename,
    plot = plot,
    device = NULL,
    path = figure_directory,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = TRUE,
  )
  ggsave(
    filename,
    plot = plot,
    device = NULL,
    path = version_directory,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = TRUE,
  )
  
  # Return the paths of the saved files
  return(list(main_path = figure_directory, versioned_path = version_directory))
  
}

## Save table
save_table <- function (filename, table, quote=FALSE, row.names=TRUE, col.names=TRUE) {
  
  # Create a dated subfolder in the versions folder
  date_subfolder <- format(Sys.Date(), "%Y%m%d")  # Current date in YYYYMMDD format
  version_directory <- file.path(versions_folder, date_subfolder)
  if (!dir.exists(version_directory)) dir.create(version_directory)
  
  path_table = file.path(table_directory, filename)
  write.table(table, path_table, sep = "\t", quote = quote, 
              row.names = row.names, col.names = col.names)
  
  path_version = file.path(version_directory, filename)
  write.table(table, path_version, sep = "\t", quote = quote, 
              row.names = row.names, col.names = col.names)
  
  # Return the paths of the saved files
  return(list(main_path = path_table, versioned_path = path_version))
  
}

# Make DESeq2 DataSet----
path <- "data/counts"

targets <- na.omit(read_csv("241019_STXBP1_BulkRNAseq_librarypreparation.csv"))
sampleFiles <- grep("*",list.files(path),value=TRUE)
sampleTable <- data.frame(sampleName = targets$experiment_id, 
                          fileName = sampleFiles,
                          genotype = targets$genotype,
                          region = targets$region,
                          litter = targets$litter,
                          #batch = targets$batch,
                          timepoint = targets$timepoint,
                          sex = targets$sex)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = path,
                                       design= ~factor(genotype) + factor(region)
                                  # + factor(litter) + factor(batch) 
                                  # + factor(timepoint) + factor(sex)
)

dds

smallestGroupSize <- 1 # Change this to the smallest group size in your dataset
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
