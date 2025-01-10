# Import libraries----
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggplotify)

# Transformation, sample clustering and visualization----
## Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
# vsd uses median of ratios method to normalize the data
# Adjust raw RNA-seq to account for differences in sequencing depth and RNA composition
# Size factors are calculated to normalize differences between samples, not within a single sample.
# For comparing expression of genes within the same sample (e.g. Gene A vs Gene B), use RPKM, TPM, TMM.
vsd_plot <- meanSdPlot(assay(vsd), plot=TRUE)
#The plot should appear flat, sd should not increase with mean
vsd_ggplot <- vsd_plot$gg + 
  labs(caption = paste0("produced on ", format(Sys.time(), "%x %X")))

save_plot("meanSdPlot.pdf", vsd_ggplot)

## Perform median of ratios method of normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)
save_table("normalized_counts.txt", normalized_counts)

## Heatmaps of count matrix

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
### Heatmap genotype, region and timepoint
df <- as.data.frame(colData(dds)[,c("genotype", "region", "timepoint")])
heatmap_obj <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                        cluster_cols=TRUE, annotation_col=df)
# We should see a clear separation of the samples by genotype and region
heatmap_gg <- as.ggplot(heatmap_obj$gtable)
heatmap_gg <- heatmap_gg + labs(caption = paste0("produced on ", format(Sys.time(), "%x %X")))
save_plot("heatmapCountMatrix_GenotypeVsRegionVsTimepoint.pdf", heatmap_gg,
          width = 7, height = 7)

## sample-sample distances
sampleDists <- dist(t(assay(vsd)))
# Sample-to-sample distances

### Heatmap genotype and region
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, vsd$region, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap_sample <- pheatmap(sampleDistMatrix,
                           clustering_distance_rows=sampleDists,
                           clustering_distance_cols=sampleDists,
                           col=colors)

heatmap_gg <- as.ggplot(heatmap_sample$gtable)
heatmap_gg <- heatmap_gg + labs(caption = paste0("produced on ", format(Sys.time(), "%x %X")))
save_plot("heatmapSample_GenotypeVsRegion.pdf", heatmap_gg)

### Heatmap genotype and timepoint
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, vsd$timepoint, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap_sample <- pheatmap(sampleDistMatrix,
                           clustering_distance_rows=sampleDists,
                           clustering_distance_cols=sampleDists,
                           col=colors)

heatmap_gg <- as.ggplot(heatmap_sample$gtable)
heatmap_gg <- heatmap_gg + labs(caption = paste0("produced on ", format(Sys.time(), "%x %X")))
save_plot("heatmapSample_GenotypeVsTimepoint.pdf", heatmap_gg)
# We should see a clear separation of the samples by genotype and region

