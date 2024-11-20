# Import libraries----

## Make gene annotations data frame----
dds

genes <- genes(EnsDb.Mmusculus.v79, columns=c("gene_id", "gene_name"))
genes <- as_tibble(genes)
genes <- dplyr::select(genes, "gene_id", "gene_name")
# Stxbp1 levels
genes_name_dds <- genes$gene_name[match(rownames(dds), genes$gene_id)]
dds_annotated <- dds
rownames(dds_annotated) <- ifelse(!is.na(genes_name_dds), genes_name_dds, rownames(dds_annotated))
dds_annotated

d <- plotCounts(dds_annotated, gene="Stxbp1", intgroup=c("genotype", "region"), returnData=TRUE, normalized = TRUE,replaced = TRUE)

d <- d %>%
  mutate(Group = paste0(genotype,region))
# Plotting the MOV10 normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = Group, y = count, shape = genotype, color=region)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0), size = 2) +
  theme_bw() +
  ggtitle("Stxbp1") +
  theme(plot.title = element_text(hjust = 0.5))
