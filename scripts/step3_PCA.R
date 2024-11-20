# Import libraries----
library(DESeq2)
library(ggplot2)

# Principal Components Analysis----
# assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch = vsd$Batch, group = vsd$group)

## PCA genotype and region
pcaData <- plotPCA(vsd,intgroup=c("genotype", "region"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = region)) +
  geom_point(size=3) +
  coord_fixed()  + 
  labs(title = "PCA: Genotype ~ Region \nVST",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance"),
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) + 
  theme_bw()

p1

## PCA genotype and region
pcaData <- plotPCA(vsd,intgroup=c("genotype", "timepoint"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, shape = genotype, color =timepoint)) +
  geom_point(size=3) +
  coord_fixed()  +
  labs(title = "PCA: Genotype ~ Timepoint \nVST",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance"),
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) + 
  theme_bw()
p2

## PCA timepoint and region
pcaData <- plotPCA(vsd,intgroup=c("timepoint", "region"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p3 <- ggplot(pcaData, aes(PC1, PC2, shape = region, color =timepoint)) +
  geom_point(size=3) +
  coord_fixed()  +
  labs(title = "PCA: Genotype ~ Timepoint \nVST",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance"),
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) + 
  theme_bw()
p3

## PCA sex and litter
pcaData <- plotPCA(vsd,intgroup=c("sex", "litter"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p4 <- ggplot(pcaData, aes(PC1, PC2, shape = sex, color =litter)) +
  geom_point(size=3) +
  coord_fixed()  +
  labs(title = "PCA: Genotype ~ Timepoint \nVST",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance"),
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) + 
  theme_bw()
p4

g <- grid.arrange(p1, p2, p3, ncol = 2)

save_plot("PCA_plot_genotypeVsregion.pdf", g)

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
data_matrix <- assay(vsd) # Extract normalized, stabilized counts
pca.res <- prcomp(t(data_matrix), scale.=F, retx=T)
summary(pca.res)

pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleTable$sampleName,
             genotype = sampleTable$genotype,
             timepoint = sampleTable$timepoint,
             region = sampleTable$region,
             litter = sampleTable$litter,
             sex = sampleTable$sex)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)
## Genotype
pm1 <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=genotype) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) +
  theme_bw() +
  coord_flip()

## Timepoint
pm2 <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=timepoint) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) +
  theme_bw() +
  coord_flip()

## Region
pm3 <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=region) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) +
  theme_bw() +
  coord_flip()

## Sex and litter
pm4 <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=sex, color = litter) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption = paste0("produced on ", format(Sys.time(), "%x %X"))) +
  theme_bw() +
  coord_flip()

g <- grid.arrange(pm1, pm2, pm3, pm4, ncol = 2)


save_plot("PCA_small_multiples_plot.pdf", g)
