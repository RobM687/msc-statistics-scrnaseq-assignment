<u>*Required packages:*</u>

```python
pip install 'scanpy[leiden]'

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
```

### Figure 2 Data Quality Assessment using Violin Plots
```python
# Calculate QC metrics
sc.pp.calculate_qc_metrics(hu, percent_top=[50],inplace=True, log1p=False, )
# Violin plot of QC metrics
sc.pl.violin(hu, ['n_genes_by_counts','total_counts'],jitter=0.4, multi_panel=True)

# Subsequent data filtering via slicing
# Remove cells that have a high number of detected genes (>6000) e.g. doublets or low-quality cells.
hu=hu[hu.obs.n_genes_by_counts <6000,:]
 # Remove cells that have a high total count of transcripts (>1,500,000)
hu=hu[hu.obs.total_counts <1500000,:]
# Violin plot of QC metrics (post QC filter)
sc.pl.violin(hu, ['n_genes_by_counts','total_counts'],jitter=0.4, multi_panel=True)
```

### Figure 3 HVG dispersion plots
```python
# Identifies genes whose variance is high relative to their mean expression across all cells
sc.pp.highly_variable_genes(hu, min_mean=0.0125, max_mean=3, min_disp=0.5)
# Plots all genes (dispersion v mean expression), provides visual confirmation of which genes passed the defined HVG threshold.
sc.pl.highly_variable_genes(hu)
```

### Figure 4 PCA Elbow Plot
```python
# PCA calculation
sc.tl.pca(hu, svd_solver='arpack')
# Plots PCA elbow plot
sc.pl.pca_variance_ratio(hu, log=True)
```

### Figure 5 UMAP visualisations of Leiden clustering and marker gene expressions
```python
# Perform leiden clustering
sc.tl.leiden(hu, resolution=0.2)
# Generates the UMAP plots based on the specified annotations
sc.pl.umap(hu, color=['Tissue', 'Age', 'Disease_stage', 'CLU'])
# Visualise leiden clustering alongside canonical marker genes:
# OVGP1 (Secretory Epithelial),
# PIFO (Ciliated Epithelial),
# FOXJ1 (Ciliated Epithelial),
# PTPRC (Immune),
# COL1A1 and DCN (Fibroblast),
# KRT17 (STIC lesion).
sc.pl.umap(
    hu,
    color=['leiden', 'OVGP1', 'PIFO', 'FOXJ1', 'PTPRC', 'COL1A1', 'DCN', 'KRT17'],
    frameon=False
)
```

### Figure 6 UMAP visualisation of annotated Leiden clusters
```python
# Map Leiden clusters to cell-type labels based on expression of canonical marker genes
celltypedict = {
        '0' : 'Secretory Epithelial',
        '1' : 'Secretory Epithelial',
        '2' : 'Ciliated Epithelial',
        '3' : 'STIC lesion',
        '4' : 'Immune',
        '5' : 'STIC lesion',
        '6' : 'Fibroblast'}

hu.obs['Celltype'] = hu.obs['leiden'].map(celltypedict)
# Visualise cell-types mapped to leiden clusters
sc.pl.umap(hu, color=['Celltype'])
```

<u>*Required packages:*</u>

```r
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BioQC")
BiocManager::install("limma")
install.packages("ggrepel")
install.packages("WebGestaltR")
install.packages("BayesFactor")
install.packages("HDInterval")
install.packages("dplyr")
install.packages("tibble")
install.packages('smotefamily')
install.packages("Boruta")
install.packages("randomForest")
install.packages("ROCR")
install.packages("caret")
library(dplyr); library(ggplot2); library(tidyr)
library(limma); library(ggrepel); library(WebGestaltR)
library(tibble); library(smotefamily); library(Boruta)
library(randomForest); library(ROCR); library(caret)
library(BioQC); library(BayesFactor); library(HDInterval)
```

### Figure 7 Volcano plot of the top 100 differentially expressed genes between SE1 and SE2 cells
```r
# Volcano plot of top 100 genes
top_df <- topTable(fit, coef = 2, number = 100)
top_df$significance <- ifelse(top_df$P.Value < 0.05 & abs(top_df$logFC) > 1, 'Significant','Not Significant')

vp <- ggplot(top_df, aes(x = logFC, y = -log10(P.Value), color = significance)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c('Not Significant'='gray','Significant'='purple')) +
  theme_minimal() +
  labs(title='Volcano Plot of Differential Expression', x='Log2 Fold Change', y='-Log10(p-value)') +
  theme(legend.title = element_blank(), legend.position='top')
# Adding gene names to singificantly differentially expressed genes
vp + geom_text_repel(aes(label = ifelse(significance=='Significant', rownames(top_df), '')), box.padding = 0.5, point.padding = 0.5, max.overlaps=10, size=3)
```

### Figure 8 Dot plot summary of GSEA results for SE1, SE2, and STIC lesion clusters
```r
# Label GSEA results with its corresponding cell type
SE1Result$Celltype <- 'Secretory Epithelial-1'
SE2Result$Celltype <- 'Secretory Epithelial-2'
STICResult$Celltype <- 'STIC lesion'

# Combine GSEA results into a single dataframe
gsea <- dplyr::bind_rows(SE1Result, SE2Result, STICResult)

# Dot plot configuration
ggplot(gsea, aes(x=Celltype, y=description)) +
  geom_point(aes(size=FDR, color=pValue)) +
  scale_color_gradient(low='red', high='blue', name='pValue') +
  scale_size_continuous(range=c(0.2,5)) +
  theme_minimal() +
  labs(
    title='Gene Set Enrichment Analysis',
    x='Cell Types',
    y='Pathway Description'
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.8))
```

### Figure 9 Boxplots of Shannon entropy for SE1, SE2, and STIC lesion cell populations
```r
# Copilot code used to colour and align plots horizontally, use common
# y-axis, display iteration points, label median.
## 1) Shared y-axis across panels
ylim_all <- range(
  c(compilation_SE1$Entropy, compilation_SE2$Entropy, compilation_STIC$Entropy),
  na.rm = TRUE
)

## 2) Layout and margins
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 0.5))  # save old par to restore later

## --- Panel 1: Secretory Epithelial-1 ---
boxplot(compilation_SE1$Entropy,
        main = "Secretory Epithelial-1",
        ylab = "Shannon Entropy",
        ylim = ylim_all,
        col  = "#A6CEE3",  # light blue
        border = "grey20")

## Add per-iteration points (jittered horizontally at x=1)
set.seed(1)  # for reproducible jitter
x_jit1 <- 1 + rnorm(length(compilation_SE1$Entropy), mean = 0, sd = 0.03)
points(x_jit1, compilation_SE1$Entropy,
       pch = 16, col = adjustcolor("grey30", alpha.f = 0.35), cex = 0.7)

## Add median label
med1 <- median(compilation_SE1$Entropy, na.rm = TRUE)
text(x = 1.35, y = med1, labels = sprintf("Median\n%.2f", med1), cex = 0.8, pos = 3)

## --- Panel 2: Secretory Epithelial-2 ---
boxplot(compilation_SE2$Entropy,
        main = "Secretory Epithelial-2",
        ylab = "",
        ylim = ylim_all,
        col  = "#B2DF8A",  # light green
        border = "grey20")

x_jit2 <- 1 + rnorm(length(compilation_SE2$Entropy), mean = 0, sd = 0.03)
points(x_jit2, compilation_SE2$Entropy,
       pch = 16, col = adjustcolor("grey30", alpha.f = 0.35), cex = 0.7)

med2 <- median(compilation_SE2$Entropy, na.rm = TRUE)
text(x = 1.35, y = med2, labels = sprintf("Median\n%.2f", med2), cex = 0.8, pos = 3)

## --- Panel 3: STIC lesion ---
boxplot(compilation_STIC$Entropy,
        main = "STIC lesion",
        ylab = "",
        ylim = ylim_all,
        col  = "#FB9A99",  # light red
        border = "grey20")

x_jit3 <- 1 + rnorm(length(compilation_STIC$Entropy), mean = 0, sd = 0.03)
points(x_jit3, compilation_STIC$Entropy,
       pch = 16, col = adjustcolor("grey30", alpha.f = 0.35), cex = 0.7)

med3 <- median(compilation_STIC$Entropy, na.rm = TRUE)
text(x = 1.35, y = med3, labels = sprintf("Median\n%.2f", med3), cex = 0.8, pos = 3)

## 3) Reset par
par(op)
```

### Figure 10 ROC curve and confusion matrix for the Random Forest classifier applied to the validation cohort
```r
# Confusion matrix
# Validation cohort
pred<-predict(rf,test.data)
confusionMatrix(pred,test.data$Disease_stage)

# ROC curve
# Validation AUC
pred1=predict(rf,newdata = test.data,type = "prob")
perf = prediction(pred1[,2],test.data$Disease_stage)
# 1. Area under curve
auc = performance(perf, "auc")
auc <- as.numeric(auc@y.values)
auc
rf$err.rate[nrow(rf$err.rate),1]
# 2. True Positive and Negative Rate
pred3 = performance(perf, "tpr","fpr")
# 3. Plot the ROC curve
plot(pred3,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
```
