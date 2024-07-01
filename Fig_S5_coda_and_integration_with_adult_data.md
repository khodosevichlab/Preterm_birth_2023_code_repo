---
title: "Fig S5"
author: "Laura Wolbeck"
date: "2024-07-01"
output: html_document
---
#setup
```{r}
#Load helper functions
source("scRNA_helper_functions_preterm.R")
library(Seurat)
library(qs)
library(pagoda2)
library(conos)
library(ggplot2)
library(magrittr)
library(ggrastr)
library(cowplot)
library(clusterProfiler)
```

# Fig S5A
## E18.5 vs full P2
```{r}
cao <- readRDS("cao_fullterm_DEG2.rds")
```

```{r}
cao$estimateCellLoadings()
```
```{r}
saveRDS(cao, "cao_fullterm_DEG2.rds" )
```

```{r, fig.height=6}
cao$plotCellLoadings(signif.threshold=0.05, show.pvals = F) + xlab("separation coefficient") 
ggsave("Compositionshifts_fullterm.pdf",width=6, height=4.5)
```

# Fig S5B
##  E18.5 vs Pre P0 (early timepoint)
```{r}
cao <- readRDS("cao_early_DEG2.rds")
```

```{r}
cao$estimateCellLoadings()
```

```{r, fig.height=6}
cao$plotCellLoadings(signif.threshold=0.05, show.pvals = F) + xlab("separation coefficient") 
ggsave("Compositionshifts_early.pdf",width=6, height=4.5)
```
```{r}
saveRDS(cao, "cao_early_DEG2.rds")
```

# Fig S5C
##  Full P2 vs Pre P3 (late timepoint)
```{r}
cao <- readRDS("cao_late_DEG2.rds")
```

```{r, fig.height=6}
cao$estimateCellLoadings()
```
```{r}
saveRDS(cao, "cao_late_DEG2.rds")
```

```{r, fig.height=6}
cao$plotCellLoadings(signif.threshold=0.05, show.pvals = F)  + xlab("separation coefficient") 
ggsave("Compositionshifts_late.pdf",width=6, height=4.5)
```

# combine our data (full-term P2) with P29/35 from Cebrian-Silla et al.
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165554
tar file of their seurat object: GSE165554_RAW.tar
untared: GSM5039270_scSeq.rds.gz
4 samples were labelled with Ab and sequenced together

Bar1 and Bar2 are males/p35
Bar3 and 4, females/p29

```{r}
seurat <- readRDS("GSM5039270_scSeq.rds")
```

get their raw count matrix from seurat object
```{r}
cms_raw <- GetAssayData(seurat, slot= "counts", assay = "RNA")
```

read in our list of all count matrices
```{r}
cms <- readRDS("cms_all.rds")
```

```{r}
# select only fullterm P2 count matrices
cmsx <- append(cms[5],cms[6])

cmsx <- append(cmsx, cms_raw)
names(cmsx)[3] <- "P29_P35"
```

## combined conos
```{r}
names = c("full_P2_1_",
          "full_P2_2_",
          "P29_P35")
```

```{r, eval=FALSE}
con02 <- quickConos(cmsx,
                  names,
                  n.cores.p2=40,
                  n.cores.con=40, get.tsne = T, alignment.strength = 0.2)
con02 <- con02$con
```

```{r, eval=FALSE}
#con02$findCommunities(method = leiden.community, resolution = 5, min.group.size = 50)
#leiden34 <- con02$clusters$leiden$groups %>% factor
#table(leiden27)
```
```{r}
qsave(con02, "con_P29_P35_fullP2.qs")
```

get annotation of adult dataset
```{r}
anno_P29_P35 <- seurat@meta.data$Cell_Type
names(anno_P29_P35) <- rownames(seurat@meta.data)
```
```{r}
anno_P29_P35 <- renameAnnotation(anno_P29_P35, "Microglia", "Microg")
```

read in annotation of our dataset
```{r}
anno_rough <- readRDS("anno_all7_rough.rds")
```

combine adult und fullterm P2 annotation
```{r}
anno_combined<- c(anno_rough, anno_P29_P35)
qsave(anno_combined, "anno_combined.qs")
```

# Fig S5D
```{r}
colours <- c("#FFAF32","#3232FF","#32FF32", "#FFFF32","#FF3232","#AF0F0F","#AF0F0F","#AF0F0F","#AF32FF","#E632E6","#FFC8FF","#DBFF00", "#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2") %>% setNames(c("Immature_ependymal_cells", "RG","Intermediate_progenitors", "Dividing_cells", "Neuroblasts", "Neurons_1","Neurons_2","Neurons_3","OPC","Microglia","Endothelial/Pericytes/VSMC","Erythrocytes", "Ependymal cells" , "Astrocytes","B cells", "C cells" , "Mitosis", "A cells","GABAergic neurons", "Neuron" , "OPC/Oligo","Microg","Endothelial cells" ,"Pericytes/VSMC","VLMC1"))
```

```{r}
plot <- con02$plotGraph(groups = anno_combined, size=0.2, label="") + scale_colour_manual(values=colours)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_combined_our3.pdf", width=7, height=4.3 )
```

# Fig S5E
```{r}
colours <- c("#FFAF32","#3232FF","#32FF32", "#FFFF32", "#FF3232","#AF0F0F","#AF0F0F","#AF32FF","#E632E6","#FFC8FF","#FFC8FF","#FFC8FF","#3232FF", "#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2","#d2d2d2") %>% setNames(c("Ependymal cells", "B cells", "C cells", "Mitosis",  "A cells","GABAergic neurons", "Neuron", "OPC/Oligo","Microg","Endothelial cells" ,"Pericytes/VSMC","VLMC1","Astrocytes", "Immature_ependymal_cells",  "RG","Intermediate_progenitors", "Dividing_cells", "Neuroblasts","Neurons_1","Neurons_2","Neurons_3","OPC","Microglia","Endothelial/Pericytes/VSMC","Erythrocytes"))
```

```{r}
plot <- con02$plotGraph(groups = anno_combined, size=0.2, label ="")+ scale_colour_manual(values=colours)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_combined_adult3.pdf", width=7, height=4.3 )
```

# Fig S5F
```{r}
#prepare condition factor
sample <- con02$getDatasetPerCell()
conditions <- gsub("full_P2_1_", "full_P2", sample)
conditions <- gsub("full_P2_2_", "full_P2", conditions)

conditions %<>% as.factor()

names(conditions) <- names(sample)
head(conditions)
```

rename some cell types in annotation for the barplot
```{r}
anno_combined_barplot <- anno_combined
anno_combined_barplot <- renameAnnotation(anno_combined_barplot, "A cells", "Neuroblasts")
anno_combined_barplot <- renameAnnotation(anno_combined_barplot, "B cells", "Radial glia/B cells")
anno_combined_barplot <- renameAnnotation(anno_combined_barplot, "RG", "Radial glia/b cells")
anno_combined_barplot <- renameAnnotation(anno_combined_barplot, "C cells", "Intermediate_progenitors")
anno_combined_barplot <- renameAnnotation(anno_combined_barplot, "Mitosis", "Dividing_cells")
```

annotation with only neurogenic lineage cell types
```{r}
anno_neurog <-anno_combined_barplot[anno_combined_barplot %in% c("Neuroblasts","Dividing_cells","Intermediate_progenitors","Radial glia/B cells")]
anno_neurog <- droplevels(anno_neurog)
```

```{r}
plotClusterBarplots(con02, groups = conditions , show.entropy = F, show.size = F, sample.factor = anno_neurog) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "cell type") + 
  theme(plot.margin = margin(0.5,0,0,2, "cm"))
ggsave("barplots_full_adult_neurogenesis_2.pdf")
```