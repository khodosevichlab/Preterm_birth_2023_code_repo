---
title: "Notebook2: all samples"
author: "Laura Wolbeck"
date: "2023-02-23"
---

In the previous notebook the samples were split into early and late timepoints. Here we will generated the UMAP embedding for all samples together and provide the code for the respective figures. 

# Setup
```{r setup, results='hide', warning=FALSE, message = FALSE}
#Load helper functions
source("scRNA_helper_functions_preterm.R")

#Load libraries
library(conos) #version 1.4.0
library(pagoda2)
library(magrittr)
library(ggplot2)
library(CRMetrics)
library(seurat)
library(cowplot)
```


```{r, eval=F, echo=F}
# cms1 and cms2 are loaded from notebook1
cms <- append(cms1, cms2)
```

```{r}
names = c("E18_5_1_",
          "E18_5_2_",
          "pre_P0_1",
          "pre_P0_2",
          "full_P2_1_",
          "full_P2_2_",
          "pre_P3_1",
          "pre_P3_2")
```

# 1. UMAP embedding and clustering 

```{r, eval=FALSE}
con <- quickConos(cms,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = TRUE, alignment.strength=0.3)
con <- con$con
```


## Additional clusters
Rerun leiden clustering to find additional clusters

```{r, eval=FALSE}
con$findCommunities(method = leiden.community, resolution = 5, min.group.size = 15)
leiden_all <- con$clusters$leiden$groups %>% factor
```

The final clusters were generated using a combination of manual selection and findSubcommunities() function of conos to split the clusters further. The clusters were annotated using known cell-type specific marker genes from the literature.

# 2. Generate seurat object
generate seurat object to use seurats plotting functions

get raw counts from filtered cms
```{r}
cm <- con$getJointCountMatrix(raw=T) %>% 
  Matrix::t() %>% 
  as.sparse()
```

Create Seurat object 
```{r}
seurat <- CreateSeuratObject(counts = cm, project = "All samples", min.cells = 0, min.features = 0)
```

add normalized cm
```{r}
cm_norm <- con$getJointCountMatrix(raw=F) %>% 
  Matrix::t() %>% 
  as.sparse()

seurat <- SetAssayData(object= seurat, slot= "data", cm_norm )
```

to make sure that cells in annotations are in same order as in cm
read in anno and anno_rough
```{r}
# order cells in "anno" as in "cm"
anno_rough %<>% .[match(colnames(cm), names(.))] 
anno %<>% .[match(colnames(cm), names(.))] 
```

add annotations as metadata
```{r}
seurat$annotation_rough <- anno_rough
seurat$annotation <- anno
```
create condition factor
```{r}
sample <- con$getDatasetPerCell()

conditions <- gsub("full_P2_1_", "fullterm", sample)
conditions <- gsub("full_P2_2_", "fullterm", conditions)
conditions <- gsub("pre_P3_1", "preterm", conditions)
conditions <- gsub("pre_P3_2", "preterm", conditions)
conditions <- gsub("E18_5_1_", "fullterm", conditions)
conditions <- gsub("E18_5_2_", "fullterm", conditions)
conditions <- gsub("pre_P0_1", "preterm", conditions)
conditions <- gsub("pre_P0_2", "preterm", conditions)

conditions %<>% as.factor()

names(conditions) <- names(sample)
```

add condition metadata to seurat object
```{r}
seurat$condition <- conditions
```

# 3. Prepare data for adata object (python)
save cm and metadata to be read into python
```{r}
cm <- con$getJointCountMatrix(raw=F) %>% 
  as.sparse()
writeMM(cm, "cm_all.mtx")
```
```{r}
genes <- colnames(cm)
write.csv(genes, 'all_genes.csv')
```

```{r}
# order cells in "anno" as in "cm"
anno_rough %<>% .[match(rownames(cm), names(.))] 
```
```{r}
#prepare condition factor
sample <- con$getDatasetPerCell())
conditions <- gsub("full_P2_1_", "full_P2", sample)
conditions <- gsub("full_P2_2_", "full_P2", conditions)
conditions <- gsub("pre_P3_1", "pre_P3", conditions)
conditions <- gsub("pre_P3_2", "pre_P3", conditions)
conditions <- gsub("E18_5_1_", "E18_5", conditions)
conditions <- gsub("E18_5_2_", "E18_5", conditions)
conditions <- gsub("pre_P0_1", "pre_P0", conditions)
conditions <- gsub("pre_P0_2", "pre_P0", conditions)
```

```{r}
#save metadata in a dataframe
metadata <- data.frame(annotation=as.character(anno_rough), cellnames=names(anno_rough), sample=con$getDatasetPerCell(),condition=conditions)
head(metadata)
write.csv(metadata, 'all_metadata.csv')
```


# Fig 1 E
```{r}
colours <- c("#FFAF32","#3232FF","#32FF32", "#FFFF32","#FF3232","#AF0F0F","#AF0F0F","#AF0F0F","#AF32FF","#E632E6","#FFC8FF","#DBFF00") %>% setNames(c("Immature_ependymal_cells",
"RG","Intermediate_progenitors", "Dividing_cells", "Neuroblasts", "Neurons_1","Neurons_2","Neurons_3","OPC","Microglia","Endothelial/Pericytes/VSMC","Erythrocytes"))
```
```{r}
con$plotGraph(groups=anno_rough, label="") + scale_colour_manual(values=colours)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```
# Fig 1 F
```{r}
cowplot::plot_grid(con$plotGraph(gene="Glul", title="Glul", size=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank()), con$plotGraph(gene="Gls", title="Gls", size=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank()))
```

# Fig 1 G
```{r}
colours <- c("#FF5D71", "#00B0F0", "#A9D18E", "#7030A0") %>% setNames(c("E18_5", "pre_P0", "full_P2", "pre_P3"))
```
```{r}
level_order <- c("E18_5", "pre_P0", "full_P2", "pre_P3")
Idents(seurat) <- "condition" 
seurat$condition <- factor(x = seurat$condition, levels = level_order)
```
```{r, fig.width=5, fig.height=8}
Idents(seurat) <- "annotation" 
VlnPlot(seurat, features="Glul",group.by = "condition", idents = "RG", cols = colours)  +
    stat_summary(fun.y = mean, geom='point', size = 3, colour = "yellow", show.legend = F) 
```
```{r, fig.width=5, fig.height=8}
VlnPlot(seurat, features="Gls",group.by = "condition", idents = "RG", cols = colours)  +
    stat_summary(fun.y = mean, geom='point', size = 3, colour = "yellow", show.legend = F) 
```

# Fig S3 A
```{r}
#defining colours of clusters
colours_high <- c("#e1af32","#FFAF32","#3232FF","#32FF32","#FFFF00","#CDCD32","#FF9b9b","#EB8C32","#C8961E","#FF3232","#AF3232","#B49191","#694646","#820000","#965050","#550505","#C80F0F","#AF32FF","#E632E6","#FFC8FF","#DBFF00") %>% setNames(c("Immature ependymal cells-2","Immature ependymal cells-2","RG","Intermediate progenitor cells","Dividing cells-2","Dividing cells-1","Neuroblasts-1", "Dividing cells-4", "Dividing cells-3","Neuroblasts-2","Neurons 1-1","Neurons 1-2","Neurons 2","Neurons 3-3","Neurons 3-4","Neurons 3-1","Neurons 3-2","OPC","Microglia", "Endothelial/Pericytes/VSMC","Erythrocytes"))
```
```{r}
con$plotGraph(groups=anno, label="") + scale_colour_manual(values=colours_high)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

# Fig S3 B
The following code is run in jupyter notebook (python)
```{python}
import scanpy as sc
import pandas as pd
import numpy as np
```
```{python}
# read in cm, genes and metadata
adata = sc.read_mtx("cm_all.mtx")
gene_df = pd.read_csv("all_genes.csv")
metadata = pd.read_csv("all_metadata.csv", index_col=0)
```
```{python}
# preparing adata object
# adding gene names and cell names and metadata
adata.var_names = gene_df.x.values
adata.obs_names = metadata.cellnames.values
adata.obs = metadata
```
```{python}
cluster_marker = ["Lrrc23", "Foxj1", "Slc1a3", "Fabp7", "Ascl1", "Mki67", "Cenpa", "Dcx","Rbfox3","Neurod6", "Nkx2-1","Ebf1","Drd2","Olig2","Mpeg1","Epas1", "Hbb-bs"]
sc.pl.dotplot(adata, cluster_marker , groupby='annotation', swap_axes=True, categories_order = ['Immature_ependymal_cells', 'RG' , 'Intermediate_progenitors','Dividing_cells', 'Neuroblasts','Neurons_1','Neurons_2','Neurons_3', 'OPC', 'Microglia',  'Endothelial/Pericytes/VSMC', 'Erythrocytes'], save="all.pdf")
```

# Fig S3 C

```{r}
genes <- c("Foxj1", "Slc1a3","Ascl1", "Mki67", "Dcx", "Rbfox3", "Olig2", "Mpeg1", "Epas1")
genes %>% lapply(function(p) con$plotGraph(gene=p, title=p, size=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank())) %>%
  cowplot::plot_grid(plotlist=., ncol=5)
```

# Fig S3 D
```{r, fig.width=5, fig.height=8}
VlnPlot(seurat, features="Glud1",group.by = "condition", idents = "RG", cols = colours)  +
    stat_summary(fun.y = mean, geom='point', size = 3, colour = "yellow", show.legend = F) 
```

# Fig S4

defining level orders
```{r}
level_order_rough <- c("Immature_ependymal_cells", "RG","Intermediate_progenitors", "Dividing_cells","Neuroblasts", "Neurons_1", "Neurons_2", "Neurons_3" , "OPC" ,"Microglia", "Endothelial/Pericytes/VSMC", "Erythrocytes" )
```
```{r}
level_order <- c("Immature ependymal cells-1", "Immature ependymal cells-2", "Radial glia","Intermediate progenitor cells", "Dividing cells-1", "Dividing cells-2","Dividing cells-3", "Dividing cells-4","Neuroblasts-1", "Neuroblasts-2","Neurons 1-1", "Neurons 1-2","Neurons 2", "Neurons 3-1" , "Neurons 3-2", "Neurons 3-3","Neurons 3-4",  "OPC" ,"Microglia", "Endothelial/Pericytes/VSMC", "Erythrocytes" )
```

```{r}
#for A, C, E and G
Idents(seurat) <- "annotation_rough" 
seurat$annotation_rough <- factor(x = seurat$annotation_rough, levels = level_order_rough)
```
```{r}
#for B, D, F and H
Idents(seurat) <- "annotation" 
seurat$annotation <- factor(x = seurat$annotation, levels = level_order)
```

## Fig S4 A and B
```{r}
VlnPlot(seurat, features= "nFeature_RNA", pt.size = 0, log = F,) +  ggtitle(" ") + ylab("genes  per cell") + xlab(" ")+ theme_bw(base_size = 16) + theme( legend.position = "none") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 12500)+
  theme(plot.margin = margin(0.5,0,0,1, "cm"))
```
## Fig S4 C and D
```{r}
VlnPlot(seurat, features= "nCount_RNA", pt.size = 0)  + ggtitle(" ") +ylab("UMIs per cell") + xlab(" ") + theme_bw(base_size = 16) + theme( legend.position = "none") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.margin = margin(0.5,0,0,1, "cm")) + ylim(0, 100000)
```



## Fig S4 E
reading in low resolution annotation "anno_rough"
```{r, fig.height=4.5}
plotClusterBarplots(con, groups = anno_rough , show.entropy = F, show.size = F, sample.factor = conditions ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "condition") + 
  scale_x_discrete(limits = level_order_rough)+ 
  theme(plot.margin = margin(0.5,0,0,2, "cm"))
```
## Fig S4 F
read in annotation file high resolution as "anno"
```{r, fig.height= 4.5}
plotClusterBarplots(con, groups = anno , show.entropy = F, show.size = F, sample.factor = conditions) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "condition") + 
  scale_x_discrete(limits = level_order)+ 
  theme(plot.margin = margin(0.5,0,0,2, "cm"))
```
## Fig S4 G
```{r, fig.height= 4.5}
plotClusterBarplots(con, groups = anno_rough , show.entropy = F, show.size = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  scale_x_discrete(limits = level_order_rough) + 
  theme(plot.margin = margin(0.5,0,0,2, "cm"))
```

## Fig S4 H
```{r,  fig.height=4.5}
plotClusterBarplots(con, groups = anno, show.entropy = F, show.size = F, ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(limits = level_order)+ 
  theme(plot.margin = margin(0.5,0,0,2, "cm"))
```



