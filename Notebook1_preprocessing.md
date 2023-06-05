---
title: "Notebook1_preprocessing"
author: "Laura Wolbeck"
date: "2023-02-21"
---

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
```


# 1. Read in count matrices
Make list of count matrices (cms), early (E 18.5 and pre P0) and late (full P2 and pre P3) timepoints will be analysed separately
```{r, eval=FALSE}
#create vector with paths to filtered cms 
paths1 <- c(E18_5_1= "/data/home/lwolbeck/data/Japan_preterm_counts/E18_5_1/outs/filtered_feature_bc_matrix", 
           E18_5_2="/data/home/lwolbeck/data/Japan_preterm_counts/E18_5_2/outs/filtered_feature_bc_matrix",
          pre_P0_1= "/data/home/lwolbeck/data/Japan_preterm_counts/pre_P0_1/outs/filtered_feature_bc_matrix",
          pre_P0_2= "/data/home/lwolbeck/data/Japan_preterm_counts/pre_P0_2/outs/filtered_feature_bc_matrix")

paths2 <- c(full_P2_1 = "/data/home/lwolbeck/data/Japan_preterm_counts/full_P2_1/outs/filtered_feature_bc_matrix",
           full_P2_2 ="/data/home/lwolbeck/data/Japan_preterm_counts/full_P2_2/outs/filtered_feature_bc_matrix",
           pre_P3_1 = "/data/home/lwolbeck/data/Japan_preterm_counts/pre_P3_1/outs/filtered_feature_bc_matrix",
          pre_P3_2= "/data/home/lwolbeck/data/Japan_preterm_counts/pre_P3_2/outs/filtered_feature_bc_matrix")

#read in data in a parallel manner 
cms1 <- pagoda2::read.10x.matrices(paths1, n.cores=10)
cms2 <- pagoda2::read.10x.matrices(paths2, n.cores=10)

#to check if any colnames are duplicates
any(duplicated(unlist(lapply(cms1,colnames))))
any(duplicated(unlist(lapply(cms2,colnames))))
```


## plot summary metrics Fig S2 A and B
combine both matrices
```{r}
cms <- append(cms1, cms2)
```

using the CRMetrics package to plot summary metrics
```{r}
metadata <- data.frame(sample=names(cms), group= c("fullterm", "fullterm", "preterm", "preterm", "fullterm", "fullterm", "preterm", "preterm"))
crm <- CRMetrics$new(data.path="/data/Japan_preterm_counts", n.cores = 30, metadata = metadata)
```

plot Fig S2 A
```{r, fig.width=5, fig.height=3}
crm$plotSummaryMetrics(comp.group = "sample", second.comp.group = "group",metrics = "estimated number of cells",plot.geom = "bar") + xlab("") + ylab("# of cells") + theme(legend.position = "right", legend.title = element_blank())
ggsave("number_of_cells.pdf")
```
```{r}
crm$addDetailedMetrics()
```

plot Fig S2 B
```{r, fig.width=6}
crm$plotDetailedMetrics(comp.group = "group",
                        metrics = "gene_count", 
                        plot.geom = "violin",hline = F)  + ylab("genes/cell") + theme(legend.position = "right", legend.title = element_blank()) +  stat_summary(fun.y=median,size=10, geom = "text", label= "-")  + theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 12))
ggsave("VlnPlot_genes_per_cell.pdf")
```


# 2. Pagoda and conos processing
Conos processing with forced alignment (alignment.strength = 0.2) to integrate samples well
```{r}
#vector with sample names 
names1 = c("E18_5_1_",
           "E18_5_2_",
           "pre_P0_1_",
           "pre_P0_2_")

names2 = c("full_P2_1_",
           "full_P2_2_",
           "pre_P3_1_",
           "pre_P3_2_")
```

```{r, eval=FALSE}
#quickConos() is a wrapper function defined in scRNA_helper_functions_preterm.R
con1 <- quickConos(cms1,
                names1,
                  n.cores.p2=10,
                  n.cores.con=20, alignment.strength = 0.2)
con1 <- con1$con

con2 <- quickConos(cms2,
                  names2,
                  n.cores.p2=10,
                  n.cores.con=20, alignment.strength = 0.2)
con2 <- con2$con

```


# 3. Filter cells
## Run scrublet 
Create new folders scrublet 1 and 2  

Write files for Srublet
```{r, eval=FALSE}
#takes quite a while
mapply(function(x,y) write.table(as.matrix(x$misc$rawCounts), file = paste0("scrublet1/",y,".csv"), dec = ".", sep=","), x=con1$samples, y=names(con1$samples))

mapply(function(x,y) write.table(as.matrix(x$misc$rawCounts), file = paste0("scrublet2/",y,".csv"), dec = ".", sep=","), x=con2$samples, y=names(con2$samples))
```

Run Scrublet on each sample separately using the python script /d0-mendel/home/lwolbeck/preterm-NCU-UCPH/scrub.py:  
python scrub.py <sample>.csv

Import doublet scores 
```{r, eval=FALSE}
doubletvec1 <- do.call("rbind", lapply(dir(path = "scrublet1/", pattern=".doubletScores", full.names = T), read.csv)) %>% .$X0 %>% setNames(con1$samples %>% sapply(function(x) rownames(x$misc$rawCounts)) %>% unlist)

doubletvec2 <- do.call("rbind", lapply(dir(path = "scrublet2/", pattern=".doubletScores", full.names = T), read.csv)) %>% .$X0 %>% setNames(con2$samples %>% sapply(function(x) rownames(x$misc$rawCounts)) %>% unlist)
```


## 3.1 Remove doublets with scrublet score > 0.3
```{r}
singlets1 <- doubletvec1[doubletvec1>0.3] #758 cells (2.46%) are removed 
cms1 %<>% lapply(function(s) s[,!colnames(s) %in% names(singlets1)])

singlets2 <- doubletvec2[doubletvec2>0.3] #940 cells (2.95%) are removed
cms2 %<>% lapply(function(s) s[,!colnames(s) %in% names(singlets2)])
```

### plot Fig S2 C and D
```{r}
plot <- con1$plotGraph(groups = doubletvec1>0.3, size=0.1, mark.groups=F, show.legend = F, title="Doublets in early timepoints") +scale_colour_manual(values=c("grey", "red"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rasterize(plot, layers='Point', dpi=1000)
ggsave("Doublets_early.pdf", width=7, height=4.3 )
```
```{r}
plot <- con2$plotGraph(groups = doubletvec2>0.3, size=0.1, mark.groups=F, show.legend = F, title="Doublets in late timepoints") +scale_colour_manual(values=c("grey", "red"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rasterize(plot, layers='Point', dpi=1000)
ggsave("Doublets_late.pdf", width=7, height=4.3 )
```


## 3.2 Remove cells with mito fraction > 5 % 
Get mito fraction (percent of mt genes per cell) with helper function
```{r}
mito1 <- mitoFraction(con1, "mouse")
mito2 <- mitoFraction(con2, "mouse")
```
```{r}
cms1 %<>% lapply(function(s) s[,!colnames(s) %in% names(mito1[mito1>0.05])])
cms2 %<>% lapply(function(s) s[,!colnames(s) %in% names(mito2[mito2>0.05])])
```

## 3.3 Remove cells with depth < 1000
get depth with helper function
```{r}
depth1 <- getConosDepth(con1)
depth2 <- getConosDepth(con2)
```

```{r}
cms1 %<>% lapply(function(s) s[,!colnames(s) %in% names(depth1[depth1<1000])])
cms2 %<>% lapply(function(s) s[,!colnames(s) %in% names(depth2[depth2<1000])])
```

# 4. Pagoda and conos 2
build new conos object with filtered cms
```{r, eval=FALSE}
con1 <- quickConos(cms1,
                  names1,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = TRUE, alignment.strength=0.2)

con1 <- con1$con


con2 <- quickConos(cms2,
                  names2,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = TRUE, alignment.strength=0.2)

con2 <- con2$con
```

## Additional clusters
Rerun leiden clustering to find additional clusters
```{r, eval=FALSE}
con1$findCommunities(method = leiden.community, resolution = 5, min.group.size = 15)
con2$findCommunities(method = leiden.community, resolution = 5, min.group.size = 15)
```

```{r}
leiden24 <- con1$clusters$leiden$groups %>% factor
leiden24_2 <- con2$clusters$leiden$groups %>% factor
```

```{r, eval=F}
con1$clusters$leiden$groups <- leiden24
con2$clusters$leiden$groups <- leiden24_2
```

in late timepoints (full P2 and pre P3) cluster 20 was removed since it could not be annotated and was spread between several clusters
```{r}
cms2 %<>% lapply(function(s) s[,!colnames(s) %in% names(leiden24_2[leiden24_2=="20"])])
```

generating new UMAP embedding with cluster 20 removed
```{r}
con2 <- quickConos(cms2,
                  names2,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = TRUE, alignment.strength=0.2)

con2 <- con2$con
```

The final clusters were generated using a combination of manual selection and findSubcommunities() function of conos to split the clusters further. The clusters were annotated using known cell-type specific marker genes from the literature.





