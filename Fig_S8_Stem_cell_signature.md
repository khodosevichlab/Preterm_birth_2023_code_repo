---
title: "stem cell signature"
author: "Laura Wolbeck"
date: "2024-07-02"
output: html_document
---


# prepare TERM2GENE data.frame

using genes from Cheung and Rando 
```{r}
quiescent_genes <- c("Ccnd3", "Pdk1","Smarca2", "Foxo3", "Ezh1", "Prdm5","Ptov1", "Zfp30", "Zbtb20","Phf1", "Ctdsp1", "Thra", "Tef", "Dicer1", "Crebrf", "Bcas3", "Ddx3y", "Gabarapl1", "Nop53", "Itm2a", "Il18", "Zyx", "Ephx1", "Clstn1", "Gstk1", "Cdip1", "Ddt", "Ivd", "Fhl1", "Ndrg2", "Grina","Pik3r1", "Fyn", "Chkb", "Pink1", "Ulk2", "Dnajb9", "Pfdn5", "Ctsf", "Crim1", "Selenop", "Gabbr1", "Grb10", "Bbs2", "Rps14", "Igf2r", "Selenbp1", "Rnf167", "Map1lc3a")
```

```{r}
active_genes <- c("Anln", "Birc5", "Ccna2","Ccnb1", "Ccne2", "Sgo1", "Mcm4", "Pcna", "Rrm2", "Top2a", "Cycs", "Mtch2", "Slc25a5", 	"H2afz", "Hat1", "Ddx39", 	"Pclaf", "Capza1", "Hadhb", "Idh3a", "Kpna2", "Pgk1")
```


TERM2GENE is a data.frame with first column of term ID and second column of corresponding mapped gene
```{r}
term2gene <- data.frame(Term = replicate(49, "Quiescent stem cell signature"), Gene = quiescent_genes)
term2gene <- rbind(term2gene, data.frame(Term = replicate(22, "Active stem cell signature"), Gene = active_genes))
saveRDS(term2gene, "term2gene.rds")
```

# Fig S8B: fullterm: E18.5 vs fullterm P2
```{r}
cao <- readRDS("cao_fullterm_DEG2.rds")
```

taking all DEG of RG
```{r}
RG_DEG <-  cao$test.results$de$RG$res
```

```{r}
RG_DEG <- RG_DEG$log2FoldChange
names(RG_DEG) <- cao$test.results$de$RG$res$Gene
#11144 genes
```

sort genes by log2fold change
```{r}
RG_DEG %<>% sort(decreasing = TRUE)
```
```{r}
set.seed(123)
res <- clusterProfiler::GSEA(RG_DEG, TERM2GENE = term2gene, pvalueCutoff = 1)
summary(res)
```
```{r}
clusterProfiler::dotplot(res)+ facet_grid(.~.sign) + ggtitle("Fullterm")
ggsave("Signature_fullterm.pdf", height=4)
```

# late timepoints: 
```{r}
cao <- readRDS("cao_late_DEG2.rds")
```
taking only DEG of RG
```{r}
RG_DEG <-  cao$test.results$de$RG$res
```

```{r}
RG_DEG <- RG_DEG$log2FoldChange
names(RG_DEG) <- cao$test.results$de$RG$res$Gene
```

```{r}
RG_DEG %<>% sort(decreasing = TRUE)
```
```{r}
set.seed(123)
res <- clusterProfiler::GSEA(RG_DEG, TERM2GENE = term2gene, pvalueCutoff = 1)
summary(res)
```

```{r}
clusterProfiler::dotplot(res)+ facet_grid(.~.sign) + ggtitle("Late timepoints")
ggsave("Signature_late.pdf", height=4)
```

# early timepoint
```{r}
cao <- readRDS("cao_early_DEG2.rds")
```
taking only DEG of RG
```{r}
RG_DEG <-  cao$test.results$de$RG$res
```
```{r}
head(RG_DEG)
```

```{r}
RG_DEG <- RG_DEG$log2FoldChange
names(RG_DEG) <- cao$test.results$de$RG$res$Gene
```
```{r}
x <- active_genes %in% names(RG_DEG) 
names(x)<- active_genes
x
```
```{r}
RG_DEG %<>% sort(decreasing = TRUE)
```
```{r}
set.seed(123)
res <- clusterProfiler::GSEA(RG_DEG, TERM2GENE = term2gene, pvalueCutoff = 1)
summary(res)
```

```{r}
clusterProfiler::dotplot(res)+ facet_grid(.~.sign) + ggtitle("Early timepoints")
ggsave("Signature_early.pdf", height=4)
```