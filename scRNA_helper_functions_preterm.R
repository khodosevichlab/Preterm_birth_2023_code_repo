
#Mitochondrial gene fraction per cell
mitoFraction <- function(con, species="human") {
  if(species=="human") lapply(con$samples, function(d) Matrix::rowSums(d$counts[,grep("MT-", colnames(d$counts))]) / Matrix::rowSums(d$counts)) %>% Reduce(c, .)
  else if(species=="mouse") lapply(con$samples, function(d) Matrix::rowSums(d$counts[,grep("mt-", colnames(d$counts))]) / Matrix::rowSums(d$counts)) %>% Reduce(c, .)
  else stop("Species must either be 'human' or 'mouse'.")
}


#estimate UMAP and find clusters
embedUMAP <- function(con,
                      min.dist=0.01,
                      spread=15,
                      min.prob.lower=1e-7,
                      method=leiden.community,
                      resolution=1,
                      min.group.size=50) {
  message("Creating UMAP embedding...")
  con$embedGraph(method="UMAP", 
                 min.dist=min.dist, 
                 spread=spread,
                 min.prob.lower=min.prob.lower)
  
  message("Estimating clusters...")
  con$findCommunities(method=leiden.community, resolution=resolution, min.group.size=min.group.size)
  
  return(con)
}

# buildConosGraph wraps the functions con$buildGraph() and embedUMAP() defined above
buildConosGraph <- function(con,
                          k.conos=15, 
                          k.self=15, 
                          space='PCA', 
                          ncomps=40,
                          n.odgenes=2e3,
                          matching.method='mNN', 
                          metric='angular', 
                          score.component.variance=T,
                          alignment.strength=0,
                          min.dist=0.01, 
                          spread=15,
                          min.prob.lower=1e-3,
                          resolution=1,
                          min.group.size=50) {
  message("Building graph...")
  con$buildGraph(k=k.conos, 
                 k.self=k.self, 
                 space=space, 
                 ncomps=ncomps, 
                 n.odgenes=n.odgenes, 
                 matching.method=matching.method, 
                 metric=metric, 
                 verbose=T, 
                 score.component.variance=score.component.variance,
                 alignment.strength=alignment.strength)
  
  embedUMAP(con=con,
            min.dist=min.dist,
            spread=spread,
            min.prob.lower=min.prob.lower,
            method=leiden.community,
            resolution=resolution,
            min.group.size=min.group.size)
  
  return(con)
}

# quickConos executes first pagoda::basicP2proc() on the list of cms, then initiates a new conos 
#object and then executes the wrapper function buildConosGraph() defined above
quickConos <- function(cms, 
                       sample.names,
                       n.cores.p2,
                       n.cores.con,
                       n.odgenes=3e3, 
                       nPcs = 50, 
                       k.p2 = 30, 
                       perplexity = 50, 
                       log.scale = TRUE, 
                       trim = 10, 
                       keep.genes = NULL, 
                       min.cells.per.gene = 3, 
                       min.transcripts.per.cell = 200, 
                       get.largevis = F, 
                       get.tsne = F, 
                       make.geneknn = TRUE,
                       k.conos=15, 
                       k.self=30, 
                       space='PCA', 
                       ncomps=40, 
                       matching.method='mNN', 
                       metric='angular', 
                       score.component.variance=T,
                       alignment.strength=0,
                       min.dist=0.01, 
                       spread=15) {
  if(length(cms)==length(sample.names)) {
    message("Performing P2 processing...")
    panel.preprocessed <- lapply(cms, function(x) basicP2proc(x, n.cores = n.cores.p2,
                                                              n.odgenes = n.odgenes, 
                                                              nPcs = nPcs,
                                                              k = k.p2, 
                                                              perplexity = perplexity, 
                                                              log.scale = log.scale, 
                                                              trim = trim, 
                                                              keep.genes = keep.genes, 
                                                              min.cells.per.gene = min.cells.per.gene, 
                                                              min.transcripts.per.cell = min.transcripts.per.cell, 
                                                              get.largevis = get.largevis, 
                                                              get.tsne = get.tsne, 
                                                              make.geneknn = make.geneknn))
    
    names(panel.preprocessed) = sample.names
    con <- Conos$new(panel.preprocessed, n.cores=n.cores.con)
    
    con <- buildConosGraph(con=con,
                           k.conos=k.conos, 
                           k.self=k.self, 
                           space=space, 
                           ncomps=ncomps, 
                           n.odgenes=n.odgenes, 
                           matching.method=matching.method, 
                           metric=metric, 
                           score.component.variance=score.component.variance,
                           alignment.strength=alignment.strength,
                           min.dist=min.dist, 
                           spread=spread)
    
    return(list(con=con, panel.preprocessed=panel.preprocessed))
  } else {
    stop("Sample names must match number of count matrices.")
  }
  
}


# extracts UMI counts per cell
getConosDepth <- function(con) {
  lapply(con$samples, function(d) d$depth) %>% unlist %>% setNames(.,(strsplit(names(.), ".", T) %>% 
                                                                                 sapply(function(d) d[2])))
}


# get conos clusters
getConosCluster <- function(con, name="leiden") {
  con$clusters[[name]]$groups
}


# rename a factor level in the annotation
renameAnnotation <- function(annotation, old, new) {
  if(!is.factor(annotation)) stop("Annotation must be a factor.")
  
  levels(annotation)[levels(annotation) %in% old] <- new
  
  return(annotation)
}

# collapse annotation
collapseAnnotation <- function(anno, label) {
  anno %<>% factor
  idx <- grepl(label,levels(anno))
  cat(paste0("Collapsing ",sum(idx)," labels containing '",label,"' in their name into one label.\n"))
  levels(anno)[idx] <- c(label)
  anno %<>% factor
  return(anno)
}