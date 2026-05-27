library(Seurat)
library(dplyr)
library(tibble)
library(reshape2)
require(scales)
library(ggpubr)
library(tidyr)
library(ggplot2)


#' 
#' @author David John
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCellsByQuantile <- function(seuratObject, lowQuantile=0.1 , highQuantile=0.95, maxMito=0.1){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  cat("FilterByQuantile\n")
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error 
  lowQuantile<<-lowQuantile
  highQuantile<<-highQuantile
  maxMito<<-maxMito
  sample<-unique(seuratObject$sample)
  Quality <- data.frame(UMI=seuratObject$nCount_RNA, nGene=seuratObject$nFeature_RNA, label = factor(seuratObject$sample), percent.mito=seuratObject$percent.mito)
  
  Quantile.low.UMI <- Quality %>% group_by(label) %>%
    summarise(UMI = list(enframe(quantile(UMI,probs = lowQuantile)))) %>%
    unnest(cols = c(UMI))
  
  Quantile.high.UMI <- Quality %>% group_by(label) %>%
    summarise(UMI = list(enframe(quantile(UMI,probs = highQuantile)))) %>%
    unnest(cols = c(UMI))
  
  Quantile.low.Gene <- Quality %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = lowQuantile)))) %>%
    unnest(cols = c(nGene))
  
  Quantile.high.Gene <- Quality %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = highQuantile)))) %>%
    unnest(cols = c(nGene))
  
  
  df<-seuratObject@meta.data
  
  gg1<- ggplot(Quality, aes(x="nUMI", y=UMI)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = Quantile.high.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.UMI$value, label=Quantile.high.UMI$value , vjust = -1)) + 
    geom_hline(yintercept = Quantile.low.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.UMI$value, label=Quantile.low.UMI$value , vjust = -1))
  
  gg2<- ggplot(Quality, aes(x="nFeature_RNA", y=nGene)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = Quantile.high.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.Gene$value, label=Quantile.high.Gene$value , vjust = -1)) +   geom_hline(yintercept = Quantile.low.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.Gene$value, label=Quantile.low.Gene$value , vjust = -1))
  
  
  gg3<- ggplot(Quality, aes(x=" % Mt Content", y=percent.mito)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = maxMito, color="red", linetype="dashed") + geom_text(aes(0.9,maxMito, label=maxMito , vjust = -1))
  
  gg<-ggarrange(gg1,gg2,gg3, ncol = 3)
  
  library(ggpubr)  
  
  gg<-annotate_figure(gg, fig.lab = sample, fig.lab.pos = "top", fig.lab.size = 15, fig.lab.face = "bold")
  
  seuratObject<- subset(x= seuratObject, subset = nCount_RNA < Quantile.high.UMI$value & nCount_RNA > Quantile.low.UMI$value & 
                          nFeature_RNA < Quantile.high.Gene$value & nFeature_RNA > Quantile.low.Gene$value & percent.mito < maxMito)
  
  
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  percent <- round(diff/sizeBefore*100,digits = 2)
  cat("Filtered ",diff, "from" , sizeBefore, " cells [",percent,"%]\n", "(minFeatures=",Quantile.low.Gene$value, "; maxFeatures=", Quantile.high.Gene$value, "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(maxMito)
  
  
  
  
  return(list(seuratObject, gg))
}





#' Import Single cell sequencing experiments into Seurat perform normalisation and scale Data and do a summary of mapping stats, optional perform Doubletfinder
#' @author David John & Mariano Ruz Jurado
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return list with Merged seurat object & statistics
Importer <- function(pathway,id, TenX=FALSE, CellBender=FALSE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE, FilterCells=TRUE, FilterByAbsoluteValues=FALSE, DeleteDoublets=FALSE, ...) {
  require(Seurat)
  require(ggplot2)
  
  if (TenX) {
    print("Read Cellranger V5")
    Matrix <- Read10X(pathway)
  }
  else if (CellBender) {
    scC <- system.file(package = "scCustomize") # Make sure package is installed
    ifelse(nzchar(scC), "", stop("Install scCustomize R package for Cellbender read in!"))
    file_path <- list.files(pathway, pattern = "*filtered.h5", full.names = TRUE)  #grab file
    Matrix <- scCustomize::Read_CellBender_h5_Mat(file_path)    
  }
  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  
  #catch optional parameters 
  optionalParameters <- list(...)
  
  print("Create SeuratObject")
  seuratObject =CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  print("Get MT Content")
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")
  
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObject), value = TRUE)
  if (length(mito.features)<10) {
    mito.features <- grep(pattern = "^mt-", x = rownames(x = seuratObject), value = TRUE)
  }
  if (length(mito.features)<1) {
    warning("Error: Could not find MT genes")
    
  }
  
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, layer = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObject, layer = 'counts'))
  seuratObject$percent.mito <- percent.mito
  
  print("Write QC images")
  #write QC to file
  p1<-VlnPlot(object = seuratObject, features = c("nFeature_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank()) #TODO This causes warnings when layer is not specified and using Seuratv5
  p2<-VlnPlot(object = seuratObject, features = c("nCount_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank()) # but it still works fine, just annoying
  p3<-VlnPlot(object = seuratObject, features = c("percent.mito"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  gg_preFiltering <- ggarrange(p1,p2,p3, nrow = 1)
  gg_preFiltering <- annotate_figure(gg_preFiltering, top = text_grob(id,face="bold",color = "darkred",size=18,hjust = 0.2))
  ggsave(filename = paste0(pathway,id,"_QC_preFiltered.svg"),device = "svg", width = 10,height = 10)
  
  if (FilterCells==TRUE) {
    print("start Filtering")
    if (FilterByAbsoluteValues==TRUE) {
      if (is.null(optionalParameters$minFeatures)) {
        stop("Please define 'minFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxFeatures)) {
        stop("Please define 'maxFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$minCounts)) {
        stop("Please define 'minCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxCounts)) {
        stop("Please define 'maxCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxMito)) {
        stop("Please define 'maxMito' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      print("Running FilterDeadCells")
      seuratObject<-FilterDeadCells(seuratObject = seuratObject,
                                    minFeatures = optionalParameters$minFeatures,
                                    maxFeatures = optionalParameters$maxFeatures,
                                    minCounts = optionalParameters$minCounts,
                                    maxCounts = optionalParameters$maxCounts,
                                    maxMito = optionalParameters$maxMito)
    }
    else {
      print("Running FilterByQuantile")
      tmp<-FilterDeadCellsByQuantile(seuratObject = seuratObject, lowQuantile = 0.1, highQuantile = 0.95)
      seuratObject<-tmp[[1]]
      svg(paste0(pathway,"QC_QuantileFiltering.svg"))
      print(tmp[[2]])
      dev.off()
      gg_preFiltering<-tmp[[2]]
    }
    
  }
  
  
  if (performNormalisation==TRUE) {
    print("Running Normalisation")
    seuratObject<-NormalizeData(object = seuratObject,verbose = FALSE)
  }
  if(performVariableGeneDetection==TRUE){
    print("Running Variable Gene Detection")
    seuratObject<-FindVariableFeatures(object = seuratObject, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  if (performScaling==TRUE) {
    print("Running ScaleData")
    seuratObject<-ScaleData(object = seuratObject)
  }
  if (DeleteDoublets==TRUE){
    print("Running scDblFinder")
    require(scDblFinder)
    require(Matrix)
    
    SCE_obj <- as.SingleCellExperiment(seuratObject)
    SCE_obj <- scDblFinder(SCE_obj)
    seuratObject$scDblFinder_score <- SCE_obj$scDblFinder.score
    seuratObject$scDblFinder_class <- SCE_obj$scDblFinder.class #singlet 60801  doublet  7029
    seuratObject <- subset(seuratObject, scDblFinder_class == "singlet") 
  }
  
  print(paste("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n"))
  
  
  return(list(seuratObject, gg_preFiltering))
}

#' @author Mariano Ruz Jurado
#' @param pathway A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param id Vector of strings that are assigned to the concordant cells
#' @return a list with summary about mapping results 
SummarizeMapping <- function(pathway,id){
  ## mapping stats, i hope everyone uses cellranger and starsolo directories for these analysis, else no summary
  if (file.exists(paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv"))) {
    metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv")
    #define your own numeric class
    setClass('myNum')
    #define conversion
    setAs("character", "myNum",
          function(from) as.numeric(gsub(",","", gsub("%","",from))))
    #read data with custom colClasses
    metrics_summ <- read.csv(metrics_summ.path,
                             header = T,   
                             stringsAsFactors=FALSE,
                             colClasses=c("myNum"))
    
    typeof(metrics_summ$Fraction.Reads.in.Cells)
    
    metrics_col <- as.data.frame(colnames(metrics_summ))
    rownames(metrics_col) <- metrics_col[,1]
    metrics_col[,1] <- as.character(as.vector(metrics_summ[1,]))
    metrics_summ <- metrics_col
    
    # warnings CELLRANGER
    if (metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),] < 70) {
      warning(paste0(id,": Reads mapped confidently to genome only ", metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),] < 30) {
      warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),]))
    }
    if (paste(unlist(strsplit(metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),],split=",")),collapse = "") < 1000) {
      warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),]), " ,maybe the 0s were cut because of CR way of displaying numbers,  if unsure check CR web_summary")
    }
    if (as.numeric(paste(unlist(strsplit(metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),],split=",")),collapse = "")) < 300) {
      warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),])," ,maybe the 0s were cut because of CR way of displaying numbers, if unsure check CR web_summary")
    }
  } else {
    metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "Gene"))[1],"Gene/Summary.csv")
    metrics_summ <- read.delim2(metrics_summ.path, header = F, sep = ",")
    rownames(metrics_summ) <- metrics_summ[,1]
    metrics_summ[,1] <- NULL
    
    # warnings STAR
    if (metrics_summ[7,] < 0.70) { # mapped to genome, no grep since same name as other row 
      warning(paste0(id,": Reads mapped confidently to genome only ",metrics_summ[7,]))
    }
    if (metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),] < 0.30) {
      warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),] < 1000) {
      warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),]))
    }
    if (metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),] < 300) {
      warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),]))
    }
  } 
  return(metrics_summ)
}

#' 
#' @author David John & Mariano Ruz Jurado
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCells <- function(seuratObject, minFeatures=300, maxFeatures=6000,minCounts=500,maxCounts=15000, maxMito=0.05){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error 
  minFeatures<<-minFeatures
  maxFeatures<<-maxFeatures
  maxMito<<-maxMito
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > minFeatures & nFeature_RNA < maxFeatures & nCount_RNA > minCounts & nCount_RNA < maxCounts & percent.mito < maxMito)
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  percent <- round(diff/sizeBefore*100,digits = 2)
  cat("Filtered ",diff, "from" , sizeBefore, " cells [",percent,"%]\n", "(minFeatures=",minFeatures, "; maxFeatures=", maxFeatures, "; minCounts=" ,minCounts,  "; maxCounts=" ,maxCounts , "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(minFeatures,maxFeatures,maxMito)
  return(seuratObject)
}

#' Import and combine several Single cell sequencing experiments into Seurat
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineSeuratObjects <- function(pathways,ids, performNormalisation = FALSE, performVariableGeneDetection=FALSE){
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  for (i in 1:length(pathways)) {
    if (i<2) { 
      seuratObject1<-Importer(pathways[i],ids[i], TenX = TRUE, performNormalisation = performNormalisation, performVariableGeneDetection = performVariableGeneDetection)
      next
    }
    seuratObject2<-Importer(pathways[i],ids[i])
    seuratObject1 <- merge(x = seuratObject1, y = seuratObject2) 
  }
  cat("Merged Seurat object contains ", length(seuratObject1@meta.data$orig.ident)," cells\n")
  return(seuratObject1)
}


#' Draw a FeatureHeatmap
#' @author Lukas Tombor and David John
#' @param object A Seurat3 object
#' @param features Genes which are plotted
#' @param group.by Split/Group TSNES by this parameter
#' @param cols Colors of visible cells
FeatureHeatmap <- function(object, features, group.by=NULL, cols=c("skyblue","red4"), assay="RNA") {
  require(reshape2);require(ggplot2)
  DefaultAssay(object)<-"RNA"
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Grouping parameter was not found in the meta data")
  }
  
  #test if Genes can be found
  for(i in 1:length(features)){
    if (!features[i] %in% rownames(object@assays$RNA@data)) {
      stop("Gene was not found in the Expression Matrix")
    }
  }
  
  
    A <- data.frame(object@meta.data)
    X <- Embeddings(object = object, reduction = "umap")
    coord = NULL
    for(i in rownames(A)){
      coord <- rbind(coord, c(X[i,], i))
    }
    nclusters<-length(table(A[,group.by]))
    A <- data.frame(A, coord)
    A$tSNE_1 <- as.numeric(levels(A$tSNE_1)[A$tSNE_1])
    A$tSNE_2 <- as.numeric(levels(A$tSNE_2)[A$tSNE_2])
    A$seurat_clusters <- factor(A$seurat_clusters, levels = 0:(nclusters-1))
    #A$sample_dummy <- as.integer(as.factor(A$sample))
    
    for(i in 1:length(features)){
      a.colnames<-colnames(A)
      a.colnames<-c(a.colnames,paste0(features[i],".exprs"))
      A <- data.frame(A, x=GetAssayData(object, assay = assay)[features[i], ])
      colnames(A)<-a.colnames
    }
    
    A.rep<-as.data.frame(lapply(A, rep, nclusters))
    f3 <- function(variables) {
      return(names(table(A.rep[,group.by])))
    }
    
    A.rep$Cluster_fate <- ave(rep(NA,nrow(A.rep)), A.rep[,group.by], FUN = f3)
    
    #add new column with NA for unfitting cluster and expression value for the correct cluster
    for (f in features) {
      A.rep[,f]<-ifelse(A.rep[,group.by]==A.rep$Cluster_fate,A.rep[,paste0(f,".exprs")],NA)
    }
    
    A.melt<-melt(A.rep, measure.vars = colnames(A.rep[,c((ncol(A.rep)-1):ncol(A.rep))]))
    
    ggplot(A.melt, aes(x=tSNE_1, y= tSNE_2, color = value))+
      geom_point(size = 0.2)+
      scale_color_continuous(low = cols[1], high = cols[2], na.value = "grey93", name="Scaled Expression")+
      facet_grid(c("variable","Cluster_fate"), scales = "free", drop = FALSE)+
      labs(y="Genes", x= group.by, title = "")+
      theme_bw()+
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", strip.text = element_text(face = "italic"))
    
}





# Draw Barplot of percentages per group
#' @author David John
#' @param SeuratObject A Seurat3 object
#' @param split.by slot for bars in each plot
#' @param group.by Split/Group Barplots by this parameter
#' @param color.by Color bars by meta.data slot
BarplotPercentages <- function(SeuratObject, split.by="orig.ident", group.by="seurat_clusters", color.by="condition", returnValues=FALSE){
  require(stringr)
  require(ggplot2)
  meta.data<- SeuratObject@meta.data
  orig.ident.ordered<-str_sort(unlist(unique(SeuratObject[[split.by]])), numeric = TRUE)
  
  V<-data.frame(split.by=factor(meta.data[[split.by]],levels = orig.ident.ordered),
                group.by=factor(meta.data[[group.by]]),
                color.by=factor(meta.data[[color.by]]))
  
  
  Summary.Celltypes <- V %>% dplyr::count(split.by,group.by, color.by,.drop = FALSE) %>% group_by(split.by) %>% 
    mutate(freq = n /sum(n)) %>% complete(group.by,fill = list(n=0,freq=0)) 
  
  if (returnValues==TRUE) {
    return(Summary.Celltypes)
  }
  
  ggplot(Summary.Celltypes, aes(x=split.by, y= freq, fill= color.by))+
    geom_col(width = 0.9, color = "black")+
    facet_wrap(~group.by, scales = "free")+
    scale_y_continuous(name = "Percent per Celltype", labels = scales::percent_format())+
    theme(panel.background = element_blank(),
          strip.text = element_text(size=12),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust= 1, size = 10))
}

## prepare data for cluster t.test from the deg list and do a cluster t-test
do_cluster_t_test <- function(seurat_subset, degs, group="condition", cluster="seurat_clusters"){
  gene_names<- names(table(rownames(degs)))
  #print(head(gene_names))
  p_values <- vector("list",length(gene_names))
  names(p_values) <- gene_names
  #gene_names <- row.names(cluster_subset)
  #if (celltype=="Adipocytes"){
  #  seurat_subset <- seurat_subset[,!seurat_subset$orig.ident=="D7"]
  #}
  group <- seurat_subset[[group]][,1]
  cluster <- seurat_subset[[cluster]][,1]
  for (gene in gene_names){
    y <- c(t(as.matrix(seurat_subset@assays$RNA[gene,])))
    test_info <- my.t.test.cluster(y, cluster, group)
    p_values[[gene]] <- test_info[nrow(test_info)]
  }
  return(p_values)
}

## added line 54-56 so that each group is tested if 
## only one observation is present and throw an error
my.t.test.cluster <- function (y, cluster, group, conf.int = 0.95) 
{
  group <- as.factor(group)
  cluster <- as.factor(cluster)
  s <- !(is.na(y) | is.na(cluster) | is.na(group))
  y <- y[s]
  cluster <- cluster[s]
  group <- group[s]
  n <- length(y)
  if (n < 2) 
    stop("n<2")
  gr <- levels(group)
  if (length(gr) != 2) 
    stop("must have exactly two treatment groups")
  n <- table(group)
  nc <- tapply(cluster, group, function(x) length(unique(x)))
  bar <- tapply(y, group, mean)
  u <- unclass(group)
  y1 <- y[u == 1]
  y2 <- y[u == 2]
  c1 <- factor(cluster[u == 1])
  c2 <- factor(cluster[u == 2])
  b1 <- tapply(y1, c1, mean)
  b2 <- tapply(y2, c2, mean)
  m1 <- table(c1)
  m2 <- table(c2)
  if (any(names(m1) != names(b1)))
    stop("logic error 1")
  if (any(names(m2) != names(b2)))
    stop("logic error 2")
  if (any(m2 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m2[m2 < 2]), collapse = " ")))
  if (any(m1 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m1[m1 < 2]), collapse = " ")))
  M1 <- mean(y1)
  M2 <- mean(y2)
  ssc1 <- sum(m1 * ((b1 - M1)^2))
  ssc2 <- sum(m2 * ((b2 - M2)^2))
  if (nc[1] != length(m1))
    stop("logic error 3")
  if (nc[2] != length(m2))
    stop("logic error 4")
  df.msc <- sum(nc) - 2
  msc <- (ssc1 + ssc2)/df.msc
  v1 <- tapply(y1, c1, var)
  v2 <- tapply(y2, c2, var)
  ssw1 <- sum((m1 - 1) * v1)
  ssw2 <- sum((m2 - 1) * v2)
  df.mse <- sum(n) - sum(nc)
  mse <- (ssw1 + ssw2)/df.mse
  na <- (sum(n) - (sum(m1^2)/n[1] + sum(m2^2)/n[2]))/(sum(nc) - 
                                                        1)
  rho <- (msc - mse)/(msc + (na - 1) * mse)
  r <- max(rho, 0)
  C1 <- sum(m1 * (1 + (m1 - 1) * r))/n[1]
  C2 <- sum(m2 * (1 + (m2 - 1) * r))/n[2]
  v <- mse * (C1/n[1] + C2/n[2])
  v.unadj <- mse * (1/n[1] + 1/n[2])
  de <- v/v.unadj
  dif <- diff(bar)
  se <- sqrt(v)
  zcrit <- qnorm((1 + conf.int)/2)
  cl <- c(dif - zcrit * se, dif + zcrit * se)
  z <- dif/se
  P <- 2 * pnorm(-abs(z))
  stats <- matrix(NA, nrow = 20, ncol = 2, dimnames = list(c("N", 
                                                             "Clusters", "Mean", "SS among clusters within groups", 
                                                             "SS within clusters within groups", "MS among clusters within groups", 
                                                             "d.f.", "MS within clusters within groups", "d.f.", "Na", 
                                                             "Intracluster correlation", "Variance Correction Factor", 
                                                             "Variance of effect", "Variance without cluster adjustment", 
                                                             "Design Effect", "Effect (Difference in Means)", "S.E. of Effect", 
                                                             paste(format(conf.int), "Confidence limits"), "Z Statistic", 
                                                             "2-sided P Value"), gr))
  stats[1, ] <- n
  stats[2, ] <- nc
  stats[3, ] <- bar
  stats[4, ] <- c(ssc1, ssc2)
  stats[5, ] <- c(ssw1, ssw2)
  stats[6, 1] <- msc
  stats[7, 1] <- df.msc
  stats[8, 1] <- mse
  stats[9, 1] <- df.mse
  stats[10, 1] <- na
  stats[11, 1] <- rho
  stats[12, ] <- c(C1, C2)
  stats[13, 1] <- v
  stats[14, 1] <- v.unadj
  stats[15, 1] <- de
  stats[16, 1] <- dif
  stats[17, 1] <- se
  stats[18, ] <- cl
  stats[19, 1] <- z
  stats[20, 1] <- P
  attr(stats, "class") <- "t.test.cluster"
  stats
}

#perform SEM Graphs
#' @author Mariano Ruz Jurado
#' @param object # combined object
#' @param Features # vector containing featurenames
#' @param ListTest # List for which conditions t-test will be performed, if NULL always against provided CTRL 
#' @param returnValues # return df.melt.sum data frame containing means and SEM for the set group
#' @param ctrl.condition # set your ctrl condition, relevant if running with empty comparison List
#' @param group.by # select the seurat object slot where your conditions can be found, default conditon
DO.Mean.SEM.Graphs.cluster.t <- function(object, Features, ListTest=NULL, returnValues=FALSE, ctrl.condition=NULL, group.by = "condition", returnPlot=FALSE){ 
  require(ggplot2)
  require(ggpubr)
  require(tidyverse)
  require(reshape2)
  print("Please use 'DO.Mean.SEM.Graphs.wilcox' for Seurat wilcox Test and Seuratv5 Support.")
  #SEM function defintion
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided object, aswell as original identifier of samples
  df<-data.frame(condition=setNames(object[[group.by]][,group.by], rownames(object[[group.by]]))
                 ,orig.ident = object$orig.ident)
  #get expression values for genes from individual cells, add to df
  for(i in Features){
    df[,i] <- expm1(object@assays$RNA$data[i,])
    
  }
  
  #melt results 
  df.melt <- melt(df)
  #group results and summarize, also add/use SEM 
  df.melt.sum <- df.melt %>% 
    dplyr::group_by(condition, variable) %>% 
    dplyr::summarise(Mean = mean(value))
  #second dataframe containing mean values for individual samples
  df.melt.orig <- df.melt %>% 
    dplyr::group_by(condition, variable, orig.ident) %>% 
    dplyr::summarise(Mean = mean(value))
  
  
  #add SEM calculated over sample means
  df.melt.sum$SEM <- NA
  for (condition in df.melt.sum$condition) {
    df.melt.orig.con <- df.melt.orig[df.melt.orig$condition %in% condition,] # condition wise
    for (gene in Features) {
      df.melt.orig.con.gen <- df.melt.orig.con[df.melt.orig.con$variable %in% gene,] #gene wise
      df.melt.sum[df.melt.sum$condition %in% condition & df.melt.sum$variable %in% gene,]$SEM <- SEM(df.melt.orig.con.gen$Mean)
    }
  }
  
  #create comparison list for t.test, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list 
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with provided control")
    conditions <- unique(object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) { 
      ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","WT","Wt","wt"),collapse ="|")
                                        ,conditions)[1]]
    }
    
    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in 1:length(conditions)) {
      cndtn <- conditions[i] 
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- as.character(c(ctrl.condition,cndtn))
      }
    }
  }
  #delete Null values, created by count index
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  #create barplot with significance
  p<-ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = Mean, ymax = Mean+SEM), width = 0.2)+
    geom_point(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1, position = "jitter")+
    #ordering, control always first
    scale_x_discrete(limits=c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
    #t-test, always against control, using means from orig sample identifier
    stat_compare_means(data=df.melt.orig, comparisons = ListTest, method = "t.test", size=3)+
    facet_wrap(~variable, ncol = 9, scales = "free") +
    scale_fill_manual(values = rep(c("royalblue" ,"forestgreen", "tomato", "sandybrown"),4) #20 colours set for more change number
                      , name = "Condition")+
    labs(title = "", y = "Mean UMI") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",angle = 45,hjust = 1, size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 14, color = "black"),
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.line = element_line(color = "black"),
          strip.text.x = element_text(size = 14, color = "black"),
          legend.position = "bottom")
  print(p)
  if (returnValues==TRUE) {
    return(df.melt.sum)
  }
  if (returnPlot==TRUE) {
    return(p)
  }
}

#perform SEM Graphs
#' @author Mariano Ruz Jurado
#' @param object # combined object
#' @param Feature # name of the feature/gene
#' @param ListTest # List for which conditions wilcox will be performed, if NULL always CTRL group against everything 
#' @param returnValues # return data frames needed for the plot, containing df.melt, df.melt.sum, df.melt.orig and wilcoxstats
#' @param ctrl.condition # set your ctrl condition, relevant if running with empty comparison List
#' @param group.by # select the seurat object slot where your conditions can be found, default conditon
#' @param bar_colours # colour vector 
#' @param plotPvalue # plot the actual p-value without adjusting for multiple tests
#' @param SeuV5 # Seuratv5
DO.Mean.SEM.Graphs.wilcox <- function(object, Feature, ListTest=NULL, 
                                      returnValues=FALSE, ctrl.condition=NULL, 
                                      group.by = "condition", wilcox_test=TRUE, 
                                      bar_colours=NULL, stat_pos_mod = 1.15, 
                                      x_label_rotation=45, plotPvalue=FALSE, y_limits = NULL,
                                      log1p_nUMI=T){ 
  require(ggplot2)
  require(ggpubr)
  require(tidyverse)
  require(reshape2)
  
  if (!(Feature %in% rownames(object)) && !(Feature %in% names(object@meta.data))) {
    stop("Feature not found in Seurat Object!")
  }
  
  if (wilcox_test == T) {
    rstat <- system.file(package = "rstatix") # Make sure package is installed
    ifelse(nzchar(rstat), "", stop("Install rstatix R package for wilcox statistic!"))    
  }
  
  #SEM function defintion
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided object, aswell as original identifier of samples
  df<-data.frame(condition=setNames(object[[group.by]][,group.by], rownames(object[[group.by]]))
                 ,orig.ident = object$orig.ident)
  #get expression values for genes from individual cells, add to df
  #if (SeuV5==F) {
  #  for(i in Feature){
  #    df[,i] <- expm1(object@assays$RNA@data[i,])
  #    
  #  }    
  #}
  #For Seuratv5 where everything is a layer now
  #if (SeuV5==T) {
  #rlang::warn("\nSeuV5 set to TRUE, if working with Seuratv4 or below change SeuV5 to FALSE", .frequency = "once", .frequency_id = "v5Mean")
  if (Feature %in% rownames(object)) {
    df[,Feature] <- expm1(FetchData(object, vars = Feature))
  }else{
    df[,Feature] <- FetchData(object, vars = Feature)
  }
  #}
  
  # stat.df$condition <- factor(stat.df$condition)
  # stat.df$variable <- factor(stat.df$variable)
  
  # stat.df.test <- data.frame(Mean = object[["RNA"]]@data["Figf",],condition = object[[group.by]][,group.by])
  
  #melt results 
  df.melt <- melt(df)
  #group results and summarize, also add/use SEM 
  df.melt.sum <- df.melt %>% 
    dplyr::group_by(condition, variable) %>% 
    dplyr::summarise(Mean = mean(value))
  #second dataframe containing mean values for individual samples
  df.melt.orig <- df.melt %>% 
    dplyr::group_by(condition, variable, orig.ident) %>% 
    dplyr::summarise(Mean = mean(value))
  
  if (Feature %in% names(object@meta.data)) {
    plot.title <- paste("Mean", Feature, sep = " ")
  }  else if (log1p_nUMI==T) {
    df.melt.sum$Mean <- log1p(df.melt.sum$Mean)
    df.melt.orig$Mean <- log1p(df.melt.orig$Mean)
    plot.title <- "Mean log(nUMI)"
  } else{
    plot.title <- "Mean nUMI"
  }
  
  
  #create comparison list for wilcox, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list 
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with each other\n")
    conditions <- unique(object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) { 
      ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","ctrl","WT","Wt","wt"),collapse ="|")
                                        ,conditions)[1]]
    }
    
    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in 1:length(conditions)) {
      cndtn <- as.character(conditions[i]) 
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }
  
  #delete Null values, created by count index also reorder for betetr p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  indices <- sapply(ListTest, function(x) match(x[2], df.melt.sum$condition))
  ListTest <- ListTest[order(indices)]
  
  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$condition %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
      }
    }
    return(lst_filtered)
  }
  
  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df.melt.sum)
  
  
  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE) {
    stat.test <- df.melt %>%
      ungroup() %>%
      rstatix::wilcox_test(value ~ condition, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
  }
  #add SEM calculated over sample means
  df.melt.sum$SEM <- NA
  for (condition in df.melt.sum$condition) {
    df.melt.orig.con <- df.melt.orig[df.melt.orig$condition %in% condition,] # condition wise
    for (gene in Feature) {
      df.melt.orig.con.gen <- df.melt.orig.con[df.melt.orig.con$variable %in% gene,] #gene wise
      df.melt.sum[df.melt.sum$condition %in% condition & df.melt.sum$variable %in% gene,]$SEM <- SEM(df.melt.orig.con.gen$Mean)
    }
  }
  
  if (is.null(bar_colours)) {
    bar_colours <- rep(c("#1f77b4","#ea7e1eff","royalblue4","tomato2","darkgoldenrod","palegreen4","maroon","thistle3"),5)#20 colours set for more change number
  }
  
  if (x_label_rotation == 45) {
    hjust <- 1
  } else{hjust <- 0.5}
  
  #create barplot with significance
  p<-ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = Mean, ymax = Mean+SEM), width = 0.2)+
    geom_point(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1, position = "jitter")+
    #ordering, control always first
    scale_x_discrete(limits=c(as.character(ctrl.condition),levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
    #t-test, always against control, using means from orig sample identifier
    #stat_compare_means(data=stat.df,aes(x=condition, y=Mean, group=variable),comparisons = ListTest,method = "wilcox",size=3, label.y = max(df.melt.sum$Mean)*1.4)+
    facet_wrap(~variable, ncol = 9, scales = "free") +
    scale_fill_manual(values =  bar_colours 
                      , name = "Condition")+
    labs(title = "", y = plot.title) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black",angle = x_label_rotation,hjust = hjust, size = 16),
          axis.text.y = element_text(color = "black", size = 16),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 16, hjust = 0.5),
          axis.line = element_line(color = "black"),
          strip.text.x = element_text(size = 16, color = "black"),
          legend.position = "none")
  if (!is.null(y_limits)) {
    p = p + ylim(y_limits) 
  }
  if (wilcox_test == TRUE) {
    
    #Adjustments when ylim is changed manually
    y_pos_test <- max(df.melt.orig$Mean)*stat_pos_mod
    if (!is.null(y_limits) && y_pos_test > max(y_limits)) {
      y_pos_test <- max(y_limits)* stat_pos_mod - 0.1 * diff(y_limits)     
    }
    if (plotPvalue==TRUE) {
      p = p + stat_pvalue_manual(stat.test, label = "p = {p}", y.position = y_pos_test, step.increase = 0.2)
    }
    else{
      p = p + stat_pvalue_manual(stat.test, label = "p = {p.adj}", y.position = y_pos_test, step.increase = 0.2)
    }
  }
  
  #print(p)
  
  if (returnValues==TRUE) {
    returnList <- list(p, df.melt, df.melt.orig, df.melt.sum, stat.test)
    names(returnList) <- c("plot","df.melt", "df.melt.orig", "df.melt.sum", "stat.test")
    return(returnList)
  }
  return(p)
}



#perform better Violins than Seurat
#' @author Mariano Ruz Jurado
#' @param object # combined object
#' @param Feature # name of the feature
#' @param ListTest # List for which conditions wilcox will be performed, if NULL always CTRL group against everything 
#' @param returnValues # return df.melt.sum data frame containing means and SEM for the set group
#' @param ctrl.condition # set your ctrl condition, relevant if running with empty comparison List
#' @param group.by # select the seurat object slot where your conditions can be found, default conditon
#' @param group.by.2 # relevant for multiple group testing, e.g. for each cell type the test between each of them in two conditions provided
#' @param geom_jitter_args # vector for dots visualisation in vlnplot: size, width, alpha value
#' @param vector_colours # specify a minimum number of colours as you have entries in your condition, default 2
#' @param wilcox_test #Bolean if TRUE a bonferoni wilcoxon test will be carried out between ctrl.condition and the rest
#' @param stat_pos_mod #value for modifiyng statistics height
#' @param SeuV5 # Seuratv5

DO.Vln.Plot.wilcox <- function(object,
                               SeuV5=T,
                               Feature,
                               ListTest=NULL,
                               returnValues=FALSE,
                               ctrl.condition=NULL,
                               group.by = "condition",
                               group.by.2 = NULL,
                               geom_jitter_args = c(0.20, 0.25, 0.25),
                               geom_jitter_args_group_by2 = c(0.1, 0.1, 1),
                               vector_colors = c("#1f77b4","#ea7e1eff","royalblue4","tomato2","darkgoldenrod","palegreen4","maroon","thistle3"),
                               wilcox_test = T,
                               stat_pos_mod = 1.15,
                               hjust.wilcox=0.8,
                               vjust.wilcox = 2.0,
                               size.wilcox=3.33,
                               step_mod=0,
                               hjust.wilcox.2=0.5,
                               vjust.wilcox.2=0,
                               width_errorbar=0.4){ 
  require(ggplot2)
  require(ggpubr)
  require(tidyverse)
  require(reshape2)  
  require(cowplot)
  
  if (!(Feature %in% rownames(object)) && !(Feature %in% names(object@meta.data))) {
    stop("Feature not found in Seurat Object!")
  }
  
  if (wilcox_test == T) {
    rstat <- system.file(package = "rstatix") # Make sure package is installed
    ifelse(nzchar(rstat), "", stop("Install rstatix R package for wilcox statistic!"))    
  }
  
  if (is.null(ctrl.condition)) {
    stop("Please specify the ctrl condition as string!")    
  }
  
  if (SeuV5==T) {
    rlang::warn("SeuV5 set to TRUE, if working with Seuratv4 or below change SeuV5 to FALSE", .frequency = "once", .frequency_id = "v5Mean")
    
    if (Feature %in% rownames(object)) {
      vln.df = data.frame(Feature = object[["RNA"]]$data[Feature,],
                          cluster = object[[group.by]])
    }else{
      vln.df = data.frame(Feature = FetchData(object, vars = Feature)[,1],
                          cluster = object[[group.by]])
    }
    
    
    
    df<-data.frame(group=setNames(object[[group.by]][,group.by], rownames(object[[group.by]]))
                   ,orig.ident = object$orig.ident)
    
    # add a second group for individual splitting and testing in the wilcoxon 
    if (!is.null(group.by.2)) {
      vln.df[group.by.2] <- object[[group.by.2]]
      df[group.by.2] <- object[[group.by.2]]
    }    
    
    #get expression values for genes from individual cells, add to df
    if (Feature %in% rownames(object)) {
      df[,Feature] <- object@assays$RNA$data[Feature,]
    }else{
      df[,Feature] <- FetchData(object, vars = Feature)[,1]
    }
    
    
    
  }
  
  if (SeuV5==F) {
    if (Feature %in% rownames(object)) {
      vln.df = data.frame(Feature = object[["RNA"]]@data[Feature,],
                          cluster = object[[group.by]])
    }else{
      vln.df = data.frame(Feature = FetchData(object, vars = Feature)[,1],
                          cluster = object[[group.by]])
    }
    
    df<-data.frame(group=setNames(object[[group.by]][,group.by], rownames(object[[group.by]]))
                   ,orig.ident = object$orig.ident)
    # add a second group for individual splitting and testing in the wilcoxon
    if (!is.null(group.by.2)) {
      vln.df[group.by.2] <- object[[group.by.2]]
      df[group.by.2] <- object[[group.by.2]]
    }
    
    #get expression values for genes from individual cells, add to df
    if (Feature %in% rownames(object)) {
      df[,Feature] <- object@assays$RNA@data[Feature,]
    }else{
      df[,Feature] <- FetchData(object, vars = Feature)[,1]
    }
  }
  
  
  df.melt <- melt(df)
  
  vln.df$group <- factor(vln.df[[group.by]]
                         ,levels = c(as.character(ctrl.condition),levels(factor(vln.df[[group.by]]))[!(levels(factor(vln.df[[group.by]])) %in% ctrl.condition)]))
  #create comparison list for wilcox, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list 
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with each other")
    group <- unique(object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) { 
      ctrl.condition <- group[grep(pattern = paste(c("CTRL","Ctrl","ctrl","WT","Wt","wt"),collapse ="|")
                                   ,group)[1]]
    }
    
    
    #create ListTest
    ListTest <- list()
    for (i in 1:length(group)) {
      cndtn <- as.character(group[i]) 
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }
  
  #delete Null values, created by count index also reorder for better p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  if (!is.null(group.by.2)) {
    indices <- sapply(ListTest, function(x) match(x[2], vln.df[[group.by.2]]))    
  } else{
    indices <- sapply(ListTest, function(x) match(x[2], vln.df[[group.by]]))
  }
  ListTest <- ListTest[order(indices)]
  
  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$group %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
      }
    }
    return(lst_filtered)
  }
  
  #group results and summarize
  if (is.null(group.by.2)) {
    df.melt.sum <- df.melt %>% 
      dplyr::group_by(group, variable) %>% 
      dplyr::summarise(Mean = mean(value))   
  } else{
    df.melt.sum <- df.melt %>% 
      dplyr::group_by(group, !!sym(group.by.2), variable) %>% 
      dplyr::summarise(Mean = mean(value))  
  }
  
  
  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df.melt.sum)  
  
  if (!is.null(group.by.2) && length(ListTest) > 1) {
    stop("The provided Seurat has more than two groups in group.by and you specified group.by.2, currently not supported (to crowded)!")
  }
  
  
  # artificially set a group to 0 in all their expression values to create the error in the test 
  # df.melt$value[df.melt$group == "Human-HFrEF"] <- 0
  # df.melt$value[df.melt$group == "Human-CTRLlv"] <- 0
  # df.melt$value[df.melt$group == "Human-HFrEF" & df.melt$cell_type == "Cardiomyocytes"] <- 0
  # df.melt$value[df.melt$group == "Human-CTRLlv" & df.melt$cell_type == "Cardiomyocytes"] <- 0
  
  #check before test if there are groups in the data which contain only 0 values and therefore let the test fail
  if (is.null(group.by.2)) {
    group_of_zero <- df.melt %>%
      dplyr::group_by(group) %>%
      summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      filter(all_zeros)
    
    if (nrow(group_of_zero) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df.melt <- df.melt %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(value = if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup() 
    }
  } else{
    
    group_of_zero <- df.melt %>%
      dplyr::group_by(group, !!sym(group.by.2)) %>%
      summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      filter(all_zeros)
    
    #check now the result for multiple entries in group.by.2
    groupby2_check <- group_of_zero %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      summarise(group_count = n_distinct(group), .groups = "drop") %>%
      filter(group_count > 1)
    
    if (nrow(groupby2_check) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df.melt <- df.melt %>%
        dplyr::group_by(group, !!sym(group.by.2)) %>%
        dplyr::mutate(value = if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup()   
    }
  }
  
  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE & is.null(group.by.2)) {
    stat.test <- df.melt %>%
      ungroup() %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p))
    
  }
  
  if (wilcox_test == TRUE & !is.null(group.by.2)) {
    stat.test <- df.melt %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p))
    
  }
  
  if (length(unique(vln.df[[group.by]])) >  length(vector_colors)) {
    stop(paste0("Only ", length(vector_colors)," colors provided, but ", length(unique(vln.df[[group.by]])), " needed!"))
  }
  
  #normal violin
  if(is.null(group.by.2)){
    p <- ggplot(vln.df, aes(x = group, y = Feature))+ 
      geom_violin(aes(fill = group), trim = T, scale = "width", )+
      geom_jitter(size = geom_jitter_args[1], width = geom_jitter_args[2], alpha = geom_jitter_args[3])+
      labs(title = Feature, y = "Expression Level")+
      xlab("")+
      ylab("")+
      theme_classic()+
      theme(plot.title = element_text(face = "bold", color = "black", hjust = 0.5, size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
            legend.position = "none")+
      scale_fill_manual(values = vector_colors)
    
    if (!Feature %in% rownames(object)) {
      p_label="p = {p.adj}"
    } else{
      p_label="p = {p}"
    }
    
    if (wilcox_test == TRUE) {
      p = p + stat_pvalue_manual(stat.test, label = p_label, y.position = max(vln.df$Feature)*1.15, step.increase = 0.2)
    }
    return(p)
  }
  
  
  if(!is.null(group.by.2)){
    #plot
    p <- ggplot(vln.df, aes(x = !!sym(group.by.2), y = Feature, fill = !!sym(group.by)))+
      geom_violin(aes(fill = group), trim = T, scale = "width")+
      labs(title = Feature, y = "log(nUMI)")+
      xlab("")+
      theme_classic()+
      theme(plot.title = element_text(face = "bold", color = "black", hjus = 0.5, size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
            legend.position = "bottom",
            panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
            axis.line = element_line(colour = "black"),
            strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
            strip.text = element_text(colour = "black", size = 12),
      )+
      scale_fill_manual(values = vector_colors)
    
    p2 <- ggplot(vln.df, aes(x = !!sym(group.by.2), y = Feature, fill = factor(!!sym(group.by),levels = levels(vln.df$group))))+
      geom_boxplot(width=.1,color="grey", position = position_dodge(width = 0.9), outlier.shape = NA)+
      xlab("")+
      scale_fill_manual(values = c("black","black"), name=group.by)+
      theme_classic()+
      theme(plot.title = element_text(face = "bold", color = "transparent", hjus = 0.5, size = 14),
            axis.title.y = element_text(face = "bold", color = "transparent", size = 14),
            axis.text.x = element_text(face = "bold", color = "transparent", angle = 45, hjust = 1, size = 14),
            axis.text.y = element_text(face = "bold", color = "transparent", hjust = 1, size = 14),
            legend.position = "bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            strip.background = element_rect(fill = "transparent", colour = NA),
            axis.ticks = element_blank(),    
            legend.background = element_rect(fill = "transparent", color = NA),  # Transparent legend background
            legend.key = element_rect(fill = "transparent", color = NA),          # Transparent legend keys
            legend.title = element_text(face = "bold", color = "transparent"),          # Legend title styling
            legend.text = element_text(color = "transparent"),
            # strip.text = element_blank(),
      )
    
    if (wilcox_test == TRUE & !is.null(group.by.2)) {
      
      if (Feature %in% rownames(object)) {
        stat.test_plot <- stat.test %>%
          mutate(y.position = seq(from = max(object@assays$RNA$data[Feature,][!is.na(object@assays$RNA$data[Feature,])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
          mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
                 xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
                 xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)                
      } else{
        stat.test_plot <- stat.test %>%
          mutate(y.position = seq(from = max(object[[Feature]][,1][!is.na(object[[Feature]][,1])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
          mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
                 xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
                 xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)  
      }
      
      # dplyr::select(x.axis, y.position, p.adj)
      
      #if Feature not a gene than use the uncorrected p
      if (Feature %in% rownames(object)) {
        p_label="p = {p.adj}"
      } else{
        p_label="p = {p}"
      }
      
      p = p + stat_pvalue_manual(stat.test_plot,
                                 label = p_label,
                                 y.position = "y.position",
                                 # x="x",
                                 xmin = "xmin",
                                 xmax = "xmax",
                                 # xend="xend",
                                 # step.increase = 0.2,
                                 inherit.aes = FALSE,
                                 size = size.wilcox,
                                 angle= 0,
                                 hjust= hjust.wilcox.2,
                                 vjust = vjust.wilcox.2,
                                 tip.length = 0.02,
                                 bracket.size = 0.8)
      
      p2 = p2 + stat_pvalue_manual(stat.test_plot,
                                   label = p_label,
                                   y.position = "y.position",
                                   # x="x",
                                   xmin = "xmin",
                                   xmax = "xmax",
                                   # xend="xend",
                                   # step.increase = 0.2,
                                   inherit.aes = FALSE,
                                   size = size.wilcox,
                                   angle= 0,
                                   hjust= hjust.wilcox.2,
                                   vjust = vjust.wilcox.2,
                                   tip.length = 0.02,
                                   bracket.size = 0.8,
                                   color = "transparent")
      
    }
    plot_p <- ggdraw() + draw_plot(p) + draw_plot(p2)
    print(plot_p)
  }
  
  
  
  if (returnValues==TRUE) {
    returnList <- list(vln.df, df.melt, stat.test)
    names(returnList) <- c("vln.df", "df.melt", "stat.test")
    return(returnList)
  }
  
  return(plot_p)
}

#' @author Mariano Ruz Jurado & David Rodriguez Morales(for celltypist questions)
#' @title DO Celltypist
#' @description Runs celltypist on a seurat object and stores the calls as metadata
#' @param object The seurat object
#' @param minCellsToRun If the input seurat object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.
#' @param runCelltypistUpdate If true, --update-models will be run for celltypist prior to scoring cells.
#' @param over_clustering Column in metadata in object with clustering assignments for cells, default seurat_clusters
#' @param assay_normalized Assay with log1p normalized expressions
#' @param returnProb will additionally return the probability matrix, return will give a list with the first element beeing the object and second prob matrix

DO.CellTypist <- function(object, modelName = "Healthy_Adult_Heart.pkl", minCellsToRun = 200, runCelltypistUpdate = TRUE, over_clustering = "seurat_clusters" , assay_normalized = "RNA", returnProb=FALSE, SeuV5=T) {
  
  # Make sure R Zellkonverter package is installed
  zk <- system.file(package = "zellkonverter") 
  ifelse(nzchar(zk), "", stop("Install zellkonverter R package for Seurat tranformation!"))
  
  # Make sure R reticulate package is installed
  rt <- system.file(package = "reticulate") 
  ifelse(nzchar(rt), "", stop("Install reticulate R package for Python usage in R!"))
  
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not correctly configured. Run "usethis::edit_r_environ()" to specify your Python instance'))
  }
  
  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed in this python environment!')
  }
  
  if (ncol(object) < minCellsToRun) {
    warning('Too few cells, will not run celltypist. NAs will be added instead')
    expectedCols <- c('predicted_labels_celltypist')
    object[[expectedCols]] <- NA
    return(object)
  }
  
  message(paste0('Running celltypist using model: ', modelName))
  
  #Create temporary folder
  outDir <- tempfile(fileext = '')
  if (endsWith(outDir, "/")){
    outDir <- gsub(outDir, pattern = "/$", replacement = "")
  }
  dir.create(outDir)
  message(paste0("Saving celltypist results to temporary folder: ", outDir))
  
  #Uppercasing gene names
  #zellkonverter h5ad
  DefaultAssay(object) <- assay_normalized
  if (SeuV5 == TRUE) {
    tmp.assay <- object
    tmp.assay[["RNA"]] <- as(tmp.assay[["RNA"]], Class = "Assay")
    tmp.sce <- Seurat::as.SingleCellExperiment(tmp.assay, assay = assay_normalized)
    rownames(tmp.sce) <- toupper(rownames(tmp.sce))
    
  } else{
    tmp.sce <- Seurat::as.SingleCellExperiment(object, assay = assay_normalized)
    rownames(tmp.sce) <- toupper(rownames(tmp.sce))
  }
  zellkonverter::writeH5AD(tmp.sce, file = paste0(outDir,"/ad.h5ad"), X_name = "logcounts")
  
  # Ensure models present:
  if (runCelltypistUpdate) {
    system2(reticulate::py_exe(), c("-m", "celltypist.command_line", "--update-models", "--quiet"))
  }
  
  # Run celltypist:
  args <- c("-m", "celltypist.command_line", "--indata",  paste0(outDir,"/ad.h5ad"), "--model", modelName, "--outdir", outDir,"--majority-voting", "--over-clustering", over_clustering)
  system2(reticulate::py_exe(), args)
  
  labelFile <- paste0(outDir, "/predicted_labels.csv")
  probFile <- paste0(outDir, "/probability_matrix.csv")
  labels <- utils::read.csv(labelFile, header = T, row.names = 1, stringsAsFactors = FALSE)
  
  #ad <- import(module = "anndata")
  #ad_obj <- ad$AnnData(X = labels)
  
  #ct <- import(module = "celltypist")
  #ct$dotplot(ad_obj, use_as_reference = "cell_type", use_as_prediction = "majority_voting")
  
  probMatrix <- utils::read.csv(probFile, header = T, row.names = 1, stringsAsFactors = FALSE)
  object@meta.data$predicted_labels_celltypist <- labels$majority_voting
  if (returnProb==TRUE) {
    returnProb <- list(object, probMatrix)
    names(returnProb) <- c("SeuratObject", "probMatrix")
    return(returnProb)
  } else {
    return(object)}
}

# WIP 
#' @author Mariano Ruz Jurado
#' @title DO Three Variable Dotplot
#' @description Creates a dotplot with three informations retrieved from a seurat object meta data, works per gene
#' @param seuratObj The seurat object
#' @param group.by.x x-axis group name to look for in meta data
#' @param group.by.y1 first y-axis group name to look for in meta data
#' @param group.by.y2 second y-axis group name to look for in meta data
#' @param dot.size Vector of dot size 
#' @param plot.margin = plot margins
#' @param Feature Gene of interest
#' @param colZ IF True calculates the Z-score of the average expression per column

DO.Three.Variable.Dotplot <- function(SeuratObject, Feature, group.by.x = NULL, group.by.y1 = NULL, group.by.y2 = NULL,
                                      dot.size = c(1,6), plot.margin = c(1, 1, 1, 1), midpoint = 0.5, colZ=FALSE, median=F){
  
  # Create Feature expression data frame with grouping information
  print("DO.Three.Variable.Dotplot is deprecated, please use DO.Dotplot!")
  geneExp <- expm1(Seurat::FetchData(object = SeuratObject, vars = Feature, layer = "data")) # This  might not work for Seuratv5...
  geneExp$id <- paste(SeuratObject@meta.data[[group.by.y1]], " (", 
                      SeuratObject@meta.data[[group.by.y2]], ")", sep = "")
  geneExp$xaxis <- SeuratObject@meta.data[[group.by.x]]
  
  
  # geneExp %<>% filter(if_any(1, ~ .x > 0)) 
  
  # Include xaxis in the overall grouping
  data.plot <- lapply(X = unique(geneExp$id), FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, ]
    
    lapply(X = unique(data.use$xaxis), FUN = function(x_axis) {
      data.cell <- data.use[data.use$xaxis == x_axis, 1:(ncol(geneExp) - 2), drop = FALSE]
      avg.exp <- apply(X = data.cell, MARGIN = 2, FUN = function(x) {
        if (median==T) {
          funct_calc <- stats::median          
        } else{
          funct_calc <- base::mean
        }
        return(expm1(funct_calc(x)))
      })
      pct.exp <- apply(X = data.cell, MARGIN = 2, FUN = PercentAbove, 
                       threshold = 0)
      
      res <- data.frame(id = ident, xaxis = x_axis, avg.exp = avg.exp, pct.exp = pct.exp * 100)
      res$gene <- rownames(res)
      return(res)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .) %>% data.frame()
  
  # get the scale dvalue for plotting 
  data.plot$avg.exp.log1p <- log1p(data.plot$avg.exp) # reapply the log transformation
  
  data.plot.res <- data.plot
  data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = sort(unique(data.plot.res$xaxis)))
  data.plot.res$group <- sapply(strsplit(as.character(data.plot.res$id), 
                                         split = "\\(|\\)"), "[", 2)
  
  # Check if Others is one of the groups to factorize as last item
  if (any(grepl("Other", unique(data.plot.res$id)))) {
    data.plot.res$id <- factor(data.plot.res$id, levels = c(grep("Other", unique(data.plot.res$id), value = T), grep("Other", unique(data.plot.res$id), value = T, invert = T)))
  }
  data.plot.res$pct.exp <- ifelse(data.plot.res$pct.exp == 0, NA, data.plot.res$pct.exp) # so fraction 0 is not displayed in plot
  data.plot.res <- data.plot.res[complete.cases(data.plot.res$pct.exp),]# remove empty lines
  
  ### Z Scoring per xaxis
  if (colZ==T) {
    data.plot.res %<>% dplyr::group_by(xaxis) %>%
      dplyr::mutate(z_avg_exp = (avg.exp - mean(avg.exp, na.rm=TRUE)) / sd(avg.exp, na.rm=TRUE)) %>%
      ungroup()
    exp.title = "Scaled expression \n in group"
    fill.values = data.plot.res$z_avg_exp
    ###   
  } else{
    exp.title = "Mean expression \n in group"
    if (median == T) {
      exp.title = "Median expression \n in group"
    }
    fill.values = data.plot.res$avg.exp.log1p
  }
  
  pmain <- ggplot2::ggplot(data.plot.res, ggplot2::aes(x = xaxis,y = id)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::labs(title = unique(data.plot.res$gene))+ ggplot2::coord_fixed(clip = "off") + 
    ggplot2::theme(plot.margin = ggplot2::margin(t = plot.margin[1], 
                                                 r = plot.margin[2],
                                                 b = plot.margin[3],
                                                 l = plot.margin[4], 
                                                 unit = "cm"),
                   axis.text = ggplot2::element_text(color = "black"),
                   legend.direction = "horizontal",
                   axis.text.x = element_text(color = "black",angle = 90,hjust = 1,vjust = 0.5, size = 14, family = "Helvetica"),
                   axis.text.y = element_text(color = "black", size = 14, family = "Helvetica"),
                   axis.title.x = element_text(color = "black", size = 14, family = "Helvetica"),
                   axis.title = element_text(size = 14, color = "black", family = "Helvetica"),
                   plot.title = element_text(size = 14, hjust = 0.5,face="bold", family = "Helvetica"),
                   plot.subtitle = element_text(size = 14, hjust = 0, family = "Helvetica"),
                   axis.line = element_line(color = "black"),
                   strip.text.x = element_text(size = 14, color = "black", family = "Helvetica"),
                   legend.text = element_text(size = 10, color = "black", family = "Helvetica"),
                   legend.title = element_text(size = 10, color = "black", family = "Helvetica", hjust =0),
                   legend.position = "right",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),)
  
  # colorbar.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = exp.title,
  #                                                                  title.position = "top",
  #                                                                  title.hjust = 0.5,
  #                                                                  barwidth = unit(4.5,"cm"),
  #                                                                  frame.colour = "black",
  #                                                                  frame.linewidth = 0.8,
  #                                                                  ticks.colour = "black",
  #                                                                  order = 1))
  # 
  # point.layer <- ggplot2::guides(size = ggplot2::guide_legend(title = "Fraction of cells \n in group (%)", 
  #                                                             title.position = "top", title.hjust = 0.5, label.position = "bottom", 
  #                                                             override.aes = list(color = "black", fill = "grey50"), 
  #                                                             keywidth = ggplot2::unit(0.3, "cm"),
  #                                                             order = 2))
  
  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = exp.title,
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.5,"cm"),
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 2),
                                  size = ggplot2::guide_legend(title = "Fraction of cells \n in group (%)", 
                                                               title.position = "top", title.hjust = 0.5, label.position = "bottom", 
                                                               override.aes = list(color = "black", fill = "grey50"), 
                                                               keywidth = ggplot2::unit(0.3, "cm"),
                                                               order = 1))
  
  dot.col = c("#fff5f0","#990000") # TODO change the scalefillgradient to +n in the else part
  gradient_colors <- c("#fff5f0", "#fcbba1", "#fc9272", "#fb6a4a", "#990000")
  # "#FFFFFF","#08519C","#BDD7E7" ,"#6BAED6", "#3182BD", 
  if (length(dot.col) == 2) {
    breaks <- scales::breaks_extended(n=5)(range(fill.values))
    limits <- c(0,max(range(fill.values)))
    if (max(breaks) > max(limits)) {
      limits[length(limits)] <- breaks[length(breaks)]
    }
    
    pmain <- pmain + ggplot2::scale_fill_gradientn(colours = gradient_colors,
                                                   breaks = breaks,
                                                   #breaks = pretty(as.vector(quantile(fill.values)), n =10),
                                                   limits = limits)
  }else
  {
    pmain <- pmain + ggplot2::scale_fill_gradient2(low = dot.col[1],
                                                   mid = dot.col[2],
                                                   high = dot.col[3],
                                                   midpoint = midpoint)
  }
  pmain <- pmain + 
    ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                     size = pct.exp),
                        color = "black",
                        shape = 21) +
    guides.layer +
    ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
    ggplot2::scale_size_continuous(breaks = setdiff(pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =5),0)[1:5],
                                   limits = c(0,100))
  # pmain <- pmain + guides(size = guide_legend(order=1),
  #                         fill = guide_legend(order=2))
  print(pmain)
  
}


# WIP works with the condition of interest vs "OTHER"
#' @author Mariano Ruz Jurado
#' @title DO MultipleFeature Dotplot
#' @description Creates a logfc dotplot, informations retrieved from a seurat object meta data, works for multiple genes
#' @param SeuratObject The seurat object
#' @param group.by.x x-axis group name to look for in meta data
#' @param group.by.y1 first y-axis group name to look for in meta data
#' @param group.by.y2 second y-axis group name to look for in meta data
#' @param dot.size Vector of dot size 
#' @param plot.margin = plot margins
#' @param Feature Gene of interest
#' @param gene_order vector of desired gene order in plot, by default alphabetically
#' @param colZ IF True calculates the Z-score of the average expression per column, DEPRICATED
#' @param breaks_extended_n value for breaks_extended to aim for legend colourscale steps
#' @param limits_mult value to multiplicate the limits max and min to adjust better for values, if you see grey dots adjust this value
#' @param return dataframe with plot data

DO.MultipleFeature.Dotplot <- function(SeuratObject, Feature, group.by.x = NULL, group.by.y1 = NULL, group.by.y2 = NULL,
                                       dot.size = c(1,6), plot.margin = c(1, 1, 1, 1), wilcox_test=T, colZ=FALSE, gene_order=NULL,
                                       breaks_extended_n = 5, limits_mult = 1.05){
  require(magrittr)
  # Create Feature expression data frame with grouping information
  geneExp <- expm1(Seurat::FetchData(object = SeuratObject, vars = Feature, layer = "data")) # This  might not work for Seuratv5...
  geneExp$id <- paste(SeuratObject@meta.data[[group.by.y1]], " (", 
                      SeuratObject@meta.data[[group.by.y2]], ")", sep = "")
  geneExp$xaxis <- SeuratObject@meta.data[[group.by.x]]
  
  # geneExp %<>% filter(if_any(1, ~ .x > 0)) 
  melt_geneExp <- melt(geneExp)
  
  ListTest <- list()
  uniques_conds <- as.vector(unique(geneExp$id))
  
  cond_non_other_name <- unlist(strsplit(group.by.y1, split = "_"))[1]
  # String manipulation
  other_conds <- uniques_conds[grep("^Other", uniques_conds)]  # Extract "Other" conditions
  non_other_conds <- uniques_conds[!grepl("^Other", uniques_conds)]  # Extract "Non Other" conditions
  # Extract the conditions by removing the "Other" and "Human" prefixes
  cond_other <- sub("Other \\((.*)\\)", "\\1", other_conds)
  cond_non_other <- sub(paste0(cond_non_other_name, " \\((.*)\\)"), "\\1", non_other_conds)
  
  # Match by cell type and combine them into pairs
  for (conds in cond_other) {
    # Find matching "Other" and "Human" conditions for the same cell type
    n_other_conds <- non_other_conds[which(cond_non_other == conds)]
    other_cond <- other_conds[which(cond_other == conds)]
    
    # Add the pair to the list if both conditions are found
    if (length(other_cond) == 1 & length(n_other_conds) == 1) {
      ListTest[[length(ListTest) + 1]] <- c(n_other_conds,other_cond)
    }
  }
  
  
  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE) {
    stat.test <- melt_geneExp %>%
      dplyr::group_by(variable, xaxis) %>%
      rstatix::wilcox_test(value ~ id, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(SeuratObject)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
    stat.test <- stat.test %>%
      mutate(species = str_extract(group1, "(?<=\\().+?(?=\\))"))
    
  }
  
  # Include xaxis in the overall grouping
  data.plot <- lapply(X = unique(geneExp$id), FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, ]
    
    lapply(X = unique(data.use$xaxis), FUN = function(x_axis) {
      data.cell <- data.use[data.use$xaxis == x_axis, 1:(ncol(geneExp) - 2), drop = FALSE]
      avg.exp <- apply(X = data.cell, MARGIN = 2, FUN = function(x) {
        return(expm1(mean(x)))
      })
      pct.exp <- apply(X = data.cell, MARGIN = 2, FUN = PercentAbove, 
                       threshold = 0)
      
      res <- data.frame(id = ident, xaxis = x_axis, avg.exp = avg.exp, pct.exp = pct.exp * 100)
      res$gene <- rownames(res)
      return(res)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .) %>% data.frame()
  
  # get the scale dvalue for plotting 
  data.plot$avg.exp.log1p <- log1p(data.plot$avg.exp) # reapply the log transformation
  
  data.plot.res <- data.plot
  data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = sort(as.vector(unique(data.plot.res$xaxis)), decreasing = T))
  data.plot.res$group <- sapply(strsplit(as.character(data.plot.res$id), 
                                         split = "\\(|\\)"), "[", 2)
  
  # Check if Others is one of the groups to factorize as last item
  if (any(grepl("Other", unique(data.plot.res$id)))) {
    data.plot.res$id <- factor(data.plot.res$id, levels = c(grep("Other", unique(data.plot.res$id), value = T), grep("Other", unique(data.plot.res$id), value = T, invert = T)))
  }
  #data.plot.res$pct.exp <- ifelse(data.plot.res$pct.exp == 0, NA, data.plot.res$pct.exp) # so fraction 0 is not displayed in plot
  #data.plot.res <- data.plot.res[complete.cases(data.plot.res$pct.exp),]# remove empty lines
  
  tmp <- data.plot.res %>% tidyr::separate(id, into = c("condition", "species"), sep = " \\(")
  tmp_diff <-  tmp %>% dplyr::group_by(gene, xaxis, group) %>% dplyr::mutate(diff_pct_exp = pct.exp[condition == unique(tmp[tmp$condition != "Other",]$condition)]-pct.exp[condition == "Other"])
  tmp_diff <-  tmp_diff %>% dplyr::group_by(gene, xaxis, group) %>% dplyr::mutate(fold_change = avg.exp.log1p[condition == unique(tmp[tmp$condition != "Other",]$condition)]-avg.exp.log1p[condition == "Other"])
  
  final_df <- tmp_diff[tmp_diff$condition == unique(tmp[tmp$condition != "Other",]$condition),]
  
  # Perform a left join to match 'p.adj' values based on 'xaxis', 'gene' (final_df), 'variable' (stat.test), and 'group'
  if (wilcox_test==TRUE) {
    final_df <- final_df %>%
      left_join(stat.test, 
                by = c("xaxis" = "xaxis", "gene" = "variable", "group" = "species")) %>%
      mutate(p_adj = p.adj)
    final_df <- final_df %>%
      select(-p.adj, -p, -n1, -n2, -statistic, -.y.)
    
    final_df$sig <- ifelse(as.numeric(final_df$p_adj) > 0.05, ">0.05", "<0.05")
  }
  
  final_df$gene_group <- final_df$gene
  final_df %<>% unite(id, gene, species, sep = " (")
  
  if (is.vector(gene_order)) { # 
    final_df$gene_group <- factor(final_df$gene_group, levels = gene_order)
  }
  
  
  ### Z Scoring per xaxis
  if (colZ==T) {
    data.plot.res %<>% dplyr::group_by(xaxis) %>%
      dplyr::mutate(z_fold_change = (fold_change - mean(fold_change, na.rm=TRUE)) / sd(fold_change, na.rm=TRUE)) %>%
      ungroup()
    exp.title = "Scaled expression \n in group"
    fill.values = data.plot.res$z_avg_exp
    ###   
  } else{
    exp.title = "LogFC expression \n in group"
    fill.values = final_df$fold_change
  }
  pmain <- ggplot2::ggplot(final_df, ggplot2::aes(x = xaxis,y = group)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::labs(title = unique(final_df$condition))+ ggplot2::coord_fixed(clip = "off") + 
    ggplot2::theme(plot.margin = ggplot2::margin(t = plot.margin[1], 
                                                 r = plot.margin[2],
                                                 b = plot.margin[3],
                                                 l = plot.margin[4], 
                                                 unit = "cm"),
                   axis.text = ggplot2::element_text(color = "black"),
                   legend.direction = "horizontal",
                   axis.text.x = element_text(color = "black",angle = 45,hjust = 1, size = 9, family = "Helvetica"),
                   axis.text.y = element_text(color = "black", size = 9, family = "Helvetica"),
                   axis.title.x = element_text(color = "black", size = 10, family = "Helvetica"),
                   axis.title = element_text(size = 10, color = "black", family = "Helvetica"),
                   plot.title = element_text(size = 10, hjust = 0.5,face="bold", family = "Helvetica"),
                   plot.subtitle = element_text(size = 10, hjust = 0, family = "Helvetica"),
                   axis.line = element_line(color = "black"),
                   strip.text.x = element_text(size = 8, color = "black", family = "Helvetica"),
                   legend.text = element_text(size = 8, color = "black", family = "Helvetica"),
                   legend.title = element_text(size = 8, color = "black", family = "Helvetica", hjust =0),
                   legend.position = "right",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),)
  
  # colorbar.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = exp.title,
  #                                                                  title.position = "top",
  #                                                                  title.hjust = 0.5,
  #                                                                  barwidth = unit(4.5,"cm"),
  #                                                                  frame.colour = "black",
  #                                                                  frame.linewidth = 0.8,
  #                                                                  ticks.colour = "black",
  #                                                                  order = 1))
  # 
  # point.layer <- ggplot2::guides(size = ggplot2::guide_legend(title = "Fraction of cells \n in group (%)", 
  #                                                             title.position = "top", title.hjust = 0.5, label.position = "bottom", 
  #                                                             override.aes = list(color = "black", fill = "grey50"), 
  #                                                             keywidth = ggplot2::unit(0.3, "cm"),
  #                                                             order = 2))
  
  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = exp.title,
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.5,"cm"),
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 2),
                                  size = ggplot2::guide_legend(title = "Difference of \n expressing cells (%)", 
                                                               title.position = "top", title.hjust = 0.5, label.position = "bottom", 
                                                               override.aes = list(color = "black", fill = "grey50"), 
                                                               keywidth = ggplot2::unit(0.3, "cm"),
                                                               order = 1),
                                  alpha = ggplot2::guide_legend(title = "Significance",
                                                                title.position = "top",
                                                                title.hjust = 0.5,
                                                                label.position = "bottom",
                                                                keywidth = ggplot2::unit(0.3, "cm"),
                                                                order = 3),
                                  color = ggplot2::guide_legend(title = "Significance",
                                                                title.position = "top",
                                                                title.hjust =0.5,
                                                                keywidth = ggplot2::unit(0, "cm"),
                                                                order = 4,
                                  ))
  
  dot.col = c("#fff5f0","#990000") # TODO change the scalefillgradient to +n in the else part
  gradient_colors <- c("#1b699e","#1f77b4","#fff5f0", "#ea7e1eff","#ea6d1e")
  
  breaks <- scales::breaks_extended(n=breaks_extended_n)(range(fill.values))
  
  if (abs(max(range(fill.values))) > abs(min(range(fill.values)))) {
    limits <- c(-max(range(fill.values)),max(range(fill.values))) 
  } else{
    limits <- c(min(range(fill.values)),-min(range(fill.values))) 
  }
  
  breaks <- unique(sort(c(-breaks,breaks), decreasing = F))[seq(1, length(c(-breaks,breaks)), by=2)]
  breaks <- breaks[!is.na(breaks)]
  # if (!0 %in% breaks) {
  #   breaks <- sort(c(breaks[breaks < 0], 0, breaks[breaks > 0]))
  # }
  if (max(breaks) > max(limits)) {
    limits[length(limits)] <- breaks[length(breaks)]
  }
  if (min(breaks) < min(limits)) {
    limits[1] <- breaks[1]
  }
  
  pmain <- pmain + ggplot2::scale_fill_gradientn(colours = gradient_colors,
                                                 breaks = breaks,
                                                 #breaks = pretty(as.vector(quantile(fill.values)), n =10),
                                                 limits = limits)
  
  if (wilcox_test == F) {
    pmain <- pmain + 
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = diff_pct_exp,
                                       # alpha = ifelse(sig == "<0.05",1,0.75), # deprecated
                                       # color = ifelse(sig == "<0.05", "black", "lightgrey") # must be commented out if run without wilcoxon
      )
      ,
      shape = 21,
      stroke=0.65) +
      guides.layer +
      coord_flip()+
      facet_grid(cols = vars(gene_group), scales = "free")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(final_df$diff_pct_exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(final_df$diff_pct_exp)*1.05,max(final_df$diff_pct_exp)*1.05))+
      ggplot2::scale_color_manual(name = "Significance",
                                  values = c("black" = "black", "lightgrey" = "lightgrey"),
                                  labels= c("black" = ">0.05", "lightgrey" = "<0.05")) +
      theme(panel.spacing = unit(0, "lines"))
  } else{
    
    
    pmain <- pmain + 
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = diff_pct_exp,
                                       # alpha = ifelse(sig == "<0.05",1,0.75), # deprecated
                                       color = ifelse(sig == "<0.05", "black", "lightgrey") # must be commented out if run without wilcoxon
      )
      ,
      shape = 21,
      stroke=0.65) +
      guides.layer +
      coord_flip()+
      # scale_y_discrete(expand = c(1.5,0))+
      facet_grid(cols = vars(gene_group), scales = "free")+
      # ggplot2::scale_alpha_continuous(name = "Significance",
      #                                 breaks = c(1,0.5),
      #                                 limits = c(0,1),
      #                                 labels = c(">0.05", "<0.05")) +
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(final_df$diff_pct_exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(final_df$diff_pct_exp)*1.05,max(final_df$diff_pct_exp)*1.05))+
      ggplot2::scale_color_manual(name = "Significance",
                                  values = c("black" = "black", "lightgrey" = "lightgrey"),
                                  labels= c("black" = ">0.05", "lightgrey" = "<0.05")) +
      theme(panel.spacing = unit(0, "lines"))
    
    # pmain <- pmain + guides(size = guide_legend(order=1),
    #                         fill = guide_legend(order=2))
    
  }
  print(pmain)
  return(pmain$data)
  # pmain <- pmain + theme(axis.ticks.x=element_blank(),axis.text.x = element_text(size = 0)) # remove x-axis and ticks
}

# Boxplot function for one or two given groups per gene, using a pseudo seurat approach 
#' @author Mariano Ruz Jurado
#' @title DO Box plot 
#' @description CReates a box plot similar to Vlnplot function but for an additional condition if wanted, wilcoxon included
#' @param object The seurat object
#' @param group.by group name to look for in meta data
#' @param group.by.2 second group name to look for in meta data
#' @param ctrl.condition select condition to compare to
#' @param outlier_removal Outlier calculation
#' @param vector_colors get the colours for the plot
#' @param wilcox_test If you want to have wilcoxon performed between ctrl.condition and given ones
#' @param stat_pos_mod modificator for where the p-value is plotted increase for higher
#' @param hjust.wilcox value for adjusting height of the text
#' @param vjust.wilcox value for vertical of text
#' @param size.wilcox value for size of text of statistical test
#' @param step_mod value for defining the space between one test and the next one
DO.Box.Plot.wilcox <- function(object,
                               Feature,
                               sample.column = "orig.ident",
                               ListTest=NULL,
                               group.by = "condition",
                               group.by.2 = NULL,
                               ctrl.condition=NULL,
                               outlier_removal = T,
                               vector_colors = c("#1f77b4","#ea7e1eff","royalblue4","tomato2","darkgoldenrod","palegreen4","maroon","thistle3"),
                               wilcox_test = T,
                               stat_pos_mod = 1.15,
                               step_mod = 0,
                               hjust.wilcox=0.5,
                               vjust.wilcox = 0.25,
                               size.wilcox=3.33,
                               hjust.wilcox.2=0.5,
                               vjust.wilcox.2=0,
                               width_errorbar=0.4){
  
  
  #aggregate expression, pseudobulk to visualize the boxplot
  if (is.null(group.by.2)) {
    pseudo_Seu <- AggregateExpression(object,
                                      assays = "RNA",
                                      return.seurat = T,
                                      group.by = c(group.by, sample.column),
                                      verbose = F)
    pseudo_Seu$celltype.con <- pseudo_Seu[[group.by]]
    
  } else{
    pseudo_Seu <- AggregateExpression(object,
                                      assays = "RNA",
                                      return.seurat = T,
                                      group.by = c(group.by, sample.column, group.by.2),
                                      verbose = F)
    pseudo_Seu$celltype.con <- paste(pseudo_Seu[[group.by]][,1], pseudo_Seu[[group.by.2]][,1], sep = "_")
    
  }
  
  
  if (Feature %in% rownames(object)) {
    df_Feature = data.frame(group=setNames(object[[group.by]][,group.by], rownames(object[[group.by]])),
                            Feature = object[["RNA"]]$data[Feature,],
                            cluster = object[[sample.column]])
    df_Feature[,Feature] <- expm1(object@assays$RNA$data[Feature,])
    
  }else{
    df_Feature = data.frame(group=setNames(object[[group.by]][,group.by], rownames(object[[group.by]])),
                            Feature = FetchData(object, vars = Feature)[,1],
                            cluster = object[[sample.column]])
  }
  
  # df_Feature<-data.frame(group=setNames(object[[group.by]][,group.by], rownames(object[[group.by]]))
  # ,orig.ident = object$orig.ident)
  # df_Feature[,Feature] <- expm1(object@assays$RNA$data[Feature,])
  
  #group results and summarize
  if (is.null(group.by.2)) {
    df_melt <- melt(df_Feature) # melt in conditon since the second group might need to get added before the melt
    df_melt_sum <- df_melt %>% 
      dplyr::group_by(group, variable) %>% 
      dplyr::summarise(Mean = mean(value))   
  } else{
    df_Feature[,{group.by.2}] <- setNames(object[[group.by.2]][,group.by.2], rownames(object[[group.by.2]]))
    df_melt <- melt(df_Feature)
    df_melt_sum <- df_melt %>% 
      dplyr::group_by(group, !!sym(group.by.2), variable) %>% #!!sym(), gets the actual variable name useable for dplyr functions
      dplyr::summarise(Mean = mean(value))  
  }
  
  #create comparison list for wilcox, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list 
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with each other")
    group <- unique(object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) { 
      ctrl.condition <- group[grep(pattern = paste(c("CTRL","Ctrl","WT","Wt","wt"),collapse ="|")
                                   ,group)[1]]
    }
    
    
    #create ListTest
    ListTest <- list()
    for (i in 1:length(group)) {
      cndtn <- as.character(group[i]) 
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }
  
  #delete Null values, created by count index also reorder for betetr p-value depiction
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  indices <- sapply(ListTest, function(x) match(x[2], df_melt_sum$group))
  ListTest <- ListTest[order(indices)]
  
  #Function to remove vectors with both elements having a mean of 0 in df.melt.sum, so the testing does not fail
  remove_zeros <- function(lst, df) {
    lst_filtered <- lst
    for (i in seq_along(lst)) {
      elements <- lst[[i]]
      if (all(df[df$group %in% elements, "Mean"] == 0)) {
        lst_filtered <- lst_filtered[-i]
        warning(paste0("Removing Test ", elements[1], " vs ", elements[2], " since both values are 0"))
      }
    }
    return(lst_filtered)
  }
  
  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df_melt_sum)  
  
  if (!is.null(group.by.2) && length(ListTest) > 1) {
    stop("The provided Seurat has more than two groups in group.by and you specified group.by.2, currently not supported (to crowded)!")
  }
  
  #check before test if there are groups in the data which contain only 0 values and therefore let the test fail
  if (is.null(group.by.2)) {
    group_of_zero <- df_melt %>%
      dplyr::group_by(group) %>%
      summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      filter(all_zeros)
    
    if (nrow(group_of_zero) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df_melt <- df_melt %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(value = if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup() 
    }
  } else{
    
    group_of_zero <- df_melt %>%
      dplyr::group_by(group, !!sym(group.by.2)) %>%
      summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      filter(all_zeros)
    
    #check now the result for multiple entries in group.by.2
    groupby2_check <- group_of_zero %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      summarise(group_count = n_distinct(group), .groups = "drop") %>%
      filter(group_count > 1)
    
    if (nrow(groupby2_check) > 0) {
      warning("Some comparisons have no expression in both groups, setting expression to minimum value to ensure test does not fail!")
      df_melt <- df_melt %>%
        dplyr::group_by(group, !!sym(group.by.2)) %>%
        dplyr::mutate(value = if_else(row_number() == 1 & all(value == 0),.Machine$double.xmin,value)) %>%
        ungroup()   
    }
  } 
  
  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE & is.null(group.by.2)) {
    stat.test <- df_melt %>%
      ungroup() %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
  }
  
  #do statistix with rstatix + stats package add second group
  if (wilcox_test == TRUE & !is.null(group.by.2)) {
    stat.test <- df_melt %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      rstatix::wilcox_test(value ~ group, comparisons = ListTest, p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p, method = "bonferroni", n = length(rownames(object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0, sprintf("%.2e",.Machine$double.xmin), sprintf("%.2e", stat.test$p.adj))
  }
  
  
  #pseudobulk boxplot
  theme_box <- function(){
    theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
        axis.line = element_line(colour = "black"),
        #facet_grid colors
        strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
        strip.text = element_text(colour = "black", size = 12),
        # legend.background = element_rect(colour = "grey", fill = "white"),
        # legend.box.background = element_rect(colour = "grey", size = 0.5),
      )
  }
  
  #TODO there might be a problem with the quantiles, recheck
  #SCpubr does not support outlier removal, therefore identify manually
  if (outlier_removal == T) {
    data_matrix <- pseudo_Seu@assays$RNA$data[Feature,]
    for (grp2 in unique(pseudo_Seu[[group.by.2]][,1])) {
      for(grp in unique(pseudo_Seu[[group.by]][,1])){
        group_cells <- pseudo_Seu@meta.data[[group.by.2]] == grp2 & pseudo_Seu@meta.data[[group.by]] == grp
        subset_mat <- data_matrix[group_cells]
        
        Q1 <- quantile(subset_mat, 0.25)
        Q3 <- quantile(subset_mat, 0.75)
        IQR <- Q3 - Q1  # interquartile range calculation
        
        lower_bound <- Q1 - 1.5 * IQR #empirical rule derived from statistics. 1.5 as a default threshold
        upper_bound <- Q3 + 1.5 * IQR
        
        data_matrix_sub <- ifelse(subset_mat >= lower_bound & subset_mat <= upper_bound,
                                  subset_mat,
                                  NA)
        pseudo_Seu@assays$RNA$data[Feature, group_cells] <- data_matrix_sub
      }
    }
  }
  
  
  if (is.null(group.by.2)) {
    p <- SCpubr::do_BoxPlot(sample = pseudo_Seu,
                            feature = Feature,
                            group.by = group.by,
                            order = F,
                            boxplot.width = 0.8,
                            legend.position = "right")
  } else {
    p <- SCpubr::do_BoxPlot(sample = pseudo_Seu,
                            feature = Feature,
                            group.by = group.by.2,
                            split.by = group.by,
                            boxplot.width = 0.8,
                            order = F)
  }
  
  p <- p + geom_point(size=2, alpha=1, position = position_dodge(width = 0.8))+
    scale_fill_manual(values = rep(vector_colors, 2))+ # 16 colours by default, 8 repeat after it
    theme_box()+
    theme(axis.text.x = element_text(color = "black",angle = 45,hjust = 1, size = 16, family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 16, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16,family = "Helvetica",face = "bold"),
          axis.title = element_text(size = 16, color = "black", family = "Helvetica"),
          plot.title = element_text(size = 16, hjust = 0.5, family = "Helvetica"),
          plot.subtitle = element_text(size = 16, hjust = 0, family = "Helvetica"),
          axis.line = element_line(color = "black"),
          strip.text.x = element_text(size = 16, color = "black", family = "Helvetica"),
          legend.text = element_text(size = 14, color = "black", family = "Helvetica"),
          legend.title = element_text(size = 14, color = "black", family = "Helvetica", face = "bold", hjust =0.5),
          legend.position = "bottom")
  #p
  # for only one group
  if (wilcox_test == TRUE & is.null(group.by.2)) {
    
    stat.test_plot <- stat.test %>%
      mutate(y.position = seq(from= max(pseudo_Seu@assays$RNA$data[Feature,][!is.na(pseudo_Seu@assays$RNA$data[Feature,])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))
    # mutate(x.axis = unique(pseudo_Seu$celltype.con)) %>%
    # dplyr::select(x.axis, y.position, p.adj)
    
    p = p + stat_pvalue_manual(stat.test_plot,
                               label = "p = {p.adj}",
                               y.position = "y.position",
                               # step.increase = 0.2,
                               inherit.aes = FALSE,
                               size = size.wilcox,
                               angle= 0,
                               hjust= hjust.wilcox,
                               vjust = vjust.wilcox)
  }
  
  if (wilcox_test == TRUE & !is.null(group.by.2)) {
    
    stat.test_plot <- stat.test %>%
      mutate(y.position = seq(from = max(pseudo_Seu@assays$RNA$data[Feature,][!is.na(pseudo_Seu@assays$RNA$data[Feature,])])*stat_pos_mod, by = step_mod, length.out = nrow(stat.test)))%>%
      mutate(x = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))),
             xmin = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) - 0.2,
             xmax = as.numeric(factor(stat.test[[group.by.2]], levels = unique(stat.test[[group.by.2]]))) + 0.2)
    # dplyr::select(x.axis, y.position, p.adj)
    
    p = p + stat_pvalue_manual(stat.test_plot,
                               label = "p = {p.adj}",
                               y.position = "y.position",
                               # x="x",
                               xmin = "xmin",
                               xmax = "xmax",
                               # xend="xend",
                               # step.increase = 0.2,
                               inherit.aes = FALSE,
                               size = size.wilcox,
                               angle= 0,
                               hjust= hjust.wilcox.2,
                               vjust = vjust.wilcox.2,
                               tip.length = 0.02,
                               bracket.size = 0.8)
  }
  
  print(p)
  
}

# Dotplot function for one or two given groups for multiple genes, using expression values
#' @author Mariano Ruz Jurado
#' @title DO Dot plot
#' @description Creates a dot plot in a similar fashion as the other DO functions but for an additional condition if wanted
#' @param object The seurat object
#' @param group.by.x group name to plot on x-axis
#' @param group.by.y group name to look for in meta data
#' @param group.by.y2 second group name to look for in meta data
#' @param across.group.by.x calculate a pseudobulk expression approach for the x-axis categories
#' @param dot.size Vector of dot size 
#' @param plot.margin = plot margins
#' @param midpoint midpoint in color gradient
#' @param Feature Genes or DF of interest, Data frame should have columns with gene and annotation information, e.g. output of FindAllMarkers
#' @param limits_colorscale Set manually colorscale limits
#' @param colZ IF True calculates the Z-score of the average expression per column
#' @param hide_zero Removes dots for genes with 0 expression
#' @param annotation_x Adds annotation on top of x axis instead on y axis
#' @param point_stroke Defines the thickness of the black stroke on the dots
#' @param ... Further arguments passed to annoSegment function if annotation_x == T
DO.Dotplot <- function(object,
                       Feature, 
                       group.by.x = NULL,
                       group.by.y = NULL,
                       group.by.y2 = NULL,
                       across.group.by.x=F,
                       dot.size = c(1,6),
                       plot.margin = c(1, 1, 1, 1),
                       midpoint = 0.5,
                       colZ=F,
                       returnValue = F,
                       log1p_nUMI=T,
                       hide_zero=T,
                       annotation_x=F,
                       annotation_x_position=0.25,
                       annotation_x_rev=F,
                       point_stroke=0.2,
                       limits_colorscale=NULL,
                       coord_flip=F,
                       ... ){ 
  require(ggtext)
  require(Seurat)
  
  if(!is.vector(Feature) && !is.data.frame(Feature)){
    stop("Feature is not a vector of strings or a data frame!")
  } 
  
  #type of Feature
  FeatureType <- mode(Feature) 
  
  #check if Feature is a vector and annotation specified as true -> no cluster information provided for annotation
  if (is.vector(Feature) && annotation_x == T) {
    stop("Feature is a vector, but annotation_x is set to TRUE. If annotation on xaxis is wanted with specific cluster names you need to provide a dataframe with a column containing cluster names for the genes, like in Seurat::FindAllMarkers!")
  }
  
  #check the input if it is a data frame 
  if (!is.vector(Feature)) {
    orig_DF <- Feature # save original df for annotation purposes
    orig_DF$cluster <- as.vector(orig_DF$cluster)
    cluster_name <- grep("cluster|group|annotation|cell", colnames(Feature), value = T) # Grep name of column, relevant for downstream assignments
    cluster <- unique(Feature[[grep("cluster|group|annotation|cell", colnames(Feature))]])
    Feature <- unique(Feature[[grep("gene|feature", colnames(Feature))]])
    
    if (is.null(cluster) || is.null(Feature)) {
      stop("Couldn't derive Cluster and Feature information from the provided Dataframe. \n Supported names for cluster: cluster|group|annotation|cell\n Supported names for Feature: gene|feature. \n Please make sure that your colnames in the provided Dataframe are supported.")
    }
  }
  
  # Create Feature expression data frame with grouping information
  geneExp <- expm1(FetchData(object = object, vars = Feature, layer = "data")) # 
  
  #catch wrong handling of arguments
  if (is.null(group.by.x) && !is.null(group.by.y) && is.null(group.by.y2)) {
    stop("If you want to make a Marker Plot with just one provided group.by then please use group.by.x!")
  }
  
  geneExp$xaxis <- object@meta.data[[group.by.x]]
  
  if (!is.null(group.by.y) && is.null(group.by.y2)) {
    geneExp$id <- paste(object@meta.data[[group.by.y]], sep = "")
  } else if(!is.null(group.by.y) && !is.null(group.by.y2)){
    geneExp$id <- paste(object@meta.data[[group.by.y]], " (", 
                        object@meta.data[[group.by.y2]], ")", sep = "")    
  } else if(is.null(group.by.y) && is.null(group.by.y2)){
    geneExp$id <- paste(object@meta.data[[group.by.x]], sep = "")
  }
  
  # Include xaxis in the overall grouping
  data.plot <- lapply(X = unique(geneExp$id), FUN = function(ident) {
    data.use <- geneExp[geneExp$id == ident, ]
    
    lapply(X = unique(data.use$xaxis), FUN = function(x_axis) {
      data.cell <- data.use[data.use$xaxis == x_axis, 1:(ncol(geneExp) - 2), drop = FALSE]
      avg.exp <- apply(X = data.cell, MARGIN = 2, FUN = function(x) {
        return(mean(x))
      })
      pct.exp <- apply(X = data.cell, MARGIN = 2, FUN = PercentAbove, 
                       threshold = 0)
      
      res <- data.frame(id = ident, xaxis = x_axis, avg.exp = avg.exp, pct.exp = pct.exp * 100)
      res$gene <- rownames(res)
      return(res)
    }) %>% do.call("rbind", .)
  }) %>% do.call("rbind", .) %>% data.frame()
  
  data.plot.res <- data.plot
  
  #add the cluster information to the plot if annotation_x is set to true and a dataframe was provided with cluster annotation
  #TODO Clean this part a bit up
  if (annotation_x==T && !is.null(cluster)) {
    data.plot.res <- purrr::map_df(seq_len(nrow(data.plot.res)), function(x){
      tmp <- data.plot.res[x,]
      tmp$celltype <- orig_DF[which(orig_DF[[grep("gene|feature", colnames(orig_DF))]] == tmp[[grep("gene|feature", colnames(tmp))]]), cluster_name][[1]]
      return(tmp)
    })
    
    # data.plot.res <- data.plot.res %>%
    #   dplyr::arrange(celltype)
  }
  
  
  data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = sort(unique(data.plot.res$xaxis)))
  
  #create grouping column for multiple grouping variables on the y-axis
  if (!is.null(group.by.y2)) {
    data.plot.res$group <- sapply(strsplit(as.character(data.plot.res$id), 
                                           split = "\\(|\\)"), "[", 2)    
  }
  
  if (hide_zero==T) {
    data.plot.res$pct.exp <- ifelse(data.plot.res$pct.exp == 0, NA, data.plot.res$pct.exp) # so fraction 0 is not displayed in plot
    data.plot.res <- data.plot.res[complete.cases(data.plot.res$pct.exp),]# remove empty lines    
  }
  
  #create bulk expression for group.by.x
  if (across.group.by.x == T) {
    bulk_tmp <- data.plot.res %>%
      dplyr::group_by(id, gene) %>%
      summarise(avg.exp = mean(avg.exp),
                pct.exp = mean(pct.exp))
    bulk_tmp$xaxis <- "Pseudobulk"
    data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
    data.plot.res$xaxis <- factor(data.plot.res$xaxis, levels = c("Pseudobulk", setdiff(sort(unique(data.plot.res$xaxis)), "Pseudobulk")))
  }
  
  # get the scale pvalue for plotting
  if (log1p_nUMI==T) {
    data.plot.res$avg.exp.plot <- log1p(data.plot.res$avg.exp) # reapply the log transformation if wanted
  } else{
    data.plot.res$avg.exp.plot <- data.plot.res$avg.exp
  }
  
  ### TODO Z Scoring per xaxis
  if (colZ==T) {
    data.plot.res %<>% dplyr::group_by(xaxis) %>%
      dplyr::mutate(z_avg_exp = (avg.exp - mean(avg.exp, na.rm=TRUE)) / sd(avg.exp, na.rm=TRUE)) %>%
      ungroup()
    exp.title = "Scaled expression \n in group"
    fill.values = data.plot.res$z_avg_exp
    ###   
  } else if(log1p_nUMI ==T){
    exp.title = "Mean log(nUMI) \n in group"
    fill.values = data.plot.res$avg.exp.plot
  } else{
    exp.title = "Mean nUMI \n in group"
    fill.values = data.plot.res$avg.exp.plot 
  }
  
  #Define which columns to take for dotplot, it should be able to correctly capture one group.by.x
  if (identical(as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis)) && FeatureType=="list") { # go over input type
    # get rid of previous factoring to set new one, first alphabetical order on y
    data.plot.res$xaxis <- as.vector(data.plot.res$xaxis)
    data.plot.res$id <- factor(data.plot.res$id, levels = sort(unique(data.plot.res$id)))
    data.plot.res$gene <- factor(data.plot.res$gene, levels = orig_DF[order(orig_DF$cluster, decreasing = F),]$gene)
    
    if (annotation_x_rev == T) {
      data.plot.res$id <- factor(data.plot.res$id, levels = rev(sort(unique(data.plot.res$id))))
    }
    
    aes_var <- c("gene", "id")
    # there is a second case here for providing just a gene list which need to be adressed with the same aes_var
  } else if(identical(as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis))){
    
    # data.plot.res$id <- factor(data.plot.res$id, levels = sort(unique(data.plot.res$id)))
    # data.plot.res$gene <- factor(data.plot.res$gene, levels = sort(unique(data.plot.res$gene)))
    
    aes_var <- c("gene", "id")
    
  } else{ # all other cases where group.by.y is specified
    data.plot.res$id <- factor(data.plot.res$id, levels = rev(sort(unique(data.plot.res$id))))
    aes_var <- c("xaxis", "id")
    
  }
  
  
  
  pmain <- ggplot2::ggplot(data.plot.res, ggplot2::aes(x = !!sym(aes_var[1]),y = !!sym(aes_var[2]))) + ggplot2::theme_bw(base_size = 14)+ 
    ggplot2::xlab("") + ggplot2::ylab("")+ ggplot2::coord_fixed(clip = "off")+ 
    ggplot2::theme(plot.margin = ggplot2::margin(t = plot.margin[1], 
                                                 r = plot.margin[2],
                                                 b = plot.margin[3],
                                                 l = plot.margin[4], 
                                                 unit = "cm"),
                   axis.text = ggplot2::element_text(color = "black"),
                   legend.direction = "horizontal",
                   axis.text.x = element_text(color = "black",angle = 90,hjust = 1,vjust = 0.5, size = 14, family = "Helvetica"),
                   axis.text.y = element_text(color = "black", size = 14, family = "Helvetica"),
                   axis.title.x = element_text(color = "black", size = 14, family = "Helvetica"),
                   axis.title = element_text(size = 14, color = "black", family = "Helvetica"),
                   plot.title = element_text(size = 14, hjust = 0.5,face="bold", family = "Helvetica"),
                   plot.subtitle = element_text(size = 14, hjust = 0, family = "Helvetica"),
                   axis.line = element_line(color = "black"),
                   strip.text.x = element_text(size = 14, color = "black", family = "Helvetica", face = "bold"),
                   legend.text = element_text(size = 10, color = "black", family = "Helvetica"),
                   legend.title = element_text(size = 10, color = "black", family = "Helvetica", hjust =0),
                   legend.position = "right",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),)
  
  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = exp.title,
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.8,"cm"), # changes the width of the color legend
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 2),
                                  size = ggplot2::guide_legend(title = "Fraction of cells \n in group (%)", 
                                                               title.position = "top", title.hjust = 0.5, label.position = "bottom", 
                                                               override.aes = list(color = "black", fill = "grey50"), 
                                                               keywidth = ggplot2::unit(0.5, "cm"), # changes the width of the precentage dots in legend
                                                               order = 1))
  
  dot.col = c("#fff5f0","#990000") # TODO change the scalefillgradient to +n in the else part
  gradient_colors <- c("#fff5f0", "#fcbba1", "#fc9272", "#fb6a4a", "#990000")
  # "#FFFFFF","#08519C","#BDD7E7" ,"#6BAED6", "#3182BD", 
  if (length(dot.col) == 2) {
    breaks <- scales::breaks_extended(n=5)(range(fill.values))
    
    if (is.null(limits_colorscale)) {
      limits_colorscale <- c(min(range(fill.values))*.99,max(range(fill.values))*1.01)      
    }
    
    if (max(breaks) > max(limits_colorscale)) {
      limits_colorscale[length(limits_colorscale)] <- breaks[length(breaks)]
    }
    
    pmain <- pmain + ggplot2::scale_fill_gradientn(colours = gradient_colors,
                                                   breaks = breaks,
                                                   #breaks = pretty(as.vector(quantile(fill.values)), n =10),
                                                   limits = limits_colorscale)
  }else{
    
    pmain <- pmain + ggplot2::scale_fill_gradient2(low = dot.col[1],
                                                   mid = dot.col[2],
                                                   high = dot.col[3],
                                                   midpoint = midpoint, name = "Gradient")
  }
  if (across.group.by.x == T) {
    
    pmain <- pmain + 
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = pct.exp),shape = 21,stroke=point_stroke)+
      guides.layer +
      facet_grid(cols = vars(gene), scales = "fixed")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(data.plot.res$pct.exp)*1.05,max(data.plot.res$pct.exp)*1.05))+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x=ggtext::element_markdown(color = "black",angle = 90,hjust = 1,vjust = 0.5, size = 14, family = "Helvetica"))+
      scale_x_discrete(labels = function(labels){
        labels <- ifelse(labels== "Pseudobulk", paste0("<b>", labels, "</b>"),labels)
        return(labels)
      }) 
    
  } else if(identical(as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis))){
    
    pmain <- pmain + 
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = pct.exp),shape = 21,stroke=point_stroke)+
      guides.layer +
      # facet_wrap(~facet_group, scales="free_x")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)])*.95,max(data.plot.res$pct.exp)*1.05))+
      theme(panel.spacing = unit(0, "lines"))
    
    if (annotation_x==T) {
      #use the annotation function from jjAnno package
      jjA <- system.file(package = "jjAnno") # Make sure package is installed
      ifelse(nzchar(jjA), "", stop("Install jjAnno CRAN package for annotation on xaxis:!"))
      
      # pmain <- pmain + 
      #   ggplot2::theme(axis.text.y = element_blank(),
      #                  axis.ticks.y = element_blank())
      
      plot_max_y <- ggplot_build(pmain)
      plot_max_y <- plot_max_y$layout$panel_params[[1]]$y.range[2] + annotation_x_position
      
      pmain <- jjAnno::annoSegment(
        object = pmain,
        annoPos = "top",
        aesGroup = T,
        aesGroName = "celltype",
        fontface = "bold",
        fontfamily = "Helvetica",
        pCol = rep("black", length(cluster)),
        textCol = rep("black", length(cluster)),
        addBranch = T,
        branDirection = -1,
        addText = T,
        yPosition = plot_max_y,
        # textSize = 14,
        # hjust = 0.5,
        # vjust = 0,
        # textRot = 0,
        # segWidth = 0.3,
        # lwd = 3
        ...
      )
      
      
    }
    
    if (coord_flip==T) {
      pmain <- pmain + ggplot2::coord_flip()
    }   
    
    if (coord_flip==T && annotation_x==T) {
      warning("Annotation_x and coord_flip set on TRUE, might result in unwanted behaviour!")
    }
    
  } else{
    
    pmain <- pmain + 
      ggplot2::geom_point(ggplot2::aes(fill = fill.values,
                                       size = pct.exp),shape = 21,stroke=point_stroke)+
      guides.layer +
      facet_grid(cols = vars(gene), scales = "fixed")+
      ggplot2::scale_size(range = c(dot.size[1],dot.size[2])) +
      ggplot2::scale_size_continuous(breaks = pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)],
                                     limits = c(min(pretty(round(as.vector(quantile(data.plot.res$pct.exp))), n =10)[seq(1, 10, by = 2)])*.95,max(data.plot.res$pct.exp)*1.05))+
      theme(panel.spacing = unit(0, "lines"))    
  }
  
  if(returnValue == T){
    return(data.plot.res)
  }
  return(pmain)
  
}


DotPlot_costum <- function (object, assay = NULL, features, cols = c("#fff5f0", 
                                                                     "#990000"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                            idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                            scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, cluster_order = NULL) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  if (!is.null(cluster_order)){
    for (i in names(cluster_order)){
      data.plot[as.numeric(as.character(data.plot$id)) %in% cluster_order[[i]],"celltype"] <- i
    }
    data.plot$celltype <- factor(data.plot$celltype, levels = names(cluster_order))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(y = "features.plot", x = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Cluster", y = ifelse(test = is.null(x = split.by), 
                                   yes = "Marker genes", no = "Split Identity")) + theme_cowplot()
  if ( (!is.null(x = feature.groups)) & (!is.null(cluster_order)) ) {
    plot <- plot + facet_grid(vars(feature.groups), vars(celltype), scales = "free", 
                              space = "free") + theme(panel.spacing = unit(x = 1, 
                                                                           units = "lines"), strip.background = element_blank())
  }
  if ( (!is.null(x = feature.groups)) & (is.null(cluster_order)) ) {
    plot <- plot + facet_grid(rows = vars(feature.groups), scales = "free", 
                              space = "free") + theme(panel.spacing = unit(x = 1, 
                                                                           units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient2(low = "#000099", mid = "#fff5f0", high = "#990000")
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}


#' @author David John
#' @title AverageExpressionHeatmap
#' @description Generate a Heatmap with the average expression values
#' @param seuratObj The seurat object
#' @param group.by meta data slot to group by 
#' @param x.axis meta data slot for individUAL COLS IN HEATMAP
#' @param color.group1 color of group 1
#' @param color.group2 color of group 2
#' @param scale sclae parameter forwarded to pheatmap
AverageExpressionHeatmap <- function(seuratObj, genes=c("Cdh5", "Ttn"), group.by="condition", 
                                     x.axis="orig.ident", colors=c(), 
                                     scale="none"){
  require(pheatmap)
  average<-AverageExpression(object = seuratObj, group.by = x.axis, features = genes, assays = "RNA")
  ann<-dplyr::distinct(seuratObj@meta.data[,c(x.axis, group.by)])
  remove_rownames(ann) %>% column_to_rownames(x.axis) -> ann
  if(is.empty(colors)){
    condition <- colorRampPalette(c("blue", "red"))( length(c(unique(ann[[group.by]])))) 
  }
  else{
    condition <- colors
  }
  names(condition) <- c(unique(ann[[group.by]]))
  anno_colors <- list(condition = condition)
  p<-pheatmap(average$RNA, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = ann, annotation_colors = anno_colors, scale = scale)
  return(p)
}

#' @title Generate CellChat object
#' @author David Rodriguez Morales
#' @date 24 - 04 - 2025
#'
#' @description Run inference of cell communication events using cellchat for sc/snRNA-seq. The input can be an AnnData,
#' SingleCellExperiment, Seurat or data matrix with logcounts. If the type is not specified it will be inferred.
#' The input object should contain only one condition.
#'
#' @param obj SingleCellExperiment, AnnData or Seurat Object
#' @param label Metadata Column wit the celltype annotation
#' @param organism human or mouse
#' @param anndata the input object is an AnnData Object
#' @param seurat the input object is a Seurat Object
#' @param sce the input object is a SingleCellExperiment Object
#' @param data.input data matrix with logcounts
#' @param meta dataframe with metadata information
#' @return cellchat object
#'
#' @examples
#' cellchat <- RunCellChat(sce, 'celltypes')
#' cellchat <- RunCellChat(ad, 'celltypes')
#' cellchat <- RunCellChat(seurat, 'celltypes')
#' cellchat <- RunCellChat(data.input = data, meta = meta)
#' @export
GenerateCellChatObject <- function(obj,
                        label,
                        organism = 'human',
                        anndata = F,
                        seurat = F,
                        sce = F,
                        data.input = NULL,
                        seurat_assay = 'RNA',
                        meta = NULL) {
  # Step 0 - Checks
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("The 'CellChat' package is required but not installed. Please install it using devtools::install_github('jinworks/CellChat').")
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required but not installed. Please install it using install.packages('reticulate').")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The 'Seurat' package is required but not installed. Please install it using install.packages('Seurat').")
  }


  if (is.null(obj) && is.null(data.input)) {
    stop('Provide an object or a matrix with logcounts')
  }
  if (!is.null(obj) && !anndata && !seurat && !sce) {  # If an object was provided and user did not to specify the type --> infer
    if (inherits(obj, 'Seurat')) {
      seurat <- T
    } else if (inherits(obj, 'SingleCellExperiment')) {
      sce <- T
    } else if (inherits(obj, "python.builtin.object") && "anndata._core.anndata.AnnData" %in% class(obj)) {
      anndata <- T
    } else {
      stop('Could not infer the type of the input object, specify the type')
    }
  }
  if (organism != 'human' && organism != 'mouse') {  # Check that human or mouse was specified
    stop('Specify the CellChatDB to use: human or mouse')
  }


  # Step 1 - Convert SingleCellExperiment/Anndata to CellChat object
  if (anndata == T) {
    warning('Assuming .X has logcounts')
    data.input <- t(reticulate::py_to_r(obj$X))
    rownames(data.input) <- rownames(reticulate::py_to_r(obj$var))
    colnames(data.input) <- rownames(reticulate::py_to_r(obj$obs))
    meta.data <- reticulate::py_to_r(obj$obs)
    meta <- meta.data
    setAs('dgRMatrix', to = 'dgCMatrix', function(from) { as(as(from, 'CsparseMatrix'), 'dgCMatrix') })
  } else if (sce == T) {
    data.input <- SingleCellExperiment::logcounts(obj)
    meta <- as.data.frame(SingleCellExperiment::colData(obj))
  } else if (seurat == T) {
    data.input <- obj[[seurat_assay]]@data # normalized data matrix
    Seurat::Idents(obj) <- label
    labels <- Seurat::Idents(obj)
    meta <- data.frame(labels = labels, row.names = names(labels))
  } else {
    NULL  # Pass
  }

  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = label)
  cellchat <- CellChat::setIdent(cellchat, ident.use = label)

  # Step 2 - Define the database to use
  if (organism == 'human') {
    CellChatDB.use <- CellChatDB.human
    cellchat@DB <- CellChatDB.use
  } else {
    CellChatDB.use <- CellChatDB.mouse
    cellchat@DB <- CellChatDB.use
  }

  # Step 3 - Cellchat computations
  cellchat <- CellChat::subsetData(cellchat)

  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
	
  tryCatch({
    cellchat <- CellChat::smoothData(cellchat)
  }, error = function(e) {
  # pass (do nothing)
  })


  cellchat <- CellChat::computeCommunProb(cellchat, population.size = T)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)

  cellchat <- CellChat::computeCommunProbPathway(cellchat)

  cellchat <- CellChat::aggregateNet(cellchat)

  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = 'netP')
  return(cellchat)
}



#' @title Extract Significant Interactions from CellChat Object
#' @author David Rodriguez Morales
#' @date 24 - 04 - 2025
#'
#' @description Extract significant cell interactions from a list of cellchat objects. The user can specify the
#' target cell type to extract the interactions for. It will extract the incoming and outgoing interactions
#'
#' @param ccc.list list of cellchat objects
#' @param targets vector with target celltypes to extract the interactions for. (Default: all celltypes)
#' @param annot_col metadata column with celltype annotation
#' @param cond_col metadata column with condition information
#' @param path path to save the excelsheet with the interactions
#' @param filename name of the excel sheet (Default: Significant_Interactions.xlsx)
#' @return None
#'
#' @examples
#' Get_CellChat_netP(cellchat.list, 'all', 'condition', './')
#' @export
Get_CellChat_netP <- function (ccc.list, targets = 'all',  annot_col, cond_col, path, filename = 'Significant_Interactions.xlsx'){

  # Step 0 - Checks
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("The 'CellChat' package is required but not installed. Please install it using devtools::install_github('jinworks/CellChat').")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("The 'openxlsx' package is required but not installed. Please install it using install.packages('openxlsx', dependencies = TRUE).")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("The 'purrr' package is required but not installed. Please install it using   install.packages('purrr').")
  }

  if (!is.list(ccc.list)) {
      message("Input is not a list. Converting to list...")
      ccc.list <- list(ccc.list)
  }

  if (!is.vector(targets)) {
      message("Targets is not a vector. Converting to vector...")
      targets <- c(targets)
   }

  # Step 1 - Get Significant Interactions
  int_all <- list()
  all_cells <-  as.character(unique(cellchat@meta[[annot_col]]))

  for (cellchat in ccc.list) {

    cellchat <- CellChat::rankNetPairwise(cellchat)

    if ('all' %in% targets){
    targets <- as.character(unique(cellchat@meta[[annot_col]]))
    }

    for (target in targets) {
      for (i in all_cells) {  # Incoming to target
        if (i == target){
          next # Skip when target and source is equal
        }

        name <- paste('from', i, 'to', target, 'in', unique(cellchat@meta[[cond_col]]))
        int <- CellChat::identifyEnrichedInteractions(cellchat, from = i, to = target, bidirection = F, pair.only = F)
        if (dim(int)[[1]] == 0) {
          next
        }
        int$pathway_index <- NULL
        int$interaction_name <- NULL
        int$interaction_name_2 <- NULL
        int$name <- name
        int$sender <- i
        int$receiver <- target
        int$condition <- unique(cellchat@meta[[cond_col]])
        # int$value <- cellchat@net$prob[target, i, rownames(int)]
        int_all[[name]] <- int
      }
    }

    for (j in all_cells) {  # Outgoing from target
      if (i == target){
        next  # Skip when target and source is equal
      }
      int <- CellChat::identifyEnrichedInteractions(cellchat, from = target, to = i, bidirection = F, pair.only = F)
      if (dim(int)[[1]]  == 0) {
          next
        }
      name <- paste('from', target, 'to', j, 'in', unique(cellchat@meta[[cond_col]]))
      int$pathway_index <- NULL
      int$interaction_name <- NULL
      int$interaction_name_2 <- NULL
      int$name <- name
      int$sender <- target
      int$receiver <- j
      int$condition <- unique(cellchat@meta[[cond_col]])
      # int$value <- cellchat@net$prob[j, target, rownames(int)]
      int_all[[name]] <- int
    }
  }

  Interaction.db <- purrr::reduce(int_all, full_join)
  wb <- openxlsx::createWorkbook()
  hs1 <- openxlsx::createStyle(fgFill = '#DCE6F1', halign='CENTER', textDecoration='bold', border='Bottom')
  openxlsx::addWorksheet(wb, 'Sig_Interactions')
  openxlsx::writeData(wb, sheet='Sig_Interactions', x=Interaction.db, rowNames=F, headerStyle=hs1)
  openxlsx::saveWorkbook (wb, file=paste0(path,'/', filename),  overwrite=T)
}


