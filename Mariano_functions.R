firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#' Volcano plot for visualizing differentially expressed features/genes, works for Human/mice
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features 
#' @param pCut.Off.Value Number for setting the p-adjust threshold
#' @param log2FC.cut.OffValue Number for setting the log2FC threshold
#' @param ylimH Vector containing limitations for y-Axis for Human data
#' @param ylimM Vector containing limitations for y-Axis for Mice data
#' @param xlimH Vector for x.axis for Human data
#' @param xlimM Vector for x.axis for Mice data
#' @param xlimMax If True x.axis is calculated by the max(log2fc) in the provided Data
#' @param labsize Label size of significant Features, set 0 for no labels
volcanoplot <- function(Differential.expression.feature.list,xlimH = 0,xlimM = 0,xlimMax,ylimH,ylimM,pCut.Off.Value,log2FC.cut.OffValue,labsize){
  require(EnhancedVolcano)
  
  for (i in names(Differential.expression.feature.list)) {
    
    i=names(Differential.expression.feature.list[i])
    species=names(Differential.expression.feature.list[i])
    # print(species)
    
    if ("Human" %in% unlist(strsplit(species, split = "_"))) {
      ylim <- ylimH
      xlim <- xlimH
    } 
    
    if ("Mice" %in% unlist(strsplit(species, split = "_"))) {
      ylim <- ylimM
      xlim <- xlimM
    }
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))){
      stop("species not supported")
    }
    
    for (j in names(Differential.expression.feature.list[[species]])) {
      j=names(Differential.expression.feature.list[[species]][j])
      #print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j) 
      contrast.name <- strsplit(marker.contrast, split = "[.]")[[1]][2]
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      if (xlimMax == T){
        xlim <- c(-((max(Differential.expression.feature.list[[species]][marker.contrast][[1]][,2]))-1),((max(Differential.expression.feature.list[[species]][marker.contrast][[1]][,2]))-1))
      }
      pdf(paste0(DEGFolder,"/Volcanoplot_Differential_expression_features_",contrast.name,".pdf"))
      volcano.plot <- EnhancedVolcano(Differential.expression.feature.list[[species]][marker.contrast][[1]],
                                      lab = rownames(Differential.expression.feature.list[[species]][marker.contrast][[1]]),
                                      x = "avg_log2FC",
                                      y = "p_val_adj",
                                      xlim = xlim,
                                      ylim = ylim,
                                      title = contrast.name,
                                      pCutoff = pCut.Off.Value,
                                      FCcutoff = log2FC.cut.OffValue,
                                      pointSize = 1.0,
                                      labSize = labsize,
                                      legendLabels = c('NS', expression(Log[2]~FC), paste("p-adj(",pCut.Off.Value,")",sep = ""), expression(p-adj~and~log[2]~FC)),
                                      subtitle = " ")
      print(volcano.plot)
      dev.off()
      ggsave(filename = paste0(DEGFolder,"/Volcanoplot_Differential_expression_features_",contrast.name,".png"), width = 16, height = 5)
      print(volcano.plot)
    }
  }
}


#' Venn plot visualizing how many differentially expressed features/genes are in a contrast versus another contrast across human/mice WORK IN PROGRESS!!!!
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features 
vennplot <- function(Differential.expression.feature.list.endothelial){
  require(venneuler)
  
  #define species
  species="Human"
  species2="Mice"
  #human gene names
  H.list <- Differential.expression.feature.list.endothelial[grep(species,names(Differential.expression.feature.list.endothelial))]
  #mice gene names
  M.list <- Differential.expression.feature.list.endothelial[grep(species2,names(Differential.expression.feature.list.endothelial))]
  for (i in names(H.list)) {
    
    #get the cell type 
    split.name <- unlist(strsplit(i, split="_"))
    #if cell type might be separated because of the usage of another underscore
    if (length(split.name) >= 4) {
      split.name <- paste0(split.name[2:(length(split.name)-1)],collapse = "_")      
    }
    cell.type <- split.name[split.name %in% unique(SeuratObject.combined.endothelial$cell_type)]
    #check if cell type is found in Seurat object
    if (!(cell.type %in% unique(SeuratObject.combined.endothelial$cell_type))) {
      stop("Cell type not in Seurat Object!")
    } 
    
    
    H.list.celltype <- H.list[grep(cell.type, names(H.list))]
    if (length(H.list.celltype) == 0) {
      cat("No significant genes found in list for",cell.type,"in Human, skipping this venn plot...")
      humanGenes <- NULL
    } else {humanGenes <- rownames(H.list.celltype[[1]][[1]])}
    
    
    M.list.celltype <- M.list[grep(cell.type, names(M.list))]
    if (length(M.list.celltype) == 0) {
      cat("No significant genes found in list for",cell.type,"in Mice, skipping this venn plot...")
      mouseGenes <- NULL
    } else {mouseGenes <- rownames(M.list.celltype[[1]][[1]])}
    
    #contrast names
    marker.contrast1 <- names(H.list.celltype)
    marker.contrast2 <- names(M.list.celltype)
    
    # Calculating dataset size for circles and genes in cross.area are defined
    length.DF.1 <- length(humanGenes)
    # print(length.DF.1)
    
    length.DF.2 <- length(mouseGenes)
    # print(length.DF.2)
    cross.area <- length(humanGenes[humanGenes %in% mouseGenes])
    # print(cross.area)
    
    # venneulers propotions work only with matrix, with colnames as circle names, and rownames should have regulated gene names
    vector.a <- c(humanGenes,mouseGenes) 
    vector.a <- vector.a[!duplicated(vector.a)]
    
    # first put 0s in the rows and coloums, binary for gene in list of sample or not
    mat <- matrix(rep(0,2*length(vector.a)), ncol = 2,dimnames = list(vector.a))
    mat[match(humanGenes,rownames(mat)),2] <- 1      
    mat[match(mouseGenes,rownames(mat)),1] <- 1
    head(mat)        
    dim(mat)    
    if(cross.area <= length.DF.1 && length.DF.1 != 0 && length.DF.2 != 0 && length.DF.1 != length.DF.2 && species != species2){
      matrix.venneuler <- venneuler(mat)
      matrix.venneuler$labels <- ""
      pdf(paste0(outputFolder,"Vennplot_",marker.contrast1,"_",marker.contrast2,".pdf"))
      plot(matrix.venneuler) 
      title(paste0(marker.contrast1," vs ", marker.contrast2))
      text(.15,.55,species)
      text(.15,.50,(length.DF.1-cross.area))
      text(.85,.55,species2)
      text(.85,.50,(length.DF.2-cross.area))
      text(.50,.50,paste0(cross.area))
      dev.off()
    }
  }  
}

# GO analysis over-representation analysis for the subontologies Biological Process, Molecular function and Cellular Component, works for mice and human
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features
#' @param SeuratObject Single Cell Object created by Seurat
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param outputFolder where to save plots
#' @return list with all GO results
GO.Analysis <- function(DEG.list,SeuratObject,pvalueCutoff,pAdjustMethod,qvalueCutoff,minGSSize,maxGSSize,OrgDb.name = "org.Hs.eg.db", outputFolder){
  require(tidyverse)
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  GO.result <- list()
  SeuratObject.combined <- SeuratObject
  for (i in 1:length(names(DEG.list))) {
    
    i=names(DEG.list[i])
    species=names(DEG.list[i])
    print(species)
    SeuratObject.universe <- SeuratObject    
    OrgDb.name <- OrgDb.name
    
    if (OrgDb.name == "org.Hs.eg.db") {
      OrgDb = org.Hs.eg.db      
    }
    
    if (OrgDb.name == "org.Mm.eg.db") {
      OrgDb = org.Mm.eg.db      
    }
    
    
    if (OrgDb.name != "org.Hs.eg.db" && OrgDb.name != "org.Mm.eg.db") {
      stop("species not supported")
    }
    
    
    for (j in names(DEG.list[[species]])) {
      j=names(DEG.list[[species]][j])
      print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      # define universe with ENTREZIDs as data frames
      egallgenes <- as.vector(rownames(SeuratObject.universe)) 
      egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
      head(egallgenes)
      
      # overmit differentially regulated genes with ENTREZIDs as data frame
      egdiffgenes <- as.vector(rownames(DEG.list[[species]][[j]]))
      egdiffgenes = clusterProfiler::bitr(egdiffgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
      head(egdiffgenes)
      write.csv(egdiffgenes, file = paste0(DEGFolder, "/",marker.contrast,"_","DEGs.csv"))
      # enrichGO on all subontologies for differentially regulated genes
      for (m in c("BP", "MF", "CC")) {
        print(m)
        
        #define geneList (ENTREZ IDs and LOG2FC)
        DFgenes <- DEG.list[[species]][[j]]
        
        
        egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))] #get log2FC for DEGs in egdiffgenes
        mygeneList <- egallgenes$stat
        names(mygeneList) <- egallgenes$ENTREZID
        mygeneList<-sort(mygeneList,decreasing = T)
        mygeneList <- mygeneList[!is.na(mygeneList)]
        head(mygeneList)
        length(mygeneList)
        
        # # top 500 from abs values
        # if(length(mygeneList) > 500){
        #   abs.list <- abs(mygeneList)
        #   abs.list <- sort(abs.list, decreasing = T)
        #   gene <- mygeneList[names(abs.list[1:500])]
        # } else{
        #   gene <- mygeneList
        # }
        
        ego <- clusterProfiler::enrichGO(gene = names(mygeneList),
                        OrgDb = OrgDb,
                        keyType = "ENTREZID",
                        ont = m,
                        universe = egallgenes[,2],
                        pvalueCutoff = pvalueCutoff, 
                        pAdjustMethod = pAdjustMethod, 
                        qvalueCutoff = qvalueCutoff, 
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        readable = TRUE,
                        pool = FALSE)
        
        
        ego <- clusterProfiler::setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID")
        openxlsx::write.xlsx(ego@result,file = paste0(DEGFolder,"/GO_",m,".xlsx"))
        GO.result[[species]][[paste0(marker.contrast,"_",m,"_result")]] <- ego@result
        
        
        word.split <- strsplit(ego@result$Description," ")
        
        # head(ego@result$Description)
        head(word.split)
        ### linebreaks for plot, so the name doesnt become too long for the individual terms
        for (b in 1:length(word.split)) 
        {
          
          for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
            if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
            {
              nohyphen <- unlist(strsplit(word.split[[b]], "-"))
              word.split[[b]] <- nohyphen
            }
          }
          
          for (c in 1:length(word.split[[b]])) {  ### checking word length
            if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
            {
              word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
              word.split2 <- word.split2[1:10]
              word.split2 <- paste(word.split2, collapse = "")
              word.split[[b]][c] <- word.split2
            }
          }
          
          
          
          if (length(word.split[[b]]) > 4) ### linebreaks
          {
            linebreaks <- word.split[[b]][1:4]
            linebreaks <- paste(linebreaks, collapse = " ")
            linebreaks2 <- word.split[[b]][5:length(word.split[[b]])] 
            if (length(linebreaks2) > 4) 
            { 
              linebreaks2 <- linebreaks2[1:4]
              
              
            }
            linebreaks2 <- paste(linebreaks2, collapse = " ")
            word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
            
          } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
          
        }
        ego@result$Description <- word.split
        
        GO_File <- ego
        GO_File.Data.Frame <- as.data.frame(GO_File)
        
        # Calculating Gene Ratio so it is numeric and usable for a bubbleplot
        if (length(GO_File.Data.Frame$GeneRatio != 0)) {
          for (i in 1:length(GO_File.Data.Frame$GeneRatio)){
            
            Gene.Ratio <- GO_File.Data.Frame$GeneRatio[i] 
            Gene.Ratio.calculated <- as.numeric(unlist(strsplit(Gene.Ratio, split = "/" ))[1]) / as.numeric(unlist(strsplit(Gene.Ratio, split = "/" ))[2])
            GO_File.Data.Frame$GeneRatio[i] <- Gene.Ratio.calculated
          }
          GO_File.Data.Frame$GeneRatio <- as.numeric(GO_File.Data.Frame$GeneRatio)
          
          
          
          
          
          # Bubble plot for first 10 results or for all results if less 20
          if (length(GO_File.Data.Frame$Description) > 20) {
            Set_Size <- GO_File.Data.Frame[1:20,]
          } else (Set_Size <- GO_File.Data.Frame[1:length(GO_File.Data.Frame$Description),])
          
          #Bubble plot
          pdf(paste0(DEGFolder, "/GO_Bubbleplot_",m,".pdf"))
          ggplot(Set_Size, aes(x = GeneRatio , y = fct_reorder(unlist(Description), GeneRatio ))) +
            geom_point(aes(size = GeneRatio * 455, color = p.adjust)) + 
            labs(size = "Counts") +
            theme_bw(base_size = 14) +
            theme(axis.text.y=element_text(size = 16), axis.text.x = element_text(size = 14))+
            scale_colour_gradient(limits=c(min(GO_File.Data.Frame$p.adjust), max(GO_File.Data.Frame$p.adjust)), low="blue") +
            ylab(NULL) +
            expand_limits(x = c(min(Set_Size$GeneRatio)-(min(Set_Size$GeneRatio)/100*15), max(Set_Size$GeneRatio)+(max(Set_Size$GeneRatio)/100*15))) +
            ggtitle(paste0(m))-> x
          x + guides(size = guide_legend(order = 1)) -> x
          print(x)
          dev.off()
          print(x)
          ggsave(filename = paste0(DEGFolder, "/GO_Bubbleplot_",m,".svg"), width = 10, height = 10)
          # barplot
          GO_File <- ego
          GO_File.Data.Frame <- as.data.frame(GO_File)
          
          if (length(rownames(GO_File.Data.Frame)) < 1) {  ### if data frame is empty
            warning("No significant terms detected")
          } else {
            for (i in 1:length(GO_File.Data.Frame$GeneRatio)){
              GO_File.Data.Frame$GeneRatio[i] <- as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[1]) /
                as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[2])
            }
            GO_File.Data.Frame$GeneRatio <- as.numeric(GO_File.Data.Frame$GeneRatio)
            
            GO_File.Data.Frame2 <- GO_File.Data.Frame
            GO_File.Data.Frame2 <- GO_File.Data.Frame2[order(GO_File.Data.Frame2$p.adjust, decreasing = F),]
            GO_File.Data.Frame <-GO_File.Data.Frame2
            
            GO_File.Data.Frame$Description <- factor(GO_File.Data.Frame$Description, levels = unique(as.character(GO_File.Data.Frame$Description)))
            GO_File.Data.Frame <- transform(GO_File.Data.Frame, Description = reorder(Description, -p.adjust))
            if (length(rownames(GO_File.Data.Frame)) < 20) {k =length(rownames(GO_File.Data.Frame))} else {k=20}
            pdf(paste0(DEGFolder, "/GO_Barplot_",m,".pdf"))
            ggplot(GO_File.Data.Frame[1:k,], aes(reorder(GO_File.Data.Frame, p.adjust), x = Description , y = -log10(p.adjust))) +
              coord_flip() +
              geom_bar(aes( fill = Count), stat="identity") +
              labs(size = "Counts") +
              theme_bw(base_size = 14) +
              theme(axis.text.y=element_text(size = 16)) +
              scale_colour_gradient(limits=c(min(GO_File.Data.Frame$Count), max(GO_File.Data.Frame$Count)), low="blue") +
              ylab("-log10(p.adjusted)") -> x
            x + guides(size = guide_legend(order = 1)) -> x
            print(x)
            dev.off()
            print(x)
            ggsave(filename = paste0(DEGFolder, "/GO_Barplot_",m,".svg"), width = 10, height = 10)
            # gene concept plot
            if (length(rownames(GO_File.Data.Frame)) >= 5) {
              pdf(paste0(DEGFolder, "/GO_GeneConcept_",m,".pdf"))          
              cnet.plot <- cnetplot(ego, foldChange=mygeneList, node_label="all",cex_label_category = 1 , cex_label_gene = 0.6)
              print(cnet.plot)
              dev.off()
              print(cnet.plot)
              ggsave(filename = paste0(DEGFolder, "/GO_GeneConcept_",m,".svg"), width = 16, height = 10)
              # circular gene concept plot
              pdf(paste0(DEGFolder, "/GO_GeneConcept_circular_",m,".pdf")) 
              cnet.plot.circ <- cnetplot(ego, foldChange = mygeneList, circular = TRUE, colorEdge = TRUE,cex_label_category = 1 , cex_label_gene = 0.6, showCategory = 4)
              print(cnet.plot.circ)
              dev.off()
              print(cnet.plot.circ)
              ggsave(filename = paste0(DEGFolder, "/GO_GeneConcept_circular_",m,".svg"), width = 16, height = 10)
              }
            }
          }
        }  
      }
    }
  return(GO.result)  
}

# Go Gene set enrichment analysis using clusterProfiler
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#' @author Mariano Ruz Jurado
#' @param Gene.list A list containing sub-lists with names for species and contrast, containing genes with avg_log2FC (e.g. Gene.list$HFrEF_Human_Cardiomyocytes$Markers$avg_log2FC)
#' @param SeuratObject Single Cell Object created by Seurat
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param outputFolder where to save plots
#' @return list with all GSEA results
GSEA.Analysis <- function(Gene.list,SeuratObject,pvalueCutoff,pAdjustMethod,qvalueCutoff,minGSSize,maxGSSize, OrgDb.name = "org.Hs.eg.db", outputFolder){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(ggnewscale)
  require(enrichplot)
  require(tidyverse)
  require(ggupset)
  GSEA.result <- list()
  for (i in 1:length(names(Gene.list))) {
    
    i=names(Gene.list[i]) # get naming from list for automated plot naming
    species=names(Gene.list[i]) 
    print(species)
    SeuratObject.universe <- SeuratObject    
    OrgDb.name <- OrgDb.name
    
    if (OrgDb.name == "org.Hs.eg.db") {
      OrgDb = org.Hs.eg.db      
    }

    if (OrgDb.name == "org.Mm.eg.db") {
      OrgDb = org.Mm.eg.db      
    }
    
    
    if (OrgDb.name != "org.Hs.eg.db" && OrgDb.name != "org.Mm.eg.db") {
      stop("species not supported")
    }
    
    # define universe with ENTREZIDs as data frames
    egallgenes <- as.vector(rownames(SeuratObject.universe)) 
    egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
    head(egallgenes)    
    
    
    
    for (j in names(Gene.list[[species]])) {
      j=names(Gene.list[[species]][j])
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      DFgenes <- Gene.list[[species]][[j]]
      
      
      egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))] 
      mygeneList <- egallgenes$stat
      names(mygeneList) <- egallgenes$ENTREZID
      mygeneList<-sort(mygeneList,decreasing = T)
      mygeneList <- mygeneList[!is.na(mygeneList)]
      head(mygeneList)
      length(mygeneList)

      # gseGO on all subontologies for differentially regulated genes
      for (m in c("BP", "MF", "CC")) {
        print(m)    
        GO.GSEA <- clusterProfiler::gseGO(geneList = mygeneList,
                         ont = m,
                         OrgDb = OrgDb,
                         keyType = "ENTREZID",
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         verbose = FALSE,
                         by= "fgsea")
        #ridgeplot first,grouped by gene set, density plots are generated by using the frequency of fold change values per gene within each set
        if (length(GO.GSEA@result$ID) >= 1) {
        rdgplot <- ridgeplot.gseaResult(GO.GSEA,showCategory = 20, label_format = 50, core_enrichment = T)+ labs(x="enrichment distribution")
        #rdgplot <- rdgplot + xlim(NA,70)
        print(rdgplot)
        ggsave(filename = paste0(DEGFolder, "/GO_GSEA_ridgeplot_",m,".svg"), width = 12, height = 15)
        dev.off()
        }
        if (length(GO.GSEA@result$ID) >= 1) {
          upsetplot <- clusterProfiler::upsetplot(GO.GSEA,n = 10, label_format = 50) + ylim(-1,2)
          # rdgplot <- rdgplot + xlim(NA,70)
          print(upsetplot)
          ggsave(filename = paste0(DEGFolder, "/GO_GSEA_upsetplot_",m,".svg"), plot = upsetplot, width = 20, height = 7)
          dev.off()
        }
        openxlsx::write.xlsx(GO.GSEA, file = paste0(DEGFolder,"/GO_GSEA",m,".xlsx"))
        GO.GSEA <- clusterProfiler::setReadable(GO.GSEA, OrgDb = OrgDb, keyType = "ENTREZID")
        openxlsx::write.xlsx(GO.GSEA, file = paste0(DEGFolder,"/GO_GSEA_readable",m,".xlsx"))
        GSEA.result[[species]][[paste0(marker.contrast,"_",m,"_result")]] <- GO.GSEA@result
        GOGSEA_File.Data.Frame <- as.data.frame(GO.GSEA)
        if (length(rownames(GOGSEA_File.Data.Frame)) < 1) {  ### if data frame is empty
          print("No significant terms detected")
        } else {
          (GOGSEA_File.Data.Frame[1:20,1:8])
          word.split <- strsplit(GO.GSEA@result$Description," ")
          
          ### linebreaks for plot, so the name doesnt become too long for the individual terms, a little bit outdated
          for (b in 1:length(word.split)) 
          {
            for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
              if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
              {
                nohyphen <- unlist(strsplit(word.split[[b]], "-"))
                word.split[[b]] <- nohyphen
              }
            }
            
            for (c in 1:length(word.split[[b]])) {  ### checking word length
              if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
              {
                word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
                word.split2 <- word.split2[1:10]
                word.split2 <- paste(word.split2, collapse = "")
                word.split[[b]][c] <- word.split2
              }
            }
            
            if (length(word.split[[b]]) > 4)
            {
              linebreaks <- word.split[[b]][1:4]
              linebreaks <- paste(linebreaks, collapse = " ")
              linebreaks2 <- word.split[[b]][5:length(word.split[[b]])]
              if (length(linebreaks2) > 4) 
              { 
                linebreaks2 <- linebreaks2[1:4]
                
                
              }
              linebreaks2 <- paste(linebreaks2, collapse = " ")
              word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
              
            } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
            
          }
          GO.GSEA@result$Description <- word.split
          
          GO_GSEA <- GO.GSEA
          GOGSEA_File.Data.Frame <- as.data.frame(GO_GSEA)
          
          if (length(GO.GSEA@result$ID) > 20) {
            Set_Size <- GOGSEA_File.Data.Frame[1:20,]
          } else (Set_Size <- GOGSEA_File.Data.Frame[1:length(GO.GSEA@result$ID),])
          #bubbleplot
          pdf(paste0(DEGFolder, "/GO_GSEA_",m,".pdf"))
          ggplot(Set_Size, aes(x = NES , y = fct_reorder(unlist(Description), NES ))) +
            geom_point(aes(size = setSize, color = p.adjust)) +
            labs(size = "SetSize") +
            theme_bw(base_size = 14) +
            scale_colour_gradient(limits=c(min(GOGSEA_File.Data.Frame$p.adjust), max(GOGSEA_File.Data.Frame$p.adjust)), low="blue") +
            ylab(NULL) +
            expand_limits(x = c(min(Set_Size$NES)-(min(Set_Size$NES)/100*10), max(Set_Size$NES)+(max(Set_Size$NES)/100*10))) +
            ggtitle(paste0(m))-> plot.GSEA
          plot.GSEA + guides(size = guide_legend(order = 1)) -> plot.GSEA
          
          print(plot.GSEA)        
          dev.off()
          print(plot.GSEA)
          ggsave(filename = paste0(DEGFolder, "/GO_GSEA_",m,".svg"), width = 10, height = 10)
          # barplot
          GO_File <- Set_Size
          GO_File.Data.Frame <- as.data.frame(GO_File)
          
          if (length(rownames(GO_File.Data.Frame)) < 1) {  ### if data frame is empty
            warning("No significant terms detected")
          } else {
            GO_File.Data.Frame$NES <- as.numeric(GO_File.Data.Frame$NES)
            
            GO_File.Data.Frame2 <- GO_File.Data.Frame
            GO_File.Data.Frame2 <- GO_File.Data.Frame2[order(GO_File.Data.Frame2$p.adjust, decreasing = F),]
            GO_File.Data.Frame <-GO_File.Data.Frame2
            
            GO_File.Data.Frame$Description <- factor(GO_File.Data.Frame$Description, levels = unique(as.character(GO_File.Data.Frame$Description)))
            GO_File.Data.Frame <- transform(GO_File.Data.Frame, Description = reorder(Description, -p.adjust))
            if (length(rownames(GO_File.Data.Frame)) < 20) {k =length(rownames(GO_File.Data.Frame))} else {k=20}
            pdf(paste0(DEGFolder, "/GO_Barplot_",m,".pdf"))
            ggplot(GO_File.Data.Frame[1:k,], aes(reorder(GO_File.Data.Frame, p.adjust), x = Description , y = -log10(p.adjust))) +
              coord_flip() +
              geom_bar(aes( fill = setSize), stat="identity") +
              labs(size = "setSize") +
              theme_bw(base_size = 14) +
              theme(axis.text.y=element_text(size = 16)) +
              scale_colour_gradient(limits=c(min(GO_File.Data.Frame$setSize), max(GO_File.Data.Frame$setSize)), low="blue") +
              ylab("-log10(p.adjusted)") -> x
            x + guides(size = guide_legend(order = 1)) -> x
            print(x)
            dev.off()
            print(x)
            ggsave(filename = paste0(DEGFolder, "/GO_GSEA_Barplot_",m,".svg"), width = 10, height = 10)
            # gene concept plot
            if (length(rownames(GO_File.Data.Frame)) >= 5) {
              pdf(paste0(DEGFolder, "/GO_GSEA_GeneConcept_",m,".pdf"))          
              cnet.plot <- clusterProfiler::cnetplot(GO.GSEA, foldChange=mygeneList, node_label="all",cex_label_category = 1 , cex_label_gene = 0.6)
              print(cnet.plot)
              dev.off()
              print(cnet.plot)
              ggsave(filename = paste0(DEGFolder, "/GO_GSEA_GeneConcept_",m,".svg"), width = 16, height = 10)
              # circular gene concept plot
              pdf(paste0(DEGFolder, "/GO_GSEA_GeneConcept_circular_",m,".pdf")) 
              cnet.plot.circ <- clusterProfiler::cnetplot(GO.GSEA, foldChange = mygeneList, circular = TRUE, colorEdge = TRUE,cex_label_category = 1 , cex_label_gene = 0.6, showCategory = 4)
              print(cnet.plot.circ)
              dev.off()
              print(cnet.plot.circ)
              ggsave(filename = paste0(DEGFolder, "/GO_GSEA_GeneConcept_circular_",m,".svg"), width = 16, height = 10)


              
            }
          }
        }  
      }
    }  
  }
 return(GSEA.result) 
}


# KEGG pathway over representation test
#' @author Mariano Ruz Jurado
#' @param DEG.list A list containing sub-lists with names for species and their expression features
#' @param SeuratObject Single Cell Object created by Seurat
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' #' @return list with all KEGG GO results
KEGG.GO.Analysis <- function(DEG.list,SeuratObject,pvalueCutoff,qvalueCutoff,minGSSize,maxGSSize,pAdjustMethod){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(ggnewscale)
  KEGG.result <- list()
  for (i in 1:length(names(DEG.list))) {
    
    i=names(DEG.list[i])
    species=names(DEG.list[i])
    print(species)

    if ("Human" %in% unlist(strsplit(species, split = "_")) && !("integrated" %in% unlist(strsplit(species, split = "_")))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject
      species.KEGG <-"hsa"
    }
    
    if ("Mice" %in% unlist(strsplit(species, split = "_")) && !("integrated" %in% unlist(strsplit(species, split = "_")))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject
      species.KEGG <- "hsa"
    }
    # 
    # if (c("Human","integrated") %in% unlist(strsplit(species, split = "_"))) {
    #   OrgDb.name <- "org.Hs.eg.db"
    #   OrgDb = org.Hs.eg.db
    #   SeuratObject.universe <- SeuratObject
    #   species.KEGG <-"hsa"
    #   
    # }   
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))) {
      stop("species not supported")
    }
    
    for (j in names(DEG.list[[species]])) {
      j=names(DEG.list[[species]][j])
      print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      
      
      # define universe with ENTREZIDs as data frames
      egallgenes <- as.vector(rownames(SeuratObject.universe)) 
      egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
      head(egallgenes)
      
      # overmit differentially regulated genes with ENTREZIDs as data frame
      egdiffgenes <- as.vector(rownames(DEG.list[[species]][[j]]))
      egdiffgenes = clusterProfiler::bitr(egdiffgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
      head(egdiffgenes)
      
      #define geneList (ENTREZ IDs and LOG2FC)
      DFgenes <- DEG.list[[species]][[j]]
      
      
      egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))] #get log2FC for DEGs in egdiffgenes
      mygeneList <- egallgenes$stat
      names(mygeneList) <- egallgenes$ENTREZID
      mygeneList<-sort(mygeneList,decreasing = T)
      mygeneList <- mygeneList[!is.na(mygeneList)]
      head(mygeneList)
      length(mygeneList)
      
      # if(length(mygeneList) > 500){
      #   abs.list <- abs(mygeneList)
      #   abs.list <- sort(abs.list, decreasing = T)
      #   gene <- mygeneList[names(abs.list[1:500])]
      # } else{
      #     gene <- mygeneList
      #     }
      
      enrich.KEGG <- clusterProfiler::enrichKEGG(gene = names(mygeneList),
                                organism = species.KEGG,
                                pvalueCutoff = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize,
                                pAdjustMethod = pAdjustMethod)
      
      
      enrich.KEGG <- clusterProfiler::setReadable(enrich.KEGG, OrgDb = OrgDb, keyType = "ENTREZID")
      KEGG.result[[species]][[paste0(marker.contrast,"_result")]] <- enrich.KEGG@result 
      word.split <- strsplit(enrich.KEGG@result$Description," ")
      openxlsx::write.xlsx(enrich.KEGG,file = paste0(DEGFolder,"/GO_KEGG_",".xlsx"))
      # head(ego@result$Description)
      head(word.split)
      ### linebreaks for plot, so the name doesnt become too long for the individual terms
      for (b in 1:length(word.split)) 
      {
        
        for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
          if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
          {
            nohyphen <- unlist(strsplit(word.split[[b]], "-"))
            word.split[[b]] <- nohyphen
          }
        }
        
        for (c in 1:length(word.split[[b]])) {  ### checking word length
          if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
          {
            word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
            word.split2 <- word.split2[1:10]
            word.split2 <- paste(word.split2, collapse = "")
            word.split[[b]][c] <- word.split2
          }
        }
        
        
        
        if (length(word.split[[b]]) > 4) ### linebreaks
        {
          linebreaks <- word.split[[b]][1:4]
          linebreaks <- paste(linebreaks, collapse = " ")
          linebreaks2 <- word.split[[b]][5:length(word.split[[b]])] 
          if (length(linebreaks2) > 4) 
          { 
            linebreaks2 <- linebreaks2[1:4]
            
            
          }
          linebreaks2 <- paste(linebreaks2, collapse = " ")
          word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
          
        } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
        
      }
      enrich.KEGG@result$Description <- word.split
      GO_File <- enrich.KEGG
      GO_File.Data.Frame <- as.data.frame(GO_File)
      
      if (length(rownames(GO_File.Data.Frame)) < 1) {  ### if data frame is empty
        warning("No significant terms detected")
      } else {
        for (i in 1:length(GO_File.Data.Frame$GeneRatio)){
          GO_File.Data.Frame$GeneRatio[i] <- as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[1]) /
            as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[2])
        }
        GO_File.Data.Frame$GeneRatio <- as.numeric(GO_File.Data.Frame$GeneRatio)
        
        GO_File.Data.Frame2 <- GO_File.Data.Frame
        GO_File.Data.Frame2 <- GO_File.Data.Frame2[order(GO_File.Data.Frame2$p.adjust, decreasing = F),]
        GO_File.Data.Frame <-GO_File.Data.Frame2
        
        GO_File.Data.Frame$Description <- factor(GO_File.Data.Frame$Description, levels = unique(as.character(GO_File.Data.Frame$Description)))
        GO_File.Data.Frame <- transform(GO_File.Data.Frame, Description = reorder(Description, -p.adjust))
        if (length(rownames(GO_File.Data.Frame)) < 20) {k =length(rownames(GO_File.Data.Frame))} else {k=20}
        # barplot
        pdf(paste0(DEGFolder, "/KEGG_GO_Barplot",".pdf"))
        ggplot(GO_File.Data.Frame[1:k,], aes(reorder(GO_File.Data.Frame, p.adjust), x = Description , y = -log10(p.adjust))) +
          coord_flip() +
          geom_bar(aes( fill = Count), stat="identity") +
          labs(size = "Counts") +
          theme_bw(base_size = 14) +
          theme(axis.text.y=element_text(size = 16)) +
          scale_colour_gradient(limits=c(min(GO_File.Data.Frame$Count), max(GO_File.Data.Frame$Count)), low="blue") +
          ylab("-log10(p.adjusted)") -> x
        x + guides(size = guide_legend(order = 1)) -> x
        
        print(x)
        dev.off()
        print(x)
        ggsave(filename = paste0(DEGFolder, "/KEGG_GO_Barplot",".svg"), width = 10, height = 10)
        # gene concept plot 
        pdf(paste0(DEGFolder, "/KEGG_GO_GeneConcept",".pdf"))          
        cnet.plot <- clusterProfiler::cnetplot(enrich.KEGG, foldChange=mygeneList, node_label="all",cex_label_category = 1 , cex_label_gene = 0.6)
        print(cnet.plot)
        dev.off()
        print(cnet.plot)
        ggsave(filename = paste0(DEGFolder, "/KEGG_GO_GeneConcept",".svg"), width = 16, height = 10)
        # circular gene concept plot
        pdf(paste0(DEGFolder, "/KEGG_GO_GeneConcept_circular",".pdf")) 
        cnet.plot.circ <- clusterProfiler::cnetplot(enrich.KEGG, foldChange = mygeneList, circular = TRUE, colorEdge = TRUE,cex_label_category = 1 , cex_label_gene = 0.6, showCategory = 4)
        print(cnet.plot.circ)
        dev.off()
        print(cnet.plot.circ)
        ggsave(filename = paste0(DEGFolder, "/KEGG_GO_GeneConcept_circular",".svg"), width = 16, height = 10)
      }
    }      
  }        
  return(KEGG.result)  
}

# KEGG pathway GSEA test
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
KEGG.GSEA.Analysis <- function(Differential.expression.feature.list,pvalueCutoff,minGSSize,maxGSSize,pAdjustMethod){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(tidyverse)
  for (i in 1:length(names(Differential.expression.feature.list))) {
    
    i=names(Differential.expression.feature.list[i])
    species=names(Differential.expression.feature.list[i])
    print(species)
    
    if ("Human" %in% unlist(strsplit(species, split = "_"))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject.human.combined
      species.KEGG <-"hsa"
    }
    
    if ("Mice" %in% unlist(strsplit(species, split = "_"))) {
      OrgDb.name <- "org.Mm.eg.db"
      OrgDb = org.Mm.eg.db
      SeuratObject.universe <- SeuratObject.mice.combined
      species.KEGG <- "mmu"
    }
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))) {
      stop("species not supported")
    }
    
    # define universe with ENTREZIDs as data frames
    egallgenes <- as.vector(rownames(SeuratObject.universe)) 
    egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
    head(egallgenes)
    
    for (j in names(Differential.expression.feature.list[[species]])) {
      j=names(Differential.expression.feature.list[[species]][j])
      print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      DFgenes <- Differential.expression.feature.list[[species]][[j]]
      
      egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))]
      mygeneList <- egallgenes$stat
      names(mygeneList) <- egallgenes$ENTREZID
      mygeneList<-sort(mygeneList,decreasing = T)
      mygeneList <- mygeneList[!is.na(mygeneList)]
      head(mygeneList)
      length(mygeneList)
      
      # gseKEGG for differentially regulated genes
      gse.KEGG <- clusterProfiler::gseKEGG(geneList  = mygeneList,
                          organism     = species.KEGG,
                          minGSSize    = minGSSize,
                          maxGSSize    = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod= pAdjustMethod,
                          exponent     = 0,
                          verbose      = FALSE)
      
      gse.KEGG <- clusterProfiler::setReadable(gse.KEGG, OrgDb = OrgDb, keyType = "ENTREZID")
      write.table(gse.KEGG, file = paste0(DEGFolder,"/gse_KEGG_",".xls"),sep="\t",quote = F, col.names = NA)
      word.split <- strsplit(gse.KEGG@result$Description," ")
      
      # head(ego@result$Description)
      head(word.split)
      ### linebreaks for plot, so the name doesnt become too long for the individual terms
      for (b in 1:length(word.split)) 
      {
        
        for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
          if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
          {
            nohyphen <- unlist(strsplit(word.split[[b]], "-"))
            word.split[[b]] <- nohyphen
          }
        }
        
        for (c in 1:length(word.split[[b]])) {  ### checking word length
          if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
          {
            word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
            word.split2 <- word.split2[1:10]
            word.split2 <- paste(word.split2, collapse = "")
            word.split[[b]][c] <- word.split2
          }
        }
        
        
        
        if (length(word.split[[b]]) > 4) ### linebreaks
        {
          linebreaks <- word.split[[b]][1:4]
          linebreaks <- paste(linebreaks, collapse = " ")
          linebreaks2 <- word.split[[b]][5:length(word.split[[b]])] 
          if (length(linebreaks2) > 4) 
          { 
            linebreaks2 <- linebreaks2[1:4]
            
            
          }
          linebreaks2 <- paste(linebreaks2, collapse = " ")
          word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
          
        } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
        
      }
      gse.KEGG@result$Description <- word.split
      GSEA_File <- gse.KEGG
      GSEA_File.Data.Frame <- as.data.frame(GSEA_File)
      
      if (length(gse.KEGG@result$ID) > 10) {
        Set_Size <- GSEA_File.Data.Frame[1:10,]
      } else (Set_Size <- GSEA_File.Data.Frame[1:length(gse.KEGG@result$ID),])
      
      pdf(paste0(DEGFolder, "/KEGG_GSEA_",".pdf"))
      ggplot(Set_Size, aes(x = NES , y = fct_reorder(unlist(Description), NES ))) +
        geom_point(aes(size = setSize, color = p.adjust)) + # Adjust * for real ratio!!!
        labs(size = "SetSize") +
        theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(min(GSEA_File.Data.Frame$p.adjust), max(GSEA_File.Data.Frame$p.adjust)), low="blue") +
        ylab(NULL) +
        expand_limits(x = c(min(Set_Size$NES)-(min(Set_Size$NES)/100*10), max(Set_Size$NES)+(max(Set_Size$NES)/100*10))) -> plot.GSEA
      plot.GSEA + guides(size = guide_legend(order = 1)) -> plot.GSEA
      
      print(plot.GSEA)        
      dev.off()
      print(plot.GSEA)
      ggsave(filename = paste0(DEGFolder, "/KEGG_GSEA_",".png"), width = 16, height = 5)
      
      
      
    }
  }
}

# KEGG pathway visualisation, for KEGG GO terms
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features
KEGG.GO.Pathway <- function(DEG.list,SeuratObject){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(xlsx)
  require(pathview)
  for (i in 1:length(names(DEG.list))) {
    
    i=names(DEG.list[i])
    species=names(DEG.list[i])
    print(species)
    SeuratObject.combined <- SeuratObject 
    if ("Human" %in% unlist(strsplit(species, split = "_"))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject.combined
      species.KEGG <-"hsa"
    }
    
    if ("Mice" %in% unlist(strsplit(species, split = "_"))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject.combined
      species.KEGG <- "hsa"
    }
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))) {
      stop("species not supported")
    }
    
    
    # define universe with ENTREZIDs as data frames
    egallgenes <- as.vector(rownames(SeuratObject.universe)) 
    egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
    head(egallgenes)
    
    
    for (j in names(DEG.list[[species]])) {
      j=names(DEG.list[[species]][j])
      print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      #define geneList (ENTREZ IDs and LOG2FC)
      DFgenes <- DEG.list[[species]][[j]]
      
      
      egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))] #get log2FC for DEGs in egdiffgenes
      mygeneList <- egallgenes$stat
      names(mygeneList) <- egallgenes$ENTREZID
      mygeneList<-sort(mygeneList,decreasing = T)
      mygeneList <- mygeneList[!is.na(mygeneList)]
      head(mygeneList)
      length(mygeneList)
      
      # if(length(mygeneList) > 500){
      #   abs.list <- abs(mygeneList)
      #   abs.list <- sort(abs.list, decreasing = T)
      #   gene <- mygeneList[names(abs.list[1:500])]
      # } else{
      #   gene <- mygeneList
      # }
      # 
      if(!dir.exists(paste0(DEGFolder,"/KEGG_Pathways"))) {
        dir.create(paste0(DEGFolder,"/KEGG_Pathways"))
        
      } else{"Pathway Directory exists!"}
      
      
      excelfile <- read.delim(paste0(DEGFolder,"/GO_KEGG_",".xls")) # last column 10
      excelfile.test <- createWorkbook(type = "xlsx")  ### creating excel sheet with hyperlinks to pathways
      sheet <- createSheet(excelfile.test, sheetName = "KEGG_PATHWAYS")
      addDataFrame(excelfile,sheet = sheet, startRow = 1,startColumn = 1)
      rows <- getRows(sheet)
      cells <- createCell(rows,11)
      pb <- txtProgressBar(min = 0, max = length(excelfile$ID), style = 3)     
      
      if (length(excelfile$ID) > 1) {
        for (n in 1:length(excelfile$ID)) {
          Sys.sleep(0.1)
          IDs <- excelfile$ID[n]
          if (IDs!="hsa05163" && IDs!= "hsa05206" && IDs!="hsa04723") {
          suppressMessages(

              
            
            x <- clusterProfiler::pathview(gene.data = names(gene),
                          pathway.id = IDs,
                          species = species.KEGG,
                          limit = list(gene=max(abs(mygeneList)), cpd = 1),
                          #kegg.dir = paste0(DEGFolder,"/KEGG_Pathways"), # this doesnt work...
                          kegg.native = T)
          )
          # move files to KEGG directory
          old_path <- list(paste(c(IDs,".png"),collapse = ""),
                           paste(c(IDs,".xml"), collapse = ""),
                           paste(c(IDs,".pathview.png"), collapse = ""))
          
          old_path <- unlist(old_path)
          new_path <-  list(paste(c(DEGFolder,"/KEGG_Pathways/", IDs,".png"),collapse = ""),
                            paste(c(DEGFolder,"/KEGG_Pathways/", IDs,".xml"), collapse = ""),
                            paste(c(DEGFolder,"/KEGG_Pathways/", IDs,".pathview.png"), collapse = ""))
          
          
          new_path <- unlist(new_path)
          file.copy(old_path, new_path, overwrite = T)
          file.remove(old_path)
          
          cell <- cells[[n+1,1]]
          address <- paste0("file:",DEGFolder,"/KEGG_Pathways/",IDs,".pathview.png")
          setCellValue(cell = cell, paste0("[Pathway]"))
          addHyperlink(cell = cell,address = address)
          saveWorkbook(excelfile.test, paste0(DEGFolder,"/GO_KEGG_Pathways_",".xls"))
          
          # update progress bar
          setTxtProgressBar(pb, n)
          }
        }
      }
      close(pb)}
  }
}


# DO analysis over-representation analysis for Diseases, works for mice and human
#' @author Mariano Ruz Jurado
#' @param Differential.expression.feature.list A list containing sub-lists with names for species and their expression features
#' @param pvalueCutoff p-value Cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing

# many mice genes are lost while transforming into Human gene names resulting in a weaker GeneRatio and therefore in weak p.adjust values, works fine for human
DO.Analysis <- function(Differential.expression.feature.list,pvalueCutoff,qvalueCutoff,minGSSize,maxGSSize,pAdjustMethod){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(ggnewscale)
  require(DOSE)
  require(biomaRt)
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  for (i in 1:length(names(Differential.expression.feature.list))) {
    
    i=names(Differential.expression.feature.list[i])
    species=names(Differential.expression.feature.list[i])
    print(species)
    
    if ("Human" %in% unlist(strsplit(species, split = "_"))) {
      OrgDb.name <- "org.Hs.eg.db"
      OrgDb = org.Hs.eg.db
      SeuratObject.universe <- SeuratObject.human.combined
    }
    
    if (!("Human" %in% unlist(strsplit(species, split = "_"))) && !("Mice" %in% unlist(strsplit(species, split = "_")))) {
      stop("species not supported")
    }
    
    for (j in names(Differential.expression.feature.list[[species]])) {
      j=names(Differential.expression.feature.list[[species]][j])
      print(j)
      marker.contrast <- j
      DEGFolder <- paste0(outputFolder,"/",species,".",j)
      
      if(!dir.exists(DEGFolder)) {
        dir.create(DEGFolder)
      } else{"Directory exists!"}
      
      #convert Mice names to Human names for DO function
      if (species =="Mice") {
        mouseUniverse <- rownames(SeuratObject.mice.combined)    
        genesV1 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouseUniverse , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        SeuratObject.universe <- unique(genesV1[, 2])
        
        mouseGenes <- rownames(Differential.expression.feature.list[[species]][[j]])
        genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouseGenes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        mouseGenes.converted <- unique(genesV2[, 2])
        
        # define universe with ENTREZIDs as data frames
        egallgenes <- SeuratObject.universe 
        egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        head(egallgenes)
        universe <- egallgenes
        
        #build mice genelist with Human Gene names for DO function below
        DFgenes <- Differential.expression.feature.list[[species]][[j]]
        DFgenes.match <- DFgenes[rownames(DFgenes) %in% genesV2[,1],]
        H.names <- genesV2[,2][match(rownames(DFgenes.match),genesV2[,1])]
        rownames(DFgenes.match) <- make.names(H.names,unique = T)
        
        log2FC <- DFgenes.match[["avg_log2FC"]]
        names(log2FC) <- H.names
        
        egdiffgenes = clusterProfiler::bitr(names(log2FC), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        head(egdiffgenes)
        
        egdiffgenes$stat <- DFgenes.match[["avg_log2FC"]][match(egdiffgenes$SYMBOL,rownames(DFgenes.match))]
        mygeneList <- egdiffgenes$stat
        names(mygeneList) <- egdiffgenes$ENTREZID
        mygeneList<-sort(mygeneList,decreasing = T)
        mygeneList <- mygeneList[!is.na(mygeneList)]
        head(mygeneList)
        length(mygeneList)
        
      } else {
        # define universe with ENTREZIDs as data frames
        egallgenes <- as.vector(rownames(SeuratObject.universe)) 
        egallgenes = clusterProfiler::bitr(egallgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb.name)
        head(egallgenes)
        universe <- egallgenes
        
        #define geneList (ENTREZ IDs and LOG2FC)
        DFgenes <- Differential.expression.feature.list[[species]][[j]]
        
        
        egallgenes$stat <- DFgenes[["avg_log2FC"]][match(egallgenes$SYMBOL,rownames(DFgenes))]
        mygeneList <- egallgenes$stat
        names(mygeneList) <- egallgenes$ENTREZID
        mygeneList<-sort(mygeneList,decreasing = T)
        mygeneList <- mygeneList[!is.na(mygeneList)]
        head(mygeneList)
        length(mygeneList)
      }
      genes <- names(mygeneList)[abs(mygeneList) > 1.5] # only genes with log2foldchange greater 1.5
      
      enrich.DO <- clusterProfiler::enrichDO(gene          = genes, 
                            ont           = "DO",
                            pvalueCutoff  = pvalueCutoff,
                            pAdjustMethod = pAdjustMethod,
                            universe      = universe[,2],
                            minGSSize     = minGSSize,
                            maxGSSize     = maxGSSize,
                            qvalueCutoff  = qvalueCutoff,
                            readable      = T)
      
      
      write.table(enrich.DO@result, file = paste0(DEGFolder,"/DO",".xls"),sep="\t",quote = F, col.names = NA)
      word.split <- strsplit(enrich.DO@result$Description," ")
      
      # head(ego@result$Description)
      head(word.split)
      ### linebreaks for plot, so the name doesnt become too long for the individual terms
      for (b in 1:length(word.split)) 
      {
        
        for (d in 1:length(word.split[[b]])) { ### checking for hyphen in words making them too long, erase them
          if (length(unlist(strsplit(word.split[[b]], " "))) < length(unlist(strsplit(word.split[[b]], "-"))))
          {
            nohyphen <- unlist(strsplit(word.split[[b]], "-"))
            word.split[[b]] <- nohyphen
          }
        }
        
        for (c in 1:length(word.split[[b]])) {  ### checking word length
          if (length(unlist(strsplit(word.split[[b]][c], split =""))) > 10)
          {
            word.split2 <- unlist(strsplit(word.split[[b]][c], split =""))
            word.split2 <- word.split2[1:10]
            word.split2 <- paste(word.split2, collapse = "")
            word.split[[b]][c] <- word.split2
          }
        }
        
        
        
        if (length(word.split[[b]]) > 4) ### linebreaks
        {
          linebreaks <- word.split[[b]][1:4]
          linebreaks <- paste(linebreaks, collapse = " ")
          linebreaks2 <- word.split[[b]][5:length(word.split[[b]])] 
          if (length(linebreaks2) > 4) 
          { 
            linebreaks2 <- linebreaks2[1:4]
            
            
          }
          linebreaks2 <- paste(linebreaks2, collapse = " ")
          word.split[[b]] <- paste(linebreaks,linebreaks2, sep = "\n" )
          
        } else (paste(word.split[[b]], collapse = " ") -> word.split[[b]])
        
      }
      enrich.DO@result$Description <- word.split
      GO_File <- enrich.DO
      GO_File.Data.Frame <- as.data.frame(GO_File@result)
      
      if (length(rownames(GO_File.Data.Frame)) < 1 || min(GO_File.Data.Frame$p.adjust) > 0.05) {  ### if data frame is empty
        warning("No significant terms detected")
      } else {
        for (i in 1:length(GO_File.Data.Frame$GeneRatio)){
          GO_File.Data.Frame$GeneRatio[i] <- as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[1]) /
            as.numeric(unlist(strsplit(GO_File.Data.Frame$GeneRatio[i], split = "/" ))[2])
        }
        GO_File.Data.Frame$GeneRatio <- as.numeric(GO_File.Data.Frame$GeneRatio)
        
        GO_File.Data.Frame2 <- GO_File.Data.Frame
        GO_File.Data.Frame2 <- GO_File.Data.Frame2[order(GO_File.Data.Frame2$p.adjust, decreasing = F),]
        GO_File.Data.Frame <-GO_File.Data.Frame2
        
        GO_File.Data.Frame$Description <- factor(GO_File.Data.Frame$Description, levels = unique(as.character(GO_File.Data.Frame$Description)))
        GO_File.Data.Frame <- transform(GO_File.Data.Frame, Description = reorder(Description, -p.adjust))
        if (length(rownames(GO_File.Data.Frame)) < 20) {k =length(rownames(GO_File.Data.Frame))} else {k=20}
        # barplot
        pdf(paste0(DEGFolder, "/DO_Barplot",".pdf"))
        ggplot(GO_File.Data.Frame[1:k,], aes(reorder(GO_File.Data.Frame, p.adjust), x = Description , y = -log10(p.adjust))) +
          coord_flip() +
          geom_bar(aes( fill = Count), stat="identity") +
          labs(size = "Counts") +
          theme_bw(base_size = 14) +
          scale_colour_gradient(limits=c(min(GO_File.Data.Frame$Count), max(GO_File.Data.Frame$Count)), low="blue") +
          ylab("-log10(p.adjusted)") -> x
        x + guides(size = guide_legend(order = 1)) -> x
        
        print(x)
        dev.off()
        print(x)
        ggsave(filename = paste0(DEGFolder, "/DO_Barplot",".png"), width = 16, height = 10)
        # gene concept plot 
        pdf(paste0(DEGFolder, "/DO_GeneConcept",".pdf"))          
        cnet.plot <- clusterProfiler::cnetplot(enrich.DO, foldChange=mygeneList, node_label="all",cex_label_category = 0.8 , cex_label_gene = 0.38)
        print(cnet.plot)
        dev.off()
        print(cnet.plot)
        ggsave(filename = paste0(DEGFolder, "/DO_GeneConcept",".png"), width = 16, height = 10)
        # circular gene concept plot
        pdf(paste0(DEGFolder, "/DO_GeneConcept_circular",".pdf")) 
        cnet.plot.circ <- clusterProfiler::cnetplot(enrich.DO, foldChange = mygeneList, circular = TRUE, colorEdge = TRUE,cex_label_category = 0.8 , cex_label_gene = 0.38)
        print(cnet.plot.circ)
        dev.off()
        print(cnet.plot.circ)
        ggsave(filename = paste0(DEGFolder, "/DO_GeneConcept_circular",".png"), width = 16, height = 10)
      }
    }      
  }        
}




# Automated Annotation of Clusters for Seuratobjects 

#' @author Mariano Ruz Jurado
#' @param SeuratObject SeuratObject to annotate
#' @param ReferenceObject Annotated SeuratObject as Reference, yes it needs to be the same species...

DO.Annotation <- function(SeuratObject,ReferenceObject){
  require(SingleR)
  require(scran)
  require(scater)
  require(ggplotify)  
  #singleR works with SingleCellExperiment like objects, Seurat does have a function for it, also remove possible NAs from Reference
  ref <- as.SingleCellExperiment(ReferenceObject)
  ref <- ref[,!is.na(ref$cell_type)]
  
  SCE.SeuratObject <- as.SingleCellExperiment(SeuratObject)
  
  #Annotation with SingleR
  Cluster.Annotation <- SingleR(test=SCE.SeuratObject, ref=ref, labels=ref$cell_type, de.method="wilcox", clusters = SCE.SeuratObject$seurat_clusters)
  
  #Switch seurat cluster numbers with cell types
  SCE.SeuratObject$cell_type <- factor(
    SCE.SeuratObject$seurat_clusters,
    levels = rownames(Cluster.Annotation),
    labels = Cluster.Annotation$labels)
  
  SeuratObject$cell_type <- SCE.SeuratObject$cell_type
  Idents(SeuratObject) <- "cell_type"
  
  collected <- list()
  all.markers <- metadata(Cluster.Annotation)$de.genes
  empirical.markers <- findMarkers(SCE.SeuratObject, SCE.SeuratObject$cell_type, direction="up")
  
  
  for (i in unique(levels(SeuratObject$cell_type))) {
    
    markers <- unique(unlist(all.markers[[i]]))
    m <- match(markers, rownames(empirical.markers[[i]]))
    m <- markers[rank(m) <= 20] #top 20 upregulated marker genes found in both objects for celltype
    collected[[i]] <- as.grob(plotHeatmap(SCE.SeuratObject, order_columns_by="cell_type", features=m, main =i))
  }
  do.call(ggpubr::ggarrange, c(collected,ncol = 1))
  
  # Dim plots with annotated data
  # p1<-DimPlot(SeuratObject, label = F, group.by = "species")
  # p2<-DimPlot(SeuratObject, label = T)
  # p3<-DimPlot(SeuratObject, label = F, group.by = "orig.ident")
  # names(table(SeuratObject$cell_type))
  # ggarrange.plot <- ggpubr::ggarrange(p1,p2,p3)
  # print(ggarrange.plot)
  return(SeuratObject)
}

#TODO NEEDS MAINTAINANCE, NOT LONGER SUPOORTED -> WE WILL SWITCH TO PACKAGE scDblfinder
# perform Doubletfinder for each sample 
#' @author Mariano Ruz Jurado
#' @param SeuratObject # combined object  
#' @return DoubletFinder analysis results

#DoubletFinder Do not apply DoubletFinder to aggregated scRNA-seq data representing multiple distinct samples (e.g., multiple 10X lanes)
DO.DoubletFinder <- function(SeuratObject){
  require(Doubletfinder)
  
  cat("Performing DoubletFinder...")
  
  #Split object by original identifier
  SeuratObject.list <- SplitObject(SeuratObject, split.by = "orig.ident")
  
  #Samplewise perform Doubletfinder  
  for (i in 1:length(SeuratObject.list)) {
    
    Seu.obj <- SeuratObject.list[[i]]
    homotypic.prop <- modelHomotypic(Seu.obj@meta.data$seurat_clusters)   
    nExp_poi <- round(0.075*nrow(Seu.obj@meta.data))  ## Assuming 7.5% doublet formation rate 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    Seu.obj <- doubletFinder_v3(Seu.obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
    
    #rename Doubletfinders specific names to a more generic one, import for merging objects
    names(Seu.obj@meta.data)[grep("DF.classification",colnames(Seu.obj@meta.data))] <- "DF.classification"
    names(Seu.obj@meta.data)[grep("pANN",colnames(Seu.obj@meta.data))] <- "pANN"
    SeuratObject.list[[i]] <- Seu.obj
  }
  
  #Merging back to one SeuratObject, all previous operations on the object are lost, return only singlet/Doublet analysis information
  SeuratObject <- S4Vectors::merge(SeuratObject.list[[1]], y=c(SeuratObject.list[c(2:length(SeuratObject.list))]))
  DF.classification<- SeuratObject$DF.classification
  
  return(DF.classification)
}

#' #perform SEM Graphs, old look for new function in Bioinformatics sheet
#' #' @author Mariano
#' #' @param SeuratObject # combined object
#' #' @param Features # vector containing featurenames
#' #' @param ListTest # List for which conditions t-test will be performed, if NULL always against provided CTRL 
#' #' @param returnValues # return df.melt.sum data frame containing means and SEM for the set group
#' #' @param ctrl.condition # set your ctrl condition, relevant if running with empty comparison List
#' #' @param group.by # select the seurat object slot where your conditions can be found, default conditon
#' DO.Mean.SEM.Graphs <- function(SeuratObject, Features, ListTest=NULL, returnValues=FALSE, ctrl.condition=NULL, group.by = "condition"){ 
#'   require(ggplot2)
#'   require(ggpubr)
#'   #require(tidyverse)
#'   require(reshape2)
#'   #SEM function definition
#'   SEM <- function(x) sqrt(var(x)/length(x))
#'   #create data frame with conditions from provided SeuratObject, aswell as original identifier of samples
#'   df<-data.frame(condition=setNames(as.character(SeuratObject[[group.by]][,group.by]), rownames(SeuratObject[[group.by]]))
#'                  ,orig.ident = SeuratObject$orig.ident)
#'   #get expression values for genes from individual cells, add to df
#'   for(i in Features){
#'     df[,i] <- SeuratObject@assays$RNA@data[i,]
#'   }
#'   
#'   #melt results 
#'   df.melt <- melt(df)
#'   #group results and summarize, also add/use SEM 
#'   df.melt.sum <- df.melt %>% 
#'     group_by(condition, variable) %>% 
#'     summarise(Mean = mean(value), SEM = SEM(value))
#'   #second dataframe containing mean values for individual samples
#'   df.melt.orig <- df.melt %>% 
#'     group_by(condition, variable, orig.ident) %>% 
#'     summarise(Mean = mean(value))
#'   
#'   #create comparison list for t.test, always against control, so please check your sample ordering
#'   # ,alternative add your own list as argument
#'   if (is.null(ListTest)) {
#'     #if ListTest is empty, so grep the ctrl conditions out of the list 
#'     # and define ListTest comparing every other condition with that ctrl condition
#'     cat("ListTest empty, comparing every sample with provided control")
#'     conditions <- unique(SeuratObject[[group.by]][,group.by])
#'     #set automatically ctrl condition if not provided
#'     if (is.null(ctrl.condition)) { 
#'       ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl","WT","Wt","wt"),collapse ="|")
#'                                         ,conditions)[1]]
#'     }
#'     
#'     df.melt.sum$condition <- factor(df.melt.sum$condition
#'                                     ,levels = c(ctrl.condition,levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
#'     #create ListTest
#'     ListTest <- list()
#'     for (i in 1:length(conditions)) {
#'       cndtn <- conditions[i] 
#'       if(cndtn!=ctrl.condition)
#'       {
#'         ListTest[[i]] <- c(ctrl.condition,cndtn)
#'       }
#'     }
#'   }
#'   #delete Null values, created by count index
#'   ListTest <- ListTest[!sapply(ListTest, is.null)]
#'   #create barplot with significance
#'   p<-ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
#'     geom_col(color = "black")+
#'     geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM), width = .2)+
#'     geom_jitter(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1)+
#'     #ordering, control always first
#'     scale_x_discrete(limits=c(ctrl.condition,levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
#'     #t-test, always against control, using means from orig sample identifier
#'     stat_compare_means(data=df.melt.orig, comparisons = ListTest, method = "t.test", size=3)+
#'     facet_wrap(~variable, ncol = 9, scales = "free") +
#'     scale_fill_manual(values = rep(c("#FFFFFF","#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C"),5)
#'                       , name = "Condition")+
#'     labs( title = "", y = "Mean UMI") +
#'     theme_classic() +
#'     theme(axis.text.x = element_text(face = "bold", color = "black",angle = 45,hjust = 1, size = 10),
#'           axis.text.y = element_text(color = "black"),
#'           axis.title.x = element_blank(),
#'           axis.title = element_text(size = 14, color = "black"),
#'           plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
#'           axis.line = element_line(color = "black", size = .6),
#'           legend.position = "bottom")
#'   print(p)
#'   if (returnValues==TRUE) {
#'     return(df.melt.sum)
#'   }
#' }



###protein sequence matching alignment
#' @author Mariano Ruz Jurado
#' @param mGene # Gene which will be matched against its possible orthologues in replacement
#' @param replacement # vector containing possible orthologues found by biomaRt
#' @param OrthologueList_allHuman # global orthologuelist which contains all Human genes and is filled with mouse orthologues
protein.matching <- function(mGene,replacement,OrthologueList_allHuman){
  
  require(mygene)
  require(UniprotR)
  require(RecordLinkage)
  require(Biostrings)
  
  outh <- try(queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="human"),silent = T)
  outm <- try(queryMany(replacement, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="mouse"),silent = T) # use try for some genes it gives errors
  
  if (!("try-error" %in% class(outm)) && !("try-error" %in% class(outh))) {
    
    #set a TrEMBL prot number or swiss prot number if TREMBL not available
    uniprot.mice <- list()
    for (i in 1:nrow(outm)) {
      if (!is.null(outm[i,]$uniprot.TrEMBL)) { # check if column is avalaible for TREMBL id if not set NA
        uniprot.mice[i] <- outm[i,]$uniprot.TrEMBL
      } else
      {
        uniprot.mice[i] <- NA
      }
      if (is.na(uniprot.mice[i]) && !is.null(outm[i,]$uniprot.Swiss.Prot)) { # check if swiss prot number available
        uniprot.mice[i] <- outm[i,]$uniprot.Swiss.Prot
      }
    }
    
    #set not found entries with NA
    for (i in 1:length(uniprot.mice)) {
      if (!is.null(uniprot.mice[[i]][1]) && !is.na(uniprot.mice[[i]][1])) {
        uniprot.mice[[i]] <- uniprot.mice[[i]][1]
      } else{
        uniprot.mice[[i]] <- NA # set NA to not found in database
        replacement[i] <- NA # set NA to these Gene Symbols as well
        
      }
    }
    #delete NAs
    replacement <- replacement[!is.na(replacement)]
    uniprot.mice <- uniprot.mice[!is.na(uniprot.mice)]
    
    #get sequences based on uniprot IDs
    #human
    sequences.h <- suppressWarnings(GetSequences(unlist(outh$uniprot.Swiss.Prot), directorypath = NULL))
    human.sequences <- sequences.h$Sequence
    if (plyr::empty(sequences.h)) {
      sequences.h <- suppressWarnings(GetSequences(unlist(outh$uniprot.TrEMBL), directorypath = NULL))
      human.sequences <- sequences.h$Sequence
    }
    #mice
    sequences.m <- suppressWarnings(GetSequences(unlist(uniprot.mice), directorypath = NULL))
    mice.sequences <- sequences.m$Sequence
    
    #pairwise alignment for checking sequence similarity, returning a score which is then used for determining best orthologue
    if (!is.null(human.sequences) && !is.null(mice.sequences)) {
      local.Align.list <- list()
      for (k in 1:length(mice.sequences)) {
        localAlign <- pairwiseAlignment(human.sequences,mice.sequences[k], substitutionMatrix="BLOSUM50" , gapOpening = 0, gapExtension = 8) #needleman wunsch
        local.Align.list[[k]] <- localAlign@score # its a global alignment not a local one, just saying
      }
      names(local.Align.list) <- replacement #set replacement names
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))] #order by score
      replacement.hit <- names(local.Align.list[1]) # first entry now has the highest score
      j<-1 # set counter for while
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
        replacement.hit <- names(local.Align.list[j+1])
        j<-j+1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement.hit # set orthologue
      
    }
  }
  else
  {
    human.sequences <- NA
    mice.sequences <- NA
  }
  return(list(OrthologueList_allHuman,human.sequences,mice.sequences))
}



###nucleotide sequence comparison alignment
#' @author Mariano Ruz Jurado
#' @param mGene # Gene which will be matched against its possible orthologues in replacement
#' @param replacement # vector containing possible orthologues found by biomaRt
#' @param OrthologueList_allHuman # global orthologuelist which contains all Human genes and is filled with mouse orthologues
nucleotide.matching <- function(mGene,replacement,OrthologueList_allHuman){
  
  require(rentrez)
  require(mygene)
  require(stringr)
  require(Biostrings)
  out <- try(queryMany(mGene, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="human"),silent = T)
  outm <- try(queryMany(replacement, scopes = "symbol", fields= c("entrezgene", "uniprot"),species ="mouse"),silent = T) # use try for some genes it gives errors
  
  if (!("try-error" %in% class(outm)) && !("try-error" %in% class(out))) {
    #human
    linked_seq_ids <- entrez_link(dbfrom = "gene", id=out$entrezgene, db="nuccore")
    linked_transcripts <- linked_seq_ids$links$gene_nuccore_refseqrna
    head(linked_transcripts)
    if (!is.null(linked_transcripts)) { #check if it hits a sequence human
      
      all_recs <- entrez_fetch(db="nuccore", id=linked_transcripts, rettype = "fasta")
      sequences <- unlist(strsplit(all_recs, split=">"))
      sequences.orth <- sequences[grep("NM",sequences)] # NM is for non predicted mRNAs
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("NR", sequences)] ## check if it is non coding RNA if mRNA is 0
      }
      if (length(sequences.orth) > 1) {
        sequences.orth <- sequences.orth[grep("transcript variant 1", sequences.orth)]
      }
      
      
      sequences.orth <- unlist(strsplit(sequences.orth, split="\n"))
      Non.variant.h <- sequences.orth[2:(length(sequences.orth)-1)]
      Non.variant.h <- paste0(Non.variant.h, collapse="")
    }
    
    #mouse
    if (!is.null(outm$entrezgene)) {
      linked_seq_ids <- entrez_link(dbfrom = "gene", id=outm$entrezgene, db="nuccore")
      linked_transcripts <- linked_seq_ids$links$gene_nuccore_refseqrna
      head(linked_transcripts)
    }
    if (!is.null(linked_transcripts)) { #check if it hits a sequence mice
      
      all_recs <- entrez_fetch(db="nuccore", id=linked_transcripts, rettype = "fasta")
      sequences <- unlist(strsplit(all_recs, split=">"))
      sequences.orth <- sequences[grep("NM",sequences)] # NM is for non predicted mRNAs
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("NR", sequences)] ## check if there is non coding RNA if mRNA is 0
      }
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("XM", sequences)] ## check if there is a predicted coding mRNA
      }
      if (identical(sequences.orth, character(0)) == T) {
        sequences.orth <- sequences[grep("XR", sequences)] ## check if there is a predicted non coding mRNA
      }
      
      orthologues.sequences <- as.list(sequences.orth)
      for (i in 1:length(orthologues.sequences)) {
        orthologues.seq <- unlist(strsplit(orthologues.sequences[[i]], split="\n"))
        orthologues.seq.2 <- orthologues.seq[2:length(orthologues.seq)]
        orthologues.sequences[[i]] <- paste0(orthologues.seq.2, collapse="")
        names(orthologues.sequences)[i] <- orthologues.seq[1]
      }
    }
    if (!is.null(linked_transcripts)) { # no transcript sequences, no alignment
      #alignment
      local.Align.list <- list()
      mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
      for (k in 1:length(orthologues.sequences)) {
        localAlign <- pairwiseAlignment(
          AAString(Non.variant.h),
          AAString(orthologues.sequences[[k]]),
          type="global",
          substitutionMatrix=mat , gapOpening = 5, gapExtension = 2)
        local.Align.list[[k]] <- localAlign@score
        #get gene name from string, stands between "()"
        names(local.Align.list)[k] <- str_match(names(orthologues.sequences)[k],"[(](.*?)[)]")[,2]
      }
      local.Align.list <- local.Align.list[order(-unlist(local.Align.list))] #order by score
      replacement.hit <- names(local.Align.list[1]) # first entry now has the highest score
      j<-1 # set counter for while
      while (replacement.hit %in% OrthologueList_allHuman$MouseGene && !is.na(replacement.hit)) { # check if already in DF, if yes then take second hit
        replacement.hit <- names(local.Align.list[j+1])
        j<-j+1
      }
      OrthologueList_allHuman[OrthologueList_allHuman$HGNC.symbol==mGene,]$MouseGene=replacement.hit # set orthologue
    }
  }
  else
  {
    linked_transcripts <- NA
  }
  return(list(OrthologueList_allHuman, linked_transcripts))
}

###modified simplifyGO function of simplfyEnricment R package, clustering of GO Terms based on GeneIDs in the terms. 
#' @author Mariano Ruz Jurado
simplifyGO2 <- function (mat, method = "binary_cut", control = list(), plot = TRUE, 
                         verbose = TRUE, column_title = qq("@{nrow(mat)} GO terms clustered by '@{method}'"), 
                         ht_list = NULL, ...) {
  require(simplifyEnrichment)
  require(ComplexHeatmap)
  require(grid)
  require(ggplot2)
  require(GetoptLong)
  require(circlize)
  results.collect <- list() 
  if (is.atomic(mat) && !is.matrix(mat)) {
    go_id = mat
    if (!all(grepl("^GO:\\d+$", go_id))) {
      stop_wrap("If you specify a vector, it should contain all valid GO IDs.")
    }
    mat = GO_similarity(go_id)
  }
  cl = do.call(cluster_terms, list(mat = mat, method = method, 
                                   verbose = verbose, control = control))
  go_id = rownames(mat)
  if (is.null(go_id)) {
    go_id = colnames(mat)
  }
  return.df = data.frame(id = go_id, cluster = cl, stringsAsFactors = FALSE)
  if (!all(grepl("^GO:\\d+$", go_id))) {
    stop_wrap("Please ensure GO IDs are the row names of the similarity matrix and should be matched to '^GO:\\\\d+$'.")
  }
  if (plot) {
    heatmap.plot <- ht_clusters(mat, cl, column_title = column_title, ht_list = ht_list,
                        ...)
    

  }
  
  results.collect[["data.frame"]] <- return.df
  results.collect[["plot"]] <- heatmap.plot
  return(results.collect)
}

###modified ridgeplot.gseaResult function of enrichplot R package, clustering of GO Terms based on GeneIDs in the terms. 
#' @author Mariano Ruz Jurado
ridgeplot.gseaResult <- function(x, showCategory=30, fill="p.adjust",
                                 core_enrichment = TRUE, label_format = 30,
                                 orderBy = "NES", decreasing = FALSE) {
  require(DOSE)
  if (!is(x, "gseaResult"))
    stop("currently only support gseaResult")
  
  ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
  if (fill == "qvalue") {
    fill <- "qvalue"
  }
  if (!fill %in% colnames(x@result)) {
    stop("'fill' variable not available ...")
  }
  
  ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')
  if (orderBy !=  'NES' && !orderBy %in% colnames(x@result)) {
    message('wrong orderBy parameter; set to default `orderBy = "NES"`')
    orderBy <- "NES"
  }
  n <- showCategory
  if (core_enrichment) {
    gs2id <- geneInCategory(x)[seq_len(n)]
  } else {
    gs2id <- x@geneSets[x$ID[seq_len(n)]]
  }
  
  if (x@readable) {
    id <- match(names(x@geneList), names(x@gene2Symbol))
    names(x@geneList) <- x@gene2Symbol[id]
  } 
  
  gs2val <- lapply(gs2id, function(id) {
    res <- x@geneList[id]
    res <- res[!is.na(res)]
  })
  
  nn <- names(gs2val)
  i <- match(nn, x$ID)
  nn <- x$Description[i]
  
  # j <- order(x$NES[i], decreasing=FALSE)
  j <- order(x@result[[orderBy]][i], decreasing = decreasing)
  len <- sapply(gs2val, length)
  gs2val.df <- data.frame(category = rep(nn, times=len),
                          color = rep(x[i, fill], times=len),
                          value = unlist(gs2val))
  
  colnames(gs2val.df)[2] <- fill
  gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])
  
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  #exclude every row with only one value for the term, cant be displayed corectly
  gs2val.df <- gs2val.df %>%
    group_by(category) %>%
    filter(n() > 2) %>%
    ungroup()
  # gs2val.df[[fill]] <- -log10(gs2val.df[[fill]])
  
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

  # Define colors with a bit of white in between
  # gradient_colors <- colorRampPalette(c("#4c7f7f","white", "#9a0000"))
  # gradient <- gradient_colors(20)
  
  gs2val.df$fill <- ifelse(gs2val.df[[fill]] > 0.05, NA, gs2val.df[[fill]])
  breaks_extended_n <- 3
  breaks <- scales::breaks_extended(n=breaks_extended_n)(range(gs2val.df$fill))

  ggplot(gs2val.df, aes_string(x="value", y="category")) +
    ggridges::geom_density_ridges(aes(fill=fill)) +
    # enrichplot::set_enrichplot_color(type = "fill", name = fill) + 
    scale_y_discrete(labels = label_func, expand = expansion(mult = c(0.05, 0.15))) +
    geom_vline(xintercept=0, color='royalblue4', linetype='dashed') +
    xlab(NULL) +
    ylab(NULL) +
    theme_box()+
    theme(plot.title = element_text(family="Helvetica",color = "black", hjust = 0.5, size = 14),
          axis.title.y = element_text(family="Helvetica",color = "black", size = 14),
          axis.text.x = element_text(family="Helvetica",color = "black", angle = 0, hjust = 0, size = 14),
          axis.text.y = element_text(family="Helvetica",color = "black", hjust = 1, size = 14),
          legend.position = "top",
          legend.text = element_text(size = 12, color = "black", family = "Helvetica"),
          legend.title = element_text(size = 12, color = "black", family = "Helvetica", face = "bold", hjust =0.5),
          legend.justification = "left",
    )+
    guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                          title.hjust = 0.5,
                                          barwidth = unit(4,"cm"),
                                          barheight = unit(0.5,"cm"),
                                          frame.colour = "black",
                                          frame.linewidth = 0.3,
                                          ticks.colour = "black",
                                          order = 1))+
    # ggplot2::scale_fill_continuous(breaks = pretty(as.vector(quantile(gs2val.df[[fill]])), n = 10)[seq(1, 10, by = 2.5)])+
    scale_fill_gradientn(colors = c("firebrick","mistyrose"),
                         values = scales::rescale(c(0,0.05)),
                         na.value = "#4c7f7f",
                         breaks = signif(breaks,digits = 2),
                         limits=c(0, max(breaks)*1.10), #TODO redefine how the limits are set, considering the actual values in the df
                         name = fill)
    # ggplot2::scale_fill_gradientn(colours = gradient,
    #                               breaks = pretty(as.vector(quantile(gs2val.df[[fill]])), n =10)[seq(1, 10, by = 2)])
}


#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `yulab.utils::str_wrap`
#' @noRd
#' @importFrom yulab.utils str_wrap
default_labeller <- function(n) {
  fun <- function(str){
    str <- gsub("_", " ", str)
    yulab.utils::str_wrap(str, n)
  }
  
  structure(fun, class = "labeller")
}


shorten_term <- function(term) {
  filler_words <- c("of", "in", "the", "on", "for", "and", "to")
  words <- strsplit(term, " ")[[1]]  # Split the term into individual words
  shortened_words <- sapply(words, function(word) {
    if (word %in% filler_words) {
      return("")  # Replace filler words with an empty string
    } else if (nchar(word) > 7) {
      paste0(substr(word, 1, 7), ".")  # Truncate to 7 characters and add a .
    } else {
      word  # Keep the word as is if it's less than or equal to 7 characters
    }
  })
  paste(shortened_words[shortened_words != ""], collapse = " ")  # Join the non-empty words
}

#Function to move the output files of pathview function to directory specified, see pathview documentation for usage
move_pathview <- function(..., path_dir = NULL){
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  wd <- getwd()
  file_copied <- file.copy(paste0(wd, "/", filename), path_dir)
  # img <- png::readPNG(filename)
  # grid::grid.raster(img)
  # if(!save_image) invisible(file.remove(filename))
  # If the copy was successful, remove the original file
  if (file_copied) {
    file.remove(paste0(wd, "/", filename))
  }
}


umap_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
  "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
  "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
)

#' @title Remove Layers from Seurat Object by Pattern
#'
#' @description This function removes layers from a Seurat object's RNA assay based on a specified regular expression pattern.
#' It first backs up the object before removing layers that match the pattern.
#'
#' @param obj A Seurat object.
#' @param pattern A regular expression pattern to match layer names.
#' @param perl A logical value indicating whether to use Perl-compatible regular expressions.
#' Default: `TRUE`.
#'
#' @importFrom CodeAndRoll2 grepv
#' @return A Seurat object with specified layers removed.
#' @export
removeLayersByPattern <- function(obj, pattern = "sc[0-9][0-9]_", perl = TRUE) {
  require(CodeAndRoll2)
  message(paste("pattern: ", pattern))
  stopifnot("obj must be a Seurat object" = inherits(obj, "Seurat"))
  
  layerNames <- Layers(obj)
  layersToRemove <- CodeAndRoll2::grepv(pattern, x = layerNames, perl = perl)
  message(paste(length(layersToRemove), "form", length(layerNames), "layers are removed."))
  obj@assays$RNA@layers[layersToRemove] <- NULL
  return(obj)
}

# Recluster function using FindSubCluster function from Seurat iterative over each cluster -> perfect for fine tuning annotation
#' @author Mariano Ruz Jurado
#' @title DO.FullRecluster
#' @description Creates a refined clustering for each major cluster found initially by FindClusters
#' @param SeuratObject The seurat object
#' @param res Resolution for the new clusters, default 0.5
#' @param algorithm Set one of the available algorithms found in FindSubCLuster function, default = 4: leiden
#' @return Seurat Object with new clustering named eurat_Recluster
DO.FullRecluster <- function(SeuratObject, res = 0.5, algorithm=4, graph.name="RNA_snn"){
  require(progress)
  require(Seurat)
  
  if (is.null(SeuratObject$seurat_clusters)) {
    stop("No seurat clusters defined, please run FindClusters before Reclustering, or fill the slot with a clustering")
  }
  Idents(SeuratObject) <- "seurat_clusters"
  
  SeuratObject$seurat_Recluster <- as.vector(SeuratObject$seurat_clusters)
  pb <- progress_bar$new(total = length(unique(SeuratObject$seurat_clusters)))
  for (cluster in unique(SeuratObject$seurat_clusters)) {
    pb$tick()
    SeuratObject <- FindSubCluster(SeuratObject,
                                   cluster = as.character(cluster),
                                   graph.name = graph.name,
                                   algorithm = algorithm,
                                   resolution = res)  
    
    cluster_cells <- rownames(SeuratObject@meta.data)[SeuratObject$seurat_clusters == cluster]
    SeuratObject$seurat_Recluster[cluster_cells] <- SeuratObject$sub.cluster[cluster_cells]
  }
  SeuratObject$sub.cluster <- NULL
  return(SeuratObject)
}

# Polished UMAP function using Dimplot or FeaturePlot function from Seurat
#' @author Mariano Ruz Jurado
#' @title DO.UMAP
#' @description Creates a polished UMAP from Dimplot or FeaturePlot function
#' @param SeuratObject The seurat object
#' @param FeaturePlot Is it going to be a Dimplot or a FeaturePlot?
#' @param features features for Featureplot
#' @param group.by grouping of plot in DImplot and defines in featureplot the labels
#' @param ... Further arguments passed to DimPlot or FeaturePlot function from Seurat
#' @return Plot with Refined colors and axes
DO.UMAP <- function(SeuratObject, FeaturePlot=F, features=NULL, group.by="seurat_clusters",umap_colors=NULL, text_size=14, label=T,order=T,plot.title=T, legend.position="none", ...){
  require(Seurat)
  
  #Dimplot
  if (FeaturePlot==F) {
    if (is.null(umap_colors)) {
      umap_colors <- rep(c(
        "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
        "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
        "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
      ),5)    
    }
  
    p <- DimPlot(SeuratObject, group.by = group.by, cols = umap_colors, ...) +
      labs(x="UMAP1",y = "UMAP2")+
      theme(plot.title = element_blank(),
            # text = element_text(face = "bold",size = 20),
            axis.title.x = element_text(size = text_size, family="Helvetica"),
            axis.title.y = element_text(size = text_size, family="Helvetica"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = legend.position,
            legend.text = element_text(face = "bold"))
    
    if (label==T) {
      p <-LabelClusters(p, id  = group.by, fontface="bold", box = F)   
    }
    return(p)
  }
  
  #FeaturePlot
  if (FeaturePlot==T) {
    
    if (is.null(features)) {
      stop("Please provide any gene names if using FeaturePlot=T.")
    }
    
    if (is.null(umap_colors)) {
      umap_colors <- c("lightgrey","red2")    
    }

    Idents(SeuratObject) <- group.by
    p <- FeaturePlot(SeuratObject,
                     features = features,
                     cols = umap_colors,
                     label = label,
                     order = order,
                     ...)&
      labs(x="UMAP1",y = "UMAP2")&
      theme(axis.title.x = element_text(size = 14, family="Helvetica"),
            axis.title.y = element_text(size = 14, family="Helvetica"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = legend.position,
            legend.text = element_text(face = "bold"))
    
    if (plot.title == F) {
      p <- p & theme(plot.title = element_blank())
    }
    
    return(p)
  }
}


# Seurat Subset function that actually works
#' @author Mariano Ruz Jurado
#' @title DO.Subset
#' @description Creates a subset from your previous Seurat object with the ident and name you want to subset
#' @param SeuratObject The seurat object
#' @param assay assay to subset by
#' @param ident meta data column to subset for
#' @param ident_name name of group of barcodes in ident of subset for
#' @param ident_thresh numeric thresholds as character, e.g ">5" or c(">5", "<200"), to subset barcodes in ident for
#' @return subsetted Seurat object
DO.Subset <- function(SeuratObject, assay="RNA", ident, ident_name=NULL, ident_thresh=NULL){
  # ident <- "annotation_refined"
  # ident_name <- c("LYVE1+MP","Ccr2+MP")
  
  reduction_names <- names(SeuratObject@reductions)
  SCE_Object <- as.SingleCellExperiment(SeuratObject)
  
  if (!is.null(ident_name) && !is.null(ident_thresh))  {
    stop("Please provide ident_name for subsetting by a name in the column or ident_thresh if it by a threshold")   
  }  
  
  #By a name in the provided column
  if (!is.null(ident_name) && is.null(ident_thresh))  {
    cat("Specified 'ident_name': expecting a categorical variable.")
    SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] %in% ident_name]   
  }
  
  #By a threshold in the provided column
  if (is.null(ident_name) && !is.null(ident_thresh))  {
    cat(paste0("Specified 'ident_thresh': expecting numeric thresholds specified as character, ident_thresh = ", paste0(ident_thresh, collapse = " ")))
    
    #Extract the numeric value and operator
    operator <- gsub("[0-9.]", "", ident_thresh)  
    threshold <- as.numeric(gsub("[^0-9.]", "", ident_thresh))
    
    #operator is valid?
    if ("TRUE" %in% !(operator %in% c("<", ">", "<=", ">="))) {
      stop("Invalid threshold operator provided. Use one of '<', '>', '<=', '>='")
    }
    

    
    #solo case
    if (length(operator) ==1) {
      if (operator == "<") {
        # filtered_cells <- ident_values[ident_values < threshold]
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] < threshold]
      } else if (operator == ">") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] > threshold]
      } else if (operator == "<=") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] <= threshold]
      } else if (operator == ">=") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] >= threshold]
      }
    }
 
    #second case
    if (length(operator) == 2) {
      if (paste(operator,collapse = "") ==  "><") {
        # filtered_cells <- ident_values[ident_values > threshold[1] & ident_values < threshold[2]]
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] > threshold[1] &
                                       SingleCellExperiment::colData(SCE_Object)[, ident] < threshold[2]]
      } else if (paste(operator,collapse = "") ==  "<>") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] < threshold[1] &
                                       SingleCellExperiment::colData(SCE_Object)[, ident] > threshold[2]]
      }
    }   
    
  }
  
  SeuratObject_sub <- as.Seurat(SCE_Object_sub) 
  SeuratObject_sub[[assay]] <- as(object = SeuratObject_sub[[assay]], Class = "Assay5")
  
  #Identify reductions that are not in uppercase -> These are the old ones, not subsetted REMOVE
  # non_uppercase_reductions <- reduction_names[!grepl("^[A-Z0-9_]+$", reduction_names)]
  # SeuratObject_sub@reductions <- SeuratObject_sub@reductions[!names(SeuratObject_sub@reductions) %in% non_uppercase_reductions]
  names(SeuratObject_sub@reductions) <- reduction_names
  SeuratObject_sub$ident <- NULL 
  
  #some checks
  ncells_interest_prior <- nrow(SeuratObject@meta.data[SeuratObject@meta.data[[ident]] %in% ident_name, ])
  ncells_interest_after <- nrow(SeuratObject_sub@meta.data[SeuratObject_sub@meta.data[[ident]] %in% ident_name, ])
  if (ncells_interest_prior != ncells_interest_after) {
    stop(paste0("Number of subsetted cell types is not equal in both objects! Before: ",ncells_interest_prior,"; After: ", ncells_interest_after))
  }
  
  return(SeuratObject_sub)
}

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
