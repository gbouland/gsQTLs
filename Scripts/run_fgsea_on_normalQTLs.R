library(fgsea)
## Load util functions
source("./utils.R")

files <- list.files("./Results/normalQTLs/")

allJustQTLs <- lapply(files, function(x){
  snp <- unlist(strsplit(x,split = "_"))[1]
  celltype <- unlist(strsplit(x,split = "_"))[2]
  out <- readRDS(sprintf("./Results/normalQTLs/%s",x))
  out$variant <- snp
  out$celltype <- celltype
  return(out)
}) |> do.call(what = "rbind")


cells <- c("Ex","Inh","Astro","Oligo","Micro","OPCs","Endo")
fgsea_KEGGAll <- lapply(cells,function(cell){
  message(cell)
  dataList <- prepareData(cell)
  expression <- dataList[[1]]
  samplesheet <- dataList[[2]]
  vars <- dataList[[3]]
  ##Specify database
  type <- "KEGG"
  new <- databaseLoader(type)
  new <- new[new$symbol %in% rownames(expression),]
  
  
  fgsea_KEGG <- lapply(unique(allJustQTLs$variant),function(variant){
    message(variant)
    tmp_results <- allJustQTLs[allJustQTLs$variant == variant & allJustQTLs$celltype == cell,]
    tmp_results <- tmp_results[tmp_results$name %in% new$symbol,]
    genes <- data.frame("gene" = tmp_results$name,"stat" = tmp_results$Estimate)
    ranks <- genes$stat
    names(ranks) <- genes$gene
    ranks <- sort(ranks)
    background <- genes$gene
    newTmp <- new[new$symbol %in% background,]
    pathways <- lapply(unique(newTmp$pathname),function(x){
      newTmp[newTmp$pathname == x,"symbol"]})
    
    names(pathways) <- unique(newTmp$pathname)
    fgseaRes <- fgsea(pathways = pathways, 
                      stats    = ranks,
                      minSize  = 5,
                      maxSize  = 500,
                      nproc=4)
    
    fgseaRes <- as.data.frame(fgseaRes)
    fgseaRes$variant <- variant
    fgseaRes$ID <- paste0(fgseaRes$variant,":",fgseaRes$pathway)
    return(fgseaRes)
  }) |> do.call(what = "rbind")
  return(fgsea_KEGG)
})
saveRDS(fgsea_KEGGAll,file = "Results/KEGG_fgsea.rds")







