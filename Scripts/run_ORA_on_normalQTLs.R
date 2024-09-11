## Load util functions
source(".utils.R")

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


KEGG_res <- lapply(cells,function(cell){
  message(cell)
  dataList <- prepareData(cell)
  expression <- dataList[[1]]
  samplesheet <- dataList[[2]]
  vars <- dataList[[3]]
  type <- "microRNA"
  new <- databaseLoader(type)
  new <- new[new$symbol %in% rownames(expression),]
  ## subetset only geneset that were calculated with gsQTLs
  genset <- readRDS(file = sprintf("./Results/gsQTLs/%s_QTLs_age_sex_diagnosis_PC1_PC_5_withMIT.rds",type))
  tested <- unique(genset[[cell]]$name)
  new <- new[new$pathname %in% tested,]
  background <- rownames(readRDS(sprintf("./Data/ROSMAP_SEAAD_MIT_%s.rds",cell))[[1]])
  ## run fisher exact
  fisherExactRes <- lapply(unique(allJustQTLs$variant),function(x){
    message(x)
    toTest <- allJustQTLs[allJustQTLs$variant == x & allJustQTLs$celltype == cell & allJustQTLs$p_adj<=0.05,"name"] |> unique()
    print(length(toTest))
    fisherRes <- fisherEnrichment(geneSet = toTest, background = background, database = new)
    fisherRes$variant <- x
    return(fisherRes)
  }) |> do.call(what = "rbind")
  fisherExactRes <- fisherExactRes[order(fisherExactRes$p),]
  rownames(fisherExactRes) <- NULL
  fisherExactRes$ID <- paste0(fisherExactRes$variant,":",fisherExactRes$path)
  return(fisherExactRes)
})
names(KEGG_res) <- cells




