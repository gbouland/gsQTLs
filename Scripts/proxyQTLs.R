## Load util functions
source("./utils.R")
## Required packages
library(reshape2)
## Cell type to run proxy QTL analyses on: Ex, Inh, Astro, Oligo, Micro, OPCs, Endo
cell <- "OPCs"

cells <- c("Ex","Inh","Astro","Oligo","Micro","OPCs","Endo")

micro_GS_QTLs <- lapply(cells,function(cell){
  message(cell)
  ## Pathway / biomolecule database to use: Metabolites, KEGG, microRNA
  database <- "microRNA"
  ## Load expression data, samplesheets and SNP data
  dataList <- prepareData(cell)
  expression <- dataList[[1]]
  samplesheet <- dataList[[2]]
  vars <- dataList[[3]]
  
  ## Load database
  new <- databaseLoader(type = database)
  ## Make proxy expression data
  Proxy_expression <- makePathwayExpression(pathwayInfo = new,
                                            expression = expression,
                                            minGenes = 5,
                                            minExplained = 0.10,
                                            linear = TRUE)#,
  ## Remove duplicates (gene sets are exactly the same, this removes duplicates)
  Proxy_expression <- unique(Proxy_expression)
  dim(Proxy_expression)
  ## Calculate PCs to adjust for unobserved confounding factors
  PCs <- prcomp(t(expression))
  PCs <- PCs$x
  colnames(PCs) <- paste0("mPCs",colnames(PCs))
  samplesheet <- cbind(samplesheet,PCs)
  toCov <- c("age","sex","diagnosis","dataset",colnames(PCs)[1:5])
  ## Calculate proxy QTLs
  proxyQTLs <- lapply(rownames(vars), function(x){
    message(x)
    res <- runLinearModel(expression = Proxy_expression,
                          samplesheet = samplesheet,
                          predictor = x,
                          cov = toCov,
                          padj = "fdr")
    res$variant <- x
    return(res)
  }) |> do.call(what = 'rbind')
  proxyQTLs <- proxyQTLs[order(proxyQTLs$p_adj),]
  return(proxyQTLs)
})

saveRDS(micro_GS_QTLs,"./microRNA_QTLs_age_sex_diagnosis_PC1_PC_5_withMIT.rds")



















