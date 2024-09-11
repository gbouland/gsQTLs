makePathwayExpression <- function(pathwayInfo, expression, minGenes, minExplained, regressOut = F, samplesheet = NULL, cov = NULL, linear = TRUE, sigma = 10){
  require(kernlab)
  pathways <- unique(pathwayInfo$pathname)
  ID_path <- pathwayInfo[pathwayInfo$symbol %in% rownames(expression),]
  if(regressOut){
    genes_test <- unique(ID_path$symbol)
    resids <- lapply(genes_test,function(gene){
      getResiduals(expression,samplesheet,gene,cov)
    }) |> do.call(what = 'rbind')
    rownames(resids) <- genes_test
    colnames(resids) <- colnames(expression)
    expression <- resids
  }
  path_expr <- lapply(pathways,function(x){
    if(match(x,pathways)%%1000 ==0){
      message(sprintf("%s of %s",match(x,pathways),length(pathways)))
    }
    genes <- ID_path[ID_path$pathname == x, "symbol"]
    genes <- unique(genes)
    if(length(genes) <= minGenes){
      return(NA)
    }
    else{
      geneSums <- rowSums(expression[genes,])
      genes <- names(geneSums)[geneSums!=0]
      if(linear){
        PCs <- prcomp(scale(t(expression[genes,])))
        PC1 <- PCs$x[,1]
        rotation <- PCs$rotation[,1]
        eigenValues <- PCs$sdev
      }else{
        PCs <- kpca(t(expression[genes,]),
                    kernel = "rbfdot",
                    kpar = list(sigma = sigma),
                    features = 0)
        PC1 <- PCs@pcv[,1]
        rotation <- cor(PC1,t(expression[genes,]))#PCs@rotated[,1] 
        eigenValues <- PCs@eig
      }
      if(eigenValues[1] / sum(eigenValues) >= minExplained){
        sumLoad <- sum(rotation)
        if(sumLoad<0){
          out <- PC1 * -1
        }else{
          out <- PC1
        }
        return(out)
      }else{
        return(NA)
      }
    }
  })|> do.call(what = "rbind")
  rownames(path_expr) <- pathways
  path_expr <- na.omit(path_expr)
  colnames(path_expr) <- colnames(expression)
  return(path_expr)
}


getResiduals <- function(expr,samplesheet,name,cov){
  tmpData <- samplesheet
  tmpData[,"toResid"] <- expr[name,]
  form <- sprintf("%s ~ %s","toResid",paste0(cov, collapse = " + "))
  resid <- lm(formula = form, data = tmpData)$residuals
  return(unname(resid))
}

runLinearModel <- function(expression, samplesheet, predictor, cov = NULL, padj, type = "linear"){
  outcomes <- rownames(expression)
  rownames(expression) <- paste0(make.names(rownames(expression)),"_expr")
  tmpData <- cbind(samplesheet,t(expression))
  tmpData <- tmpData[!is.na(tmpData[,predictor]),]
  res <- lapply(rownames(expression),function(x){
    #if(match(x,rownames(expression))%%500 ==0){message(match(x,rownames(expression)))}
    
    if(!is.null((cov))){
      form <- sprintf("%s ~ %s + %s",x,predictor, paste0(cov, collapse = " + "))
    }
    else{
      form <- sprintf("%s ~ %s",x,predictor)
    }
    
    if(type == "LR"){
      form <- sprintf("%s ~ %s + %s",predictor,x, paste0(cov, collapse = " + "))
      res <- glm(formula = form, data = tmpData,family = "binomial")
      summ <- summary(res)
    }else{
      res <- lm(formula = form,data = tmpData)
      summ <- summary(res)
    }
    out <- data.frame(summ$coefficients[2,1],
                      summ$coefficients[2,2],
                      summ$coefficients[2,3],
                      summ$coefficients[2,4])
    
    colnames(out) <- switch(type,
                            "LR" = c("Estimate","stdError","z","p"),
                            "linear" = c("Estimate","stdError","t","p"))
    
    
    return(out)
  }) |> do.call(what = "rbind")
  
  res$name <- outcomes
  res$p_adj <- p.adjust(res$p,method = padj)
  res <- res[order(res$p),]
  return(res)
}

fisherEnrichment <- function(geneSet, background, database){
  background <- intersect(background,unique(database$symbol))
  database <- database[database$symbol %in% background,]
  set_to_test <- geneSet[geneSet %in% background]
  pathways <- lapply(unique(database$pathname),function(x){
    genes <- database[database$pathname == x,"symbol"]
    out <- (background %in% genes)*1
    return(out)
  }) |> do.call(what = "cbind")
  colnames(pathways) <- unique(database$pathname)
  rownames(pathways) <- background
  set_to_test <- (background %in% set_to_test)*1
  
  enrich_res <- apply(pathways,2,function(x){
    tab <- table(x,set_to_test)
    if(prod(dim(tab)) != 4){
      out <- out <- data.frame("OR" = NA,
                               "p" = NA,
                               "pathwaySize" = NA,
                               "geneSetSize" = NA,
                               "intersect" = NA)
    }else{
      res <- fisher.test(tab)
      out <- data.frame("OR" = res$estimate,
                        "p" = res$p.value,
                        "pathwaySize" = sum(x),
                        "geneSetSize" = sum(set_to_test),
                        "intersect" = sum(x*set_to_test))
    }
    return(out)
  }) |> do.call(what = "rbind")
  enrich_res$path <- rownames(enrich_res)
  enrich_res <- enrich_res[order(enrich_res$p),]
  enrich_res$fdr <- p.adjust(enrich_res$p,method = "BH")
  return(enrich_res)
}

metaFisher <- function(results, database, expression){
  require(poolr)
  database <- database[database$symbol %in% results$name,]
  fisherPs <- sapply(unique(database$pathname),function(x){
    genesIn <- database[database$pathname == x,"symbol"]
    genesIn <- genesIn[genesIn %in% results$name]
    dependence <- cor(t(expression[genesIn,]))
    if(is.na(dependence[1,1])){
      out <- NA
    }else{
      Ps <- results[match(genesIn,results$name),"p"]
      res <- fisher(p = Ps,adjust = "galwey",R = dependence)
      out <- res$p
    }
    return(out)
  }) 
  fisher_enrich <- data.frame("path" = unique(database$pathname),
                              "combinedP" = fisherPs)
  fisher_enrich <- fisher_enrich[order(fisher_enrich$combinedP),]
  fisher_enrich$fdr <- p.adjust(fisher_enrich$combinedP,method = "BH")
  return(fisher_enrich)
}

sampleBasedCorrelation <- function(expr, database, method = "pearson", minGenes = 10, minExplained = 0.00){
  require(reshape2)
  expression <- makePathwayExpression(pathwayInfo = database,
                                      expression = expr,
                                      minGenes = minGenes,
                                      minExplained = minExplained)
  cor_matr <- cor(t(expression),method = method)
  cor_matr[upper.tri(cor_matr)] <- NA
  diag(cor_matr) <- NA
  molten <- melt(cor_matr)
  molten <- molten[!is.na(molten$value),]
  molten$name <- make.names(paste0(molten$Var1,"_",molten$Var2))
  out <- data.frame("Individual" = molten$value)
  rownames(out) <- molten$name
  rownames(out) <- molten$name
  return(out)
}

multiCore_SBCM <- function(individualExpressions, database, method, minGenes, minExplained, cores){
  require(foreach)
  require(doParallel)
  require(parallel)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  exports <- c("sampleBasedCorrelation","makePathwayExpression")
  sampleCorsList <- foreach(ID_expr=individualExpressions,
                            .export = exports) %dopar% {
                              require(Matrix)
                              sampleBasedCorrelation(expr = ID_expr,
                                                     database = database,
                                                     method = method,
                                                     minGenes = minGenes,
                                                     minExplained = minExplained)
                            }
  stopCluster(cl)
  return(sampleCorsList)
}

singleCore_SBCM <- function(individualExpressions, database, method, minGenes, minExplained){
  sampleCorsList <- lapply(individualExpressions,function(ID_expr){
    sampleCor <- sampleBasedCorrelation(expr = ID_expr,
                                        database = database,
                                        method = method,
                                        minGenes = minGenes,
                                        minExplained = minExplained)
    
  })
  return(sampleCorsList)
}


sampleBasedCorrelationMatrix <- function(seuratObj, sampleColumn, database, method = "pearson", minGenes = 10, minExplained = 0.00, cores = 1){
  require(Seurat)
  require(SeuratObject)
  if(seuratObj@assays$RNA@data[1,][seuratObj@assays$RNA@data[1,]!=0][1]%%1 == 0){
    stop("Your data is not normalized. Please normalize your data!")
  }
  allSampleIDs <- unique(seuratObj@meta.data[,sampleColumn])
  
  
  individualExpressions <- vector("list", length = length(allSampleIDs))
  
  for(i in 1:length(allSampleIDs)){
    cellIDs <- rownames(seuratObj@meta.data[seuratObj@meta.data[,sampleColumn] == allSampleIDs[i],])
    individualExpressions[[i]] <- seur@assays$RNA@data[,cellIDs]
  }
  
  if(cores == 1){
    message("Running on 1 core")
    sampleCorsList <- singleCore_SBCM(individualExpressions,
                                      database,
                                      method,
                                      minGenes,
                                      minExplained)
  }else{
    message(sprintf("Running on %s cores",cores))
    sampleCorsList<- multiCore_SBCM(individualExpressions,
                                    database,
                                    method,
                                    minGenes,
                                    minExplained,
                                    cores)
  }
  
  calculatedPerSample <- lapply(sampleCorsList,rownames) |> unlist() |> table()
  keep <- names(calculatedPerSample)[calculatedPerSample == max(calculatedPerSample)]
  sampleCors <- lapply(sampleCorsList,function(x)x[keep,]) |> do.call(what = "cbind")
  rownames(sampleCors) <- keep
  colnames(sampleCors) <- allSampleIDs
  return(sampleCors)
}


makeNetwork <- function(res){
  require(igraph)
  links <- data.frame(source = as.character(res$GeneA),
                      target = as.character(res$GeneB))
  
  genes <- unique(c(links$source, links$target))
  nodes <- data.frame(nodes = genes)
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  return(network)
}

getCentrality <- function(network){
  centra <- as.data.frame(degree(network))
  centra$genes <- rownames(centra)
  colnames(centra) <- c("degree","genes")
  centra <- centra[order(centra$degree,decreasing = T),]
  centra$index <- 1:nrow(centra)
  return(centra)
}

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = colnames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
diffCor <- function(rho1, n1, rho2, n2,method="pearson"){
  zr1 <- atanh(rho1)
  zr2 <- atanh(rho2)
  if(method == "pearson"){
    div <- sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
    zdiff <- (zr2 - zr1)
    output <- zdiff / div
  }
  if(method == "spearman"){
    div <- sqrt((1.06/(n1 - 3)) + (1.06/(n2 - 3)))
    zdiff <- (zr2 - zr1)
    output <- zdiff / div
  }
  return(output)
}

getPvalue <- function(ZscoresDiff){
  out <- 2*pnorm(-abs(ZscoresDiff))
  return(out)
}


permute_mtc2 <- function(data, samplesheet, genesA,genesB,contrasCol, contrast,method){
  require(reshape2)
  #Permute labels
  samplesheet[,contrasCol] <- sample(samplesheet[,contrasCol])
  ref_cor <- cor(t(data[genesA,samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"]]),
                 t(data[genesB,samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"]]),
                 method = method)
  
  alt_cor <- cor(t(data[genesA,samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"]]),
                 t(data[genesB,samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"]]),
                 method = method)
  
  ref_cor <- melt(ref_cor)
  alt_cor <- melt(alt_cor)
  colnames(ref_cor) <- c("GeneA","GeneB","ct_cor")
  ref_cor$AD_cor <- alt_cor$value
  ref_cor$ct_N <- length(samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"])
  ref_cor$ad_N <- length(samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"])
  ref_cor$zdiff <- diffCor(ref_cor$ct_cor,ref_cor$ct_N,ref_cor$AD_cor,ref_cor$ad_N,method = method)
  ref_cor$P <- getPvalue(ref_cor$zdiff)
  return(ref_cor$P)
  
}

permute_mtc <- function(data, samplesheet,contrasCol, contrast,method){
  #Permute labels
  samplesheet[,contrasCol] <- sample(samplesheet[,contrasCol])
  ref_cor <- cor(t(data[,samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"]]),method = method)
  alt_cor <- cor(t(data[,samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"]]),method = method)
  ref_cor <- flattenCorrMatrix(ref_cor)
  alt_cor <- flattenCorrMatrix(alt_cor)
  colnames(ref_cor) <- c("GeneA","GeneB","ct_cor")
  ref_cor$AD_cor <- alt_cor$cor
  ref_cor$ct_N <- length(samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"])
  ref_cor$ad_N <- length(samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"])
  ref_cor$zdiff <- diffCor(ref_cor$ct_cor,ref_cor$ct_N,ref_cor$AD_cor,ref_cor$ad_N,method = method)
  ref_cor$P <- getPvalue(ref_cor$zdiff)
  return(ref_cor$P)
  
}

difcorRun <- function(data, samplesheet,contrasCol, contrast,method, padjust, Nperm=NULL){
  ref_cor <- cor(t(data[,samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"]]),method = method)
  alt_cor <- cor(t(data[,samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"]]),method = method)
  ref_cor <- flattenCorrMatrix(ref_cor)
  alt_cor <- flattenCorrMatrix(alt_cor)
  colnames(ref_cor) <- c("GeneA","GeneB","ct_cor")
  ref_cor$AD_cor <- alt_cor$cor
  ref_cor$ct_N <- length(samplesheet[samplesheet[,contrasCol] == contrast[1],"colnames"])
  ref_cor$ad_N <- length(samplesheet[samplesheet[,contrasCol] == contrast[2],"colnames"])
  ref_cor$zdiff <- diffCor(ref_cor$ct_cor,ref_cor$ct_N,ref_cor$AD_cor,ref_cor$ad_N,method = method)
  ref_cor$P <- getPvalue(ref_cor$zdiff)
  
  if(padjust == "fdr"){
    ref_cor$fdr <- p.adjust(ref_cor$P,method = "fdr")
    ref_cor <- ref_cor[order(ref_cor$fdr),]
  }else if(padjust == "perm"){
    ref_cor <- ref_cor[order(ref_cor$P),]
    
    nominal_pvalues <- ref_cor$P
    perm_corrected <- vector(mode="integer", length=length(nominal_pvalues))
    
    for(i in 1:Nperm){
      if(i %% 10 == 0){message(sprintf("%s of %s permutations done...",i,Nperm))}
      Permuted_pvalues <- permute_mtc(data,samplesheet,contrasCol,contrast,method)
      Permuted_pvalues <- sort(Permuted_pvalues)
      perm_corrected <- perm_corrected + findInterval(nominal_pvalues, Permuted_pvalues)
    }
    perm_corrected <- perm_corrected / (length(nominal_pvalues)*Nperm)
    ref_cor$perm_P <- perm_corrected
  }
  return(ref_cor)
}

variantsToTest <- function(cell){
  vars <- readRDS("./data/ROSMAP_SEAAD_AD_107.rds")[[2]]
  samplesheet <- readRDS(sprintf("./data/ROSMAP_SEAAD_MIT_%s.rds",cell))[[2]]
  vars <- vars[,colnames(vars)%in% samplesheet$IDs]
  dim(vars)
  minVar <- lapply(rownames(vars),function(x){
    dosages <- c(round(vars[x,])) |> unlist()
    min(table(dosages) / sum(table(dosages)))
  }) |> unlist()
  nVars <- lapply(rownames(vars),function(x){
    dosages <- c(round(vars[x,])) |> unlist()
    length(table(dosages))
  }) |> unlist()
  minVar <- data.frame(vars = rownames(vars),min = minVar,nVars = nVars)
  toTest <- minVar[minVar$min>=0.05 & minVar$nVars ==3,"vars"]
  vars <- as.matrix(vars)
  var_cors <- cor(t(vars))
  var_cors[upper.tri(var_cors)] <- NA
  var_cors <- melt(var_cors)
  var_cors <- var_cors[!is.na(var_cors$value),]
  var_cors <- var_cors[var_cors$Var1 != var_cors$Var2,]
  var_cors <- var_cors[order(abs(var_cors$value),decreasing = T),]
  removeVars <- as.character(unique(var_cors[abs(var_cors$value)>=0.8,"Var2"]))
  toTest <- toTest[!toTest%in%removeVars]
  return(toTest)
}

prepareData <- function(celltype){
  vars <- readRDS("./Data/ROSMAP_SEAAD_AD_107.rds")[[2]]
  expression <- readRDS(sprintf("./Data/ROSMAP_SEAAD_MIT_%s.rds",celltype))[[1]]
  samplesheet <- readRDS(sprintf("./Data/ROSMAP_SEAAD_MIT_%s.rds",celltype))[[2]]
  dim(expression)
  inBoth <- intersect(colnames(vars),colnames(expression))
  vars <- vars[variantsToTest(celltype),inBoth]
  expression <- expression[,inBoth]
  samplesheet <- samplesheet[match(inBoth,samplesheet$IDs),]
  samplesheet <- cbind(samplesheet,t(vars))
  samplesheet <- samplesheet[!is.na(samplesheet$age),]
  samplesheet[samplesheet$age == "90+","age"] <- 90
  samplesheet$age <- as.numeric(samplesheet$age)
  expression <- expression[,samplesheet$IDs]
  return(list(expression,samplesheet,vars))
}

databaseLoader <- function(type){
  if(type == "microRNA"){
    new <- rio::import("C:/Users/gabouland/OneDrive - Delft University of Technology/Documents/004 PhD/024 Proxy method/Databases/hsa_MTI.xlsx")
    table(new$`Support Type`)
    new <- new[grepl("Functional MTI",new$`Support Type`),]
    new <- new[,c("miRNA","Target Gene")]
    colnames(new) <- c("pathname","symbol")
  }else if(type == "KEGG"){
    x <- scan("./Databases/KEGG_2021_Human", what="", sep="\n")
    y <- strsplit(x, "\t")
    new <- lapply(y,function(z){
      data.frame("TF" = z[1],
                 "targets" = z[3:length(z)])
    }) |> do.call(what = "rbind")
    colnames(new) <- c("pathname","symbol")
  }else if(type == "Reactome"){
    x <- scan("./Databases/Reactome_2022.txt", what="", sep="\n")
    y <- strsplit(x, "\t")
    new <- lapply(y,function(z){
      data.frame("TF" = z[1],
                 "targets" = z[3:length(z)])
    }) |> do.call(what = "rbind")
    colnames(new) <- c("pathname","symbol")
  }else if (type == "Metabolites"){
    x <- scan("./Databases/Metabolomics_Workbench_Metabolites_2022", what="", sep="\n")
    y <- strsplit(x, "\t")
    new <- lapply(y,function(z){
      data.frame("TF" = z[1],
                 "targets" = z[3:length(z)])
    }) |> do.call(what = "rbind")
    colnames(new) <- c("pathname","symbol")
  }
  return(new)
}