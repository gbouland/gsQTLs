library(ggplot2)
library(ggbeeswarm)
library(reshape2)

source("./utils.R")


#"Endo"
#"Ex"
#"Inh"#
#"Astro"
#"Micro"
#"OPCs"
#"Oligo"

## specify cell type
cells <- "Endo"

## Run QTL analysis per major cell type
lapply(cells,function(cell){
  ## Load genetic variants
  vars <- readRDS("./ROSMAP_SEAAD_AD_107.rds")[[2]]
  ## Load samplesheet 
  samplesheet <- readRDS(sprintf("./ROSMAP_SEAAD_%s.rds",cell))[[2]]
  vars <- vars[,colnames(vars)%in% samplesheet$IDs]
  
  ## Filter variants
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
  
  vars <- readRDS("./ROSMAP_SEAAD_AD_107.rds")[[2]]
  expression <- readRDS(sprintf("./ROSMAP_SEAAD_%s.rds",cell))[[1]]
  samplesheet <- readRDS(sprintf("./ROSMAP_SEAAD_%s.rds",cell))[[2]]
  
  ## Match IDs between datasets
  inBoth <- intersect(colnames(vars),colnames(expression))
  vars <- vars[toTest,inBoth]
  expression <- expression[,inBoth]
  samplesheet <- samplesheet[match(inBoth,samplesheet$IDs),]
  samplesheet <- cbind(samplesheet,t(vars))
  samplesheet <- samplesheet[!is.na(samplesheet$age),]
  samplesheet[samplesheet$age == "90+","age"] <- 90
  samplesheet$age <- as.numeric(samplesheet$age)
  expression <- expression[,samplesheet$IDs]
  
  ## Calculate PCs to correct for
  PCs <- prcomp(t(expression))
  PCs <- PCs$x
  colnames(PCs) <- paste0("mPCs",colnames(PCs))
  samplesheet <- cbind(samplesheet,PCs)
  ## Specify covariates
  toCov <- c("age","sex","diagnosis","dataset",colnames(PCs)[1:5])
  ## Run QTLs
  justQTLs <- lapply(rownames(vars), function(x){
    message(x)
    res <- runLinearModel(expression = expression,
                          samplesheet = samplesheet,
                          predictor = x,
                          cov = toCov,
                          padj = "fdr",
                          type = "linear")
    res$variant <- x
    saveRDS(res, file = sprintf("./justQTLs/%s_%s_age_sex_diagnosis_PC1_PC_5.rds",x,cell))
    return(NA)
  })
})






