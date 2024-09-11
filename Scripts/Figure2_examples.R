## Load util functions
source("./utils.R")

library(ggplot2)
library(ggbeeswarm)
library(reshape2)

## Specify gene set, variant, cell type, and database
name <- "Taurine and hypotaurine metabolism"
variant <- "rs10933431"
cell <- "Astro"
type <- "KEGG"

## Load data
new <- databaseLoader(type)
dataList <- prepareData(cell)
expression <- dataList[[1]]
samplesheet <- dataList[[2]]
vars <- dataList[[3]]
rm(dataList)
varINFO <- readRDS("./ROSMAP_SEAAD_AD_107.rds")[[1]]

## make GS variable
new <- new[new$pathname == name,]
new <- new[new$symbol %in% rownames(expression),]
evalData <- makePathwayExpression(pathwayInfo = new,
                                  expression = expression,
                                  minGenes = 5,
                                  minExplained = 0.10,
                                  linear = TRUE)

## Prepare genotype data
plotData <- data.frame(t(evalData),round(t(vars[variant,colnames(evalData)])))
colnames(plotData) <- c(name,"variant")
plotData$variant <- as.factor(plotData$variant)
counted <- varINFO[variant,"COUNTED"]
alt <- varINFO[variant,"ALT"]
genotypes <- c("0" = sprintf("%s/%s",alt,alt),
               "1" = sprintf("%s/%s",alt,counted),
               "2" = sprintf("%s/%s",counted,counted))
plotData$genotype <- unname(genotypes[plotData$variant])


## load normal QTLs to select most associated individual genes 
justQTLs <- readRDS(sprintf("U:/brainQTLs/Gerard/runProxy/normalQTLs_2.0/%s_%s_age_sex_diagnosis_dataset.rds",variant,cell))
justPaths <- justQTLs[match(new$symbol,justQTLs$name),]
justPaths <- justPaths[order(justPaths$p),]
justPaths <- unique(justPaths)
## select top 4 genes
genes <- justPaths[order(justPaths$p),"name"][1:4]
plotData <- cbind(plotData,t(expression[genes,]))
plotData <- melt(plotData,id.vars = c("variant","genotype"))

## plot
ggplot(plotData,aes(genotype,value, col = genotype)) + geom_boxplot(outlier.size = 0) +
  geom_quasirandom(size = 0.75,alpha = 0.5) + facet_wrap(~variable,scales = "free_y",nrow = 1) + theme_minimal() + labs(x = variant, y = "Expression")

