library(Hmisc)
library(ggplot2)
library(reshape2)
library(ggrepel)
source("C:/Users/gabouland/OneDrive - Delft University of Technology/Documents/004 PhD/024 Proxy method/Github/Scripts/utils.R")

## Does proxy value represent real value
cell <- "Ex"

expression <- readRDS(sprintf("C:/Users/gabouland/OneDrive - Delft University of Technology/Documents/004 PhD/024 Proxy method/Github/Data/ROSMAP_SEAAD_MIT_%s.rds",cell))[[1]]
samplesheet <- readRDS(sprintf("C:/Users/gabouland/OneDrive - Delft University of Technology/Documents/004 PhD/024 Proxy method/Github/Data/ROSMAP_SEAAD_MIT_%s.rds",cell))[[2]]


samplesheet <- samplesheet[!is.na(samplesheet$IDs),]
expression <- expression[,samplesheet$IDs]

## Load TF database
x <-scan("C:/Users/gabouland/OneDrive - Delft University of Technology/Documents/004 PhD/024 Proxy method/Github/Databases/TRRUST_Transcription_Factors_2019", what="", sep="\n")
y <- strsplit(x, "\t")
new <- lapply(y,function(z){
  data.frame("TF" = z[1],
             "targets" = z[3:length(z)])
}) |> do.call(what = "rbind")
colnames(new) <- c("pathname","symbol")
## remove mouse TFs
new <-new[!grepl("mouse",new$pathname),]
new$pathname <- gsub(" human","",new$pathname)
new <- new[new$pathname %in% rownames(expression),]
new <- new[new$symbol %in% rownames(expression),]
new<- new[new$pathname!=new$symbol,]
toKeep <- unique(c(new$pathname,new$symbol))
expression <- expression[toKeep,]
rownames(expression)[1]

## regress out
expr2 <- lapply(rownames(expression),function(x){
  getResiduals(expression,samplesheet,x,c("age","sex","dataset","braak"))
}) |> do.call(what = "rbind")
rownames(expr2) <- rownames(expression)
colnames(expr2) <- colnames(expression)

expression <- expr2
## Linear ##
TF_expression <- makePathwayExpression(pathwayInfo = new,
                                       expression = expression,
                                       minGenes = 5,
                                       minExplained = 0.05,
                                       linear = TRUE)



expr_sub <- expression[rownames(TF_expression),]
rownames(TF_expression) <- paste0("proxy_",rownames(TF_expression))
corRes <- rcorr(x = t(rbind(expr_sub,TF_expression)),type = "spearman")
cors <- data.frame("TF" = rownames(corRes$r[rownames(expr_sub),rownames(TF_expression)]),
                   "abs_cor" = abs(diag(corRes$r[rownames(expr_sub),rownames(TF_expression)])),
                   "P" = abs(diag(corRes$P[rownames(expr_sub),rownames(TF_expression)])))


cors <- cors[order(cors$abs_cor,decreasing = T),]
cors$n_targets <- sapply(cors$TF,function(x){nrow(new[new$pathname == x,])})
cors$var_explained <- sapply(cors$TF,function(x){
  genes <- new[new$pathname == x,"symbol"]
  PCs <- prcomp(scale(t(expression[genes,])))
  return(PCs$sdev[1] / sum(PCs$sdev))
})
cors$max_individual_gene <- sapply(cors$TF,function(x){
  genes <- c(x,new[new$pathname == x,"symbol"])
  ordered <- sort(abs(cor(t(expression[genes,]),method = "spearman")[x,]),decreasing = T)
  max(ordered[names(ordered) != x])
})

cors$mean_cor <- sapply(cors$TF,function(x){
  genes <- c(x,new[new$pathname == x,"symbol"])
  ordered <- sort(abs(cor(t(expression[genes,]),method = "spearman")[x,]),decreasing = T)
  mean(ordered[names(ordered) != x])
})
cors$fdr <- p.adjust(cors$P,method = "BH") 
cors$type <- "linear"
corsLin <- cors
corsLin$proxBetter <- ifelse(corsLin$abs_cor >= corsLin$max_individual_gene,1,0)
corsLin[corsLin$proxBetter==1,] |> dim()

cor(t(expr[new[new$pathname == "ING1","symbol"],]), expr["ING1",],method = "spearman")
plotData <- data.frame(t(expr[new[new$pathname == "ING1","symbol"],]),"TF" = expr["ING1",],"proxy" = TF_expression["proxy_ING1",])

plotData <- melt(plotData,id.vars = "TF")
ggplot(plotData,aes(TF,value)) + geom_point() + facet_wrap(~variable,scales = "free_y",nrow = 2) + geom_smooth(method = "lm") + theme_minimal()


proxBetter <- corsLin[corsLin$proxBetter==1,]
proxBetter[order(abs(proxBetter$abs_cor - proxBetter$max_individual_gene),decreasing = T),]

TF_expression
TFto <- "RB1"

true <- expression[TFto,]
proxy <- TF_expression[sprintf("proxy_%s",TFto),]
targeted_genes <- new[new$pathname == TFto,"symbol"]


plotDataNew <- cor(true,t(expression[targeted_genes,]),method = "spearman") |> abs() |> melt()
plotDataNew <- rbind(data.frame("Var1"=1,"Var2" = "gs_TF","value" = cor(true,proxy,method = "spearman") |> abs()),plotDataNew)
plotDataNew <- plotDataNew[order(plotDataNew$value,decreasing = F),]
plotDataNew$Var2 <- factor(plotDataNew$Var2,levels = plotDataNew$Var2)
p0 <- ggplot(plotDataNew,aes(value,Var2)) + geom_point() + theme_minimal()


plotData2 <- data.frame(true,proxy,t(expression[c("ATF2","RBL1","IGF1"),]))

p1 <- ggplot(plotData2,aes(true,proxy)) + geom_point(size = 0.5) + geom_smooth(method = "lm") + theme_minimal()
p2 <- ggplot(plotData2,aes(true,ATF2)) + geom_point(size = 0.5) + geom_smooth(method = "lm")+ theme_minimal()
p3 <- ggplot(plotData2,aes(true,RBL1)) + geom_point(size = 0.5) + geom_smooth(method = "lm")+ theme_minimal()
p4 <- ggplot(plotData2,aes(true,IGF1)) + geom_point(size = 0.5) + geom_smooth(method = "lm")+ theme_minimal()


library(patchwork)
TFexample1 <- p0 + (p1 / p2 / p3 / p4) + plot_layout(widths = c(4,2))




