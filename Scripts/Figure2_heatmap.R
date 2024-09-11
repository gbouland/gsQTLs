## Load main gsQTL results
files <- list.files("../Results/gsQTLs/")
overview <- lapply(files,function(y){
  results <- readRDS(sprintf("../Results/gsQTLs/%s",y))
  sum <- lapply(results,function(x){
    out <- x[x$p_adj<=0.05,]
    if(nrow(out) == 0){
      out <- NA
      return(out)
    }else{
      out$type <- y
      return(out)
    }
  }) |> do.call(what = "rbind")
}) |> do.call(what = "rbind")
overview <- na.omit(overview)
overview$celltype <- strsplit(rownames(overview),split = "[.]") |> sapply(FUN = function(x)x[1])
overview$type <- strsplit(overview$type,split = "_") |> sapply(FUN = function(x)x[1])


library(ggplot2)
overview$celltype <- factor(overview$celltype, levels = c("Ex","Inh","Astro","Oligo","OPCs",'Micro','Endo'))
ggplot(overview,aes(variant,name,fill = Estimate)) + geom_tile(col = "black") +
  facet_grid(type~celltype,scales = "free",space = "free") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.5,0.5), space = "Lab", 
                       name="Estimate")




overview$prox <- ifelse(grepl("KEGG",overview$type),"KEGG",NA)
overview$prox <- ifelse(grepl("microRNA",overview$type),"microRNA",overview$prox)
overview$prox <- ifelse(grepl("Metabolites",overview$type),"Metabolites",overview$prox)
overview$PCA <- ifelse(grepl("Kernel",overview$type),"Kernel","Linear")
overview$type <- NULL
overview$celltype <- strsplit(rownames(overview),split = "[.]") |> sapply(FUN = function(x)x[1])
overview <- overview[order(overview$p_adj),]
rownames(overview) <- NULL
write.csv2(overview, file = "./Supplementary/Supplementary_Table01.csv",row.names = FALSE)