library(faux)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(patchwork)

## Simulate data
dat <- rnorm_multi(n = 100, 
                   mu = c(1, 1, 1),
                   sd = c(1, 1, 1),
                   r = c( 1,   0.3, 0.3, 
                          0.3, 1,    0.01, 
                          0.3, 0.01,  1), 
                   varnames = c("A", "B", "C"),
                   empirical = T)

## Calculate PCs
PC <- princomp(dat[,c("C","B")])
loadings_PC1 <- PC$loadings[, 1]

## Calculate slopes and intercept of PC axis
pc1_slope <- loadings_PC1[2] / loadings_PC1[1] *-1
pc1_intercept <- mean(dat$C) - pc1_slope * mean(dat$B)


##Prepare data for plotting
plotData <- data.frame("A" = dat$A,
                       "B" = dat$B,
                       "C" = dat$C,
                       "PC1_B_C" = PC$scores[,1])

## Make continous variable A into three genotypes
plotData$A_category <- category_vector <- cut(plotData$A,
                                              breaks = c(-Inf, quantile(plotData$A, c(1/3)) , quantile(plotData$A, c(2/3)) , Inf),
                                              labels = c(0, 1, 2),
                                              right = TRUE)
plotData$A_category <- as.factor(plotData$A_category)

## Plots
main <- ggplot(plotData,aes(B,C)) +  geom_point(aes(col = A_category)) + theme_minimal() +
  geom_abline(slope = pc1_slope, intercept = pc1_intercept, color = 'red', linetype = "dashed") + theme(legend.position = "none") 
bottom <- ggplot(plotData,aes(B,A_category, col = A_category)) + geom_boxplot(outlier.shape = NA) + geom_quasirandom() + theme_minimal() + theme(legend.position = "none")
left <- ggplot(plotData,aes(A_category,C,col = A_category)) + geom_boxplot(outlier.shape = NA) + geom_quasirandom() + theme_minimal() + theme(legend.position = "none")
pcPl <- ggplot(plotData,aes(A_category,PC1_B_C,col = A_category)) + geom_boxplot(outlier.shape = NA) + geom_quasirandom() + theme_minimal() + theme(legend.position = "none")
left + main + plot_spacer() + bottom + plot_layout(ncol = 2,widths = c(1,5),heights = c(5,1))





