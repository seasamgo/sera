
<!-- README.md is generated from README.Rmd. Please edit that file -->
README
------

This package contains wrapper functions for a Spatial Expression Regression Algorithm (SERA) to infer upon unmeasured spatial patterns of gene expression through the integration of single-cell RNA-seq (scRNA-seq) and sequential fluorescence *in situ* hybridization (seqFISH) data.

Installation
------------

``` r
install.packages("devtools", repos="http://cran.rstudio.com/")
library(devtools)
devtools::install_github("seasamgo/sera")
library(sera)
```

Methods
-------

Methods in the `sera` package expect processed cell by gene seqFISH and scRNA-seq expression matrices, having been z-score standardized by rows and columns. SERA is implemented in three steps

1.  `quantileNormalize`: ensure statistically similar distributions
2.  `predictExpression`: multi-response elastic net regression
3.  `determineEstimableGenes`: local polynomial regression of estimate variances against the L~1 norm of corresponding model coefficient estimates

``` r
quantileNormalize(
  seqFISH_expression,        # gene by cell expression matrix to be normalized
  scRNAseq_expression,       # gene by cell expression matrix for target
  normalize_scRNAseq = TRUE  # quantile-normalize scRNA-seq expression first
)

predictExpression <- function(
  training_x,                # gene by cell training predictor matrix
  training_y,                # gene by cell training outcome matrix
  prediction_x,              # gene by cell testing predictor matrix
  lambda_mse = TRUE,         # lambda minimizing the MSE of cross-validation
  family = 'mgaussian',      # multiple responses for regression
  ...                        # other arguments that may be passed to glmnet
)
  
determineEstimableGenes <- function(
  outcome,                   # prediction matrix from predictExpression result
  fit,                       # model fit from predictExpression result
  quantile_threshold = .75   # percentile to threshold estimate variances
)
```

The `predictExpression` function uses cross-validation implemented in the `glmnet` package to select lambda. Register parallel for this, e.g.

``` r
library(foreach)
library(doParallel)
registerDoParallel(10) 
```

Example
-------

To illustrate, we'll use published [seqFISH+](https://www.nature.com/articles/s41586-019-1049-y) and [scRNA-seq](https://www.nature.com/articles/nn.4216) data from mouse visual cortex (MVC).

``` r
library(Matrix)
library(methods)
data(mvc)
summary(mvc)
#>             Length Class  Mode
#> seqFISHplus 5      -none- list
#> scRNAseq    1      -none- list
```

We'll also need some additional packages to visualize what's going on along the way.

``` r
library(ggplot2)
library(pheatmap)
```

To keep this simple, subset the data to marker genes only and randomly sample some of them from each cell type for validation purposes. Then we'll add in some random genes.

``` r
set.seed(0)
markers <- as.character(mvc$seqFISHplus$markers$gene)
markers <- markers[markers%in%colnames(mvc$seqFISHplus$expression)]
type <- mvc$seqFISHplus$markers$cluster
validation_genes <- c()
for(i in unique(type))
  validation_genes <- c(validation_genes, sample(markers[type == i], size = 2))
training_genes <- markers[!markers %in% validation_genes]
validation_genes <- c(validation_genes, sample(mvc$seqFISHplus$genes[!mvc$seqFISHplus$genes%in%validation_genes], size = 50))

mvc$seqFISHplus$expression <- mvc$seqFISHplus$expression[,c(validation_genes, training_genes)]
mvc$scRNAseq$expression <- mvc$scRNAseq$expression[,c(validation_genes, training_genes)]
```

Now apply SERA to estimate our sampled genes using the remaining genes as predictors. First, we quantile-normalize the seqFISH+ genes to the scRNA-seq genes.

``` r
qn_data <- sera::quantileNormalize(
  seqFISH_expression = methods::as(mvc$seqFISHplus$expression, 'matrix'),
  scRNAseq_expression = methods::as(mvc$scRNAseq$expression, 'matrix'),
)

stats::qqplot(
  qn_data$scRNAseq_expression[,1], 
  qn_data$seqFISH_expression[,1],
  xlab = 'scRNA-seq', ylab = 'seqFISH', main = paste(markers[1], 'quantiles')
)
abline(a = 0, b = 1, col = 2)
```

![](README-unnamed-chunk-8-1.png)

Then estimate gene expression using multi-response elastic net regression. We can select the `alpha` between 0 and 1 to favor either ridge or LASSO regression. As many cell type markers are highly correlated, use a mixture that favors LASSO.

``` r
estimates <- sera::predictExpression(
  training_x = qn_data$scRNAseq_expression[, training_genes],
  training_y = qn_data$scRNAseq_expression[, validation_genes],
  prediction_x = qn_data$seqFISH_expression[, training_genes],
  alpha = .9
)
#> computing lambda 
#> training model 
#> predicting expression
```

Let's see how these estimates compare to the observed expression values.

``` r
correlation_score <- data.frame(
  Correlation = round(diag(cor(estimates$outcome, qn_data$seqFISH_expression[, validation_genes])), 2)
)

ggplot2::ggplot(data = correlation_score, aes(x = Correlation)) + geom_histogram(binwidth = .2) + ylab('Genes')
```

![](README-unnamed-chunk-10-1.png)

``` r
correlation_score[correlation_score$Correlation > .5, , drop = FALSE]
#>          Correlation
#> Hs3st4          0.51
#> Zfp385a         0.60
#> Acsbg1          0.68
#> Gja1            0.80
#> Slc27a1         0.57
#> Gsn             0.59
#> Arhgef10        0.59
#> Pld4            0.91
#> Fgfr3           1.00
#> Klhl5           1.00
#> S1pr1           0.94
```

The Pearson correlation is pretty good for several of the estimated expression vectors. We'll utilize the `pheatmap` package to visualize what's going on.

``` r
fontsize = 8
paletteLength = 100
colors <- colorRampPalette(c("darkblue", 'white', "red"))(paletteLength)
lower_bound <- -3
upper_bound <- 3
breaks <- unique(c(
  seq(lower_bound, 0, length.out=ceiling(paletteLength/2) + 1),
  seq(abs(lower_bound)/paletteLength, abs(lower_bound), length.out=floor(paletteLength/2))
))

## Observed

ph <- pheatmap::pheatmap(
  t(qn_data$seqFISH_expression[, validation_genes]), 
  show_colnames = F, show_rownames = T, 
  cluster_rows = T, cluster_cols = T,
  clustering_distance_rows = 'correlation', 
  clustering_distance_cols = 'correlation',
  scale = 'none', treeheight_row = 0, treeheight_col = 0, 
  color = colors, fontsize = fontsize, breaks = breaks
)
```

![](README-unnamed-chunk-11-1.png)

``` r

## Estimated

pheatmap::pheatmap(
  t(estimates$outcome[, validation_genes])[ph$tree_row$order, ph$tree_col$order], 
  show_colnames = F, show_rownames = T, 
  cluster_rows = F, cluster_cols = F,
  clustering_distance_rows = 'correlation', 
  clustering_distance_cols = 'correlation',
  scale = 'none', treeheight_row = 0, treeheight_col = 0, 
  color = colors, fontsize = fontsize, breaks = breaks
)
```

![](README-unnamed-chunk-11-2.png)

SERA appears to have captured some of the overall expression patterns. Consider filtering the estimates using the third step of SERA.

``` r
decision_rule <- sera::determineEstimableGenes(
  outcome = estimates$outcome,
  fit = estimates$fit,
  quantile_threshold = .75
)

decision_rule$components$correlation_score <- correlation_score$Correlation
threshold = stats::quantile(decision_rule$components$stdev, probs = .75)

ggplot2::ggplot(decision_rule$components, aes(x = l1norm, y = stdev, colour = correlation_score )) +
  geom_point() + geom_hline(yintercept = threshold, lty = 2) +
  geom_line(aes(l1norm, local_mean), colour = 'black', lty = 1, size = 1) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  ylab('Estimate Std. Dev.') + xlab('L1 Norm of Beta')
```

![](README-unnamed-chunk-12-1.png)

``` r

decision_rule$estimable_genes
#>  [1] "She"      "Hs3st4"   "Vamp1"    "Kcnab3"   "Olfml3"   "Cx3cr1"  
#>  [7] "Acsbg1"   "Gja1"     "Vcan"     "Matn4"    "Gsn"      "Arhgef10"
#> [13] "Pld4"     "Tspan2"   "Lyn"      "Fgfr3"    "Klhl5"    "S1pr1"
```

Nearly all of the well-performing genes were selected according to the decision rule, while still removing most of the poor estimates. Note that SERA was designed to estimate and filter hundreds of genes and this threshold may be adjusted as desired to remove only those with low variation. The filtered results capture much of the overall expression pattern.

``` r
## Observed

ph2 <- pheatmap::pheatmap(
  t(qn_data$seqFISH_expression[, decision_rule$estimable_genes]), 
  show_colnames = F, show_rownames = T, 
  cluster_rows = T, cluster_cols = T,
  clustering_distance_rows = 'correlation', 
  clustering_distance_cols = 'correlation',
  scale = 'none', treeheight_row = 0, treeheight_col = 0, 
  color = colors, fontsize = fontsize, breaks = breaks
)
```

![](README-unnamed-chunk-13-1.png)

``` r

## Estimated

pheatmap::pheatmap(
  t(estimates$outcome[, decision_rule$estimable_genes])[ph2$tree_row$order, ph2$tree_col$order], 
  show_colnames = F, show_rownames = T, 
  cluster_rows = F, cluster_cols = F,
  clustering_distance_rows = 'correlation', 
  clustering_distance_cols = 'correlation',
  scale = 'none', treeheight_row = 0, treeheight_col = 0, 
  color = colors, fontsize = fontsize, breaks = breaks
)
```

![](README-unnamed-chunk-13-2.png)

To perform your own analysis, consider predicting domain specific genes and analyze their estimated spatial pattern using the spatial coordinates.

``` r
summary(mvc$seqFISHplus$spatial)
#>  Field.of.View    Cell.ID            X                  Y         
#>  Min.   :100   Min.   :  1.0   Min.   :-1777.00   Min.   :  40.0  
#>  1st Qu.:100   1st Qu.:131.5   1st Qu.:   54.77   1st Qu.: 651.4  
#>  Median :100   Median :262.0   Median : 2193.52   Median :1169.5  
#>  Mean   :100   Mean   :262.0   Mean   : 2069.21   Mean   :1250.2  
#>  3rd Qu.:100   3rd Qu.:392.5   3rd Qu.: 4197.05   3rd Qu.:1659.9  
#>  Max.   :100   Max.   :523.0   Max.   : 5505.12   Max.   :3240.8

gene_expression <- data.frame(Estimated = estimates$outcome[, 'Hs3st4'], Observed = qn_data$seqFISH_expression[, 'Hs3st4'])
colours <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFE9", "#FFFFCF", "#FFEBB6", "#FFD59E", "#FFC084", "#FFAA6C", "#FF9454", "#FF8040", "#FF6A2B", "#FF5215", "#FF3B05", "#FF2100", "#FF0900", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#F70000", "#E10000", "#CF0000", "#B00000", "#900000", "#6E0000", "#4E0000")

Estimated <- data.frame(x = mvc$seqFISHplus$spatial$X, y = mvc$seqFISHplus$spatial$Y, Expression = (gene_expression$Estimated))
Estimated$Data <- rep('Estimated', length(mvc$seqFISHplus$spatial$X))
Estimated$Expression <- cut(Estimated$Expression, breaks = c(-10, seq(0, 3, (3+0)/(length(colours)-2)), 10))

rescale <- function(x, max_new = 101, min_new = 0){ max_old <- max(x); min_old = min(x); v <- (max_new-min_new)/(max_old-min_old)*(x-max_old)+max_new; return(round(v))}
gene_expression$rescaled <- gene_expression$Estimated
gene_expression$rescaled[gene_expression$rescaled > 0 ] <- rescale(gene_expression$rescaled[gene_expression$rescaled > 0 ], max_new = max(gene_expression$Observed))
Scaled <- data.frame(x = mvc$seqFISHplus$spatial$X, y = mvc$seqFISHplus$spatial$Y, Expression = gene_expression$rescaled)
Scaled$Data <- rep('Estimated (scaled)', length(mvc$seqFISHplus$spatial$X))
Scaled$Expression <- cut(Scaled$Expression, breaks = c(-10, seq(0, 3, (3+0)/(length(colours)-2)), 10))
  
Observed <- data.frame(x = mvc$seqFISHplus$spatial$X, y = mvc$seqFISHplus$spatial$Y, Expression = (gene_expression$Observed))
Observed$Data <- rep('Observed', length(mvc$seqFISHplus$spatial$X))
Observed$Expression <- cut(Observed$Expression, breaks = c(-10, seq(0, 3, (3+0)/(length(colours)-2)), 10))
  
ggplot(data = rbind(Estimated, Scaled, Observed)) + geom_point(aes(x = x, y = y, colour = Expression), size = 1.5) + scale_color_manual(values = colours) + ggtitle('Hs3st4') + facet_wrap(~Data, nrow = 1) + theme(legend.position="none")
```

![](README-unnamed-chunk-14-1.png)
