
## -----------------------------------------------------------------------------
#| label: xkcd-1725-linear-regression-2x
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Strictly speaking, we don't sequence the RNA but the complementary DNA (cDNA) obtained from reverse transcription. The pool of all RNA might be reduced to a subset of interest (e.,g., messenger RNA) by biochemical means, such as poly-A selection or ribosomal RNA depletion. Sensitive variants of RNA-Seq exist that enable assaying single cells, and large numbers of them."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: loadpas
#| results: hide
fn = system.file("extdata", "pasilla_gene_counts.tsv",
                  package = "pasilla", mustWork = TRUE)
counts = as.matrix(read.csv(fn, sep = "\t", row.names = "gene_id"))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "In the code shown here, we use the function `system.file` to locate a file that is shipped together with the **[pasilla](https://bioconductor.org/packages/pasilla/)** package. When you work with your own data, you will need to prepare the matrix `counts` yourself."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: counts
dim(counts)
counts[ 2000+(0:3), ]

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "There are important conceptual and practical differences between experiments and studies – see also @sec-design."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: checkClaimMadeAbove
#| echo: false
conditionNames = (sub("[[:digit:]]$", "", colnames(counts))) #$
stopifnot(length(unique(conditionNames)) == 2,
  sum(conditionNames=="untreated") == 4,
  sum(conditionNames=="treated")   == 3)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "In principle, we should consider **sampling without replacement** and the multinomial distribution here: the probability of sampling a read for the $i^{\\text{th}}$ gene depends on how many times the same gene, and other genes, have already been sampled. However, these dependencies are so negligibly small that we'll ignore them. This is because $n$ is so much larger than $r$, the number of genes is large, and each individual $n_i$ is small compared to $n$."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-countdata-normalization
#| echo: false
#| column: margin
#| results: hide
#| fig-width: 2.5
#| fig-height: 2.5
#| fig-cap: "Size factor estimation. The points correspond to hypothetical genes whose counts in two samples are indicated by their $x$- and $y$-coordinates. The lines indicate two different ways of size factor estimation explained in the text."
szfcDemo = data.frame(
  x = c(2, 4, 6, 6,  8) * 10,
  y = c(3, 6, 2, 9, 12) * 10,
  name = LETTERS[1:5],
  check.names = FALSE)
slopes =  c(
  blue = with(szfcDemo, sum(y) / sum(x)),
  red = szfcDemo[, c("x", "y")] |> as.matrix() |>
    (DESeq2::estimateSizeFactorsForMatrix)() |> (\(x) x[2]/x[1])() |> as.vector()
)
library("ggplot2")
ggplot(szfcDemo, aes(x = x, y = y, label = name)) + geom_point() +
  coord_fixed() + xlim(c(0, 128)) + ylim(c(0, 128)) + xlab("sample 1") + ylab("sample 2") +
  geom_text(hjust= 0.5, vjust = -0.6) +
  geom_abline(slope = slopes[1], col = names(slopes)[1]) +
  geom_abline(slope = slopes[2], col = names(slopes)[2])

## -----------------------------------------------------------------------------
#| label: fig-countdata-sfvssum
#| out-width: 50%
#| fig-width: 3
#| fig-height: 2.5
#| fig-cap: "Size factors versus sums for the pasilla data."
library("tibble")
library("ggplot2")
library("DESeq2")
ggplot(tibble(
  `size factor` = estimateSizeFactorsForMatrix(counts),
  `sum` = colSums(counts)), aes(x = `size factor`, y = `sum`)) +
  geom_point()

## -----------------------------------------------------------------------------
#| label: fig-countdata-varmean
#| warning.known: !expr c("rows containing non-finite values")
#| out-width: 50%
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "Variance versus mean for the (size factor adjusted) `counts` data. The axes are logarithmic. Also shown are lines through the origin with slopes 1 (green) and 2 (red)."
library("matrixStats")
sf = estimateSizeFactorsForMatrix(counts)
ncounts  = counts / matrix(sf,
   byrow = TRUE, ncol = ncol(counts), nrow = nrow(counts))
uncounts = ncounts[, grep("^untreated", colnames(ncounts)),
                     drop = FALSE]
ggplot(tibble(
        mean = rowMeans(uncounts),
        var  = rowVars( uncounts)),
     aes(x = log(mean), y = log(var))) +
  geom_hex() + coord_fixed() + theme(legend.position = "none") +
  geom_abline(slope = 1:2, color = c("forestgreen", "red"))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "In the code shown here, we load the file `pasilla_sample_annotation.csv` that comes with the **[pasilla](https://bioconductor.org/packages/pasilla/)** package. We locate it with the function `system.file`. When you work with your own data, you will need to prepare an analogous file, or directly a dataframe like `pasillaSampleAnno`."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: annotationFile
annotationFile = system.file("extdata",
  "pasilla_sample_annotation.csv",
  package = "pasilla", mustWork = TRUE)
pasillaSampleAnno = readr::read_csv(annotationFile)
pasillaSampleAnno

## -----------------------------------------------------------------------------
#| label: factors
library("dplyr")
pasillaSampleAnno = mutate(pasillaSampleAnno,
condition = factor(condition, levels = c("untreated", "treated")),
type = factor(sub("-.*", "", type), levels = c("single", "paired")))

## -----------------------------------------------------------------------------
#| label: checkfacs
#| echo: false
stopifnot(
  !any(is.na(pasillaSampleAnno$condition)),
  !any(is.na(pasillaSampleAnno$type)),
  sum(pasillaSampleAnno$type == "single") == 3,
  sum(pasillaSampleAnno$type == "paired") == 4)

## -----------------------------------------------------------------------------
#| label: condvstype
with(pasillaSampleAnno,
       table(condition, type))

## -----------------------------------------------------------------------------
#| label: DESeq2
mt = match(colnames(counts), sub("fb$", "", pasillaSampleAnno$file))
stopifnot(!any(is.na(mt)))

pasilla = DESeqDataSetFromMatrix(
  countData = counts,
  colData   = pasillaSampleAnno[mt, ],
  design    = ~ condition)
class(pasilla)
is(pasilla, "SummarizedExperiment")

## -----------------------------------------------------------------------------
#| label: ispasillaSummarizedExperiment
#| echo: false
stopifnot(is(pasilla, "SummarizedExperiment"))

## -----------------------------------------------------------------------------
#| label: deseq
pasilla = DESeq(pasilla)

## -----------------------------------------------------------------------------
#| label: theresults
res = results(pasilla)
res[order(res$padj), ] |> head()

## -----------------------------------------------------------------------------
#| label: fig-countdata-hist1
#| warning.known: !expr c("rows containing non-finite values")
#| column: margin
#| fig-width: 4.5
#| fig-height: 4.5
#| fig-cap: "Histogram of p-values of a differential expression analysis."
ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

## -----------------------------------------------------------------------------
#| label: hist2
#| echo: false
thehist = hist(res$pvalue, breaks = 100, plot=FALSE)
thehist$bgl = median(thehist$counts)

## -----------------------------------------------------------------------------
#| label: fig-countdata-MA
#| column: margin
#| fig-width: 4.5
#| fig-height: 4.5
#| fig-cap: "[MA plot](https://en.wikipedia.org/wiki/MA_plot): fold change versus mean of size-factor normalized counts. Logarithmic scaling is used for both axes. By default, points are colored red if the adjusted p-value is less than 0.1. Points which fall out of the $y$-axis range are plotted as triangles."
plotMA(pasilla, ylim = c( -2, 2))

## -----------------------------------------------------------------------------
#| label: fig-countdata-PCA
#| column: margin
#| fig-width: 4
#| fig-height: 2.5
#| fig-cap: !expr paste("PCA plot. The", ncol(pasilla), "samples are shown in the 2D plane spanned by their first two principal components.")
pas_rlog = rlogTransformation(pasilla)
plotPCA(pas_rlog, intgroup=c("condition", "type")) + coord_fixed()

## -----------------------------------------------------------------------------
#| label: fig-figHeatmap-1
#| column: margin
#| fig-width: 4
#| fig-height: 6
#| fig-cap: !expr paste("Heatmap of regularized log transformed data of the top", length(select), "genes.")
library("pheatmap")
select = order(rowMeans(assay(pas_rlog)), decreasing = TRUE)[1:30]
pheatmap( assay(pas_rlog)[select, ],
     scale = "row",
     annotation_col = as.data.frame(
        colData(pas_rlog)[, c("condition", "type")] ))

## -----------------------------------------------------------------------------
#| label: writecsv
#| eval: false
## write.csv(as.data.frame(res), file = "treated_vs_untreated.csv")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "For the normalization, although not for the dispersion estimation, one can slightly relax this assumption: it is still valid if many genes are changing, but in a way that is balanced between up- and downward directions."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Sometimes @eq-countdata-basiclm is written with an additional term $x_0$ that is multiplied with $\\beta_0$, where it is understood that $x_0=1$ always. It turns out that this makes subsequent notation and bookkeeping easier since then the intercept can be handled consistently together with the other $\\beta$s, instead of being a separate case."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Note that the addition is on the logarithmic scale, which corresponds to multiplication on the original scale."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Remember that since $\\beta_0$ is the intercept, $x_{j0}=1$ for all $j$."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "The distinction between noise and systematic variability is in the eye of the beholder, and depends on our model, not on reality."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-countdata-mestimator
#| out-width: 50%
#| fig-width: 3
#| fig-height: 2.5
#| fig-cap: "Graph of $\\rho_s(\\varepsilon)$, for a choice of $s=2$."
rho = function(x, s)
  ifelse(abs(x) < s, x^2 / 2,  s * abs(x) - s^2 / 2)

df = tibble(
  x        = seq(-7, 7, length.out = 100),
  parabola = x ^ 2 / 2,
  Huber    = rho(x, s = 2))

ggplot(reshape2::melt(df, id.vars = "x"),
  aes(x = x, y = value, col = variable)) + geom_line()

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "It is good to remember that, while we can use the normal distribution as a convenient argument to motivate least sum of squares regression through the maximum likelihood principle, the data do not have to be distributed according to the normal for least sum of squares regression to provide a useful result. In fact, least sum of squares fitting often provides useful estimates for the $\\beta$s even when the data are non-normal, although that depends on the specific circumstances."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: replaceDesign
pasillaTwoFactor = pasilla
design(pasillaTwoFactor) = formula(~ type + condition)
pasillaTwoFactor = DESeq(pasillaTwoFactor)

## -----------------------------------------------------------------------------
#| label: multiResults
res2 = results(pasillaTwoFactor)
head(res2, n = 3)

## -----------------------------------------------------------------------------
#| label: multiTypeResults
resType = results(pasillaTwoFactor,
  contrast = c("type", "single", "paired"))
head(resType, n = 3)

## -----------------------------------------------------------------------------
#| label: fig-countdata-scpres1res2
#| column: margin
#| fig-width: 4
#| fig-height: 3.5
#| fig-cap: "Comparison of p-values from the models with a single factor (condition) and with two factors (type + condition). The axes correspond to $(-\\log_{10}p)^{\\frac{1}{6}}$, an arbitrarily chosen monotonically decreasing transformation that compresses the dynamic range of the p-values for the purpose of visualization. We can see a trend for the joint distribution to lie above the bisector, indicating that the small p-values in the two-factor analysis are generally smaller than those in the one-factor analysis."
trsf = function(x) ifelse(is.na(x), 0, (-log10(x)) ^ (1/6))
ggplot(tibble(pOne = res$pvalue,
              pTwo = res2$pvalue),
    aes(x = trsf(pOne), y = trsf(pTwo))) +
    geom_hex(bins = 75) + coord_fixed() +
    xlab("Single factor analysis (condition)") +
    ylab("Two factor analysis (type + condition)") +
    geom_abline(col = "orange")

## -----------------------------------------------------------------------------
#| label: compareRes
compareRes = table(
   `simple analysis` = res$padj < 0.1,
   `two factor` = res2$padj < 0.1 )
addmargins( compareRes )

## -----------------------------------------------------------------------------
#| label: checkClaim
#| echo: false
#| error: true
#| results: hide
stopifnot(compareRes[1, 2] > compareRes[2, 1])

## -----------------------------------------------------------------------------
#| label: shrink1
#| echo: false
# The following code is based on guidance from Mike Love
res1  = results(pasilla, cooksCutoff = FALSE)
res2  = lfcShrink(pasilla, coef = "condition_treated_vs_untreated", type="normal", res = res1)

# Something like the two lines commented out below may be used to reproduce   
# how Mike Love selected the two genes: they should have similar intercepts,
# large unshrunken fold change and very different Wald statistic (i.e., have
# small / large dispersion, respectively):
#
# with(res1,
#  plot(baseMean, log2FoldChange, log = "x", ylim = c(0, 3), xlim = c(10, 1e5),
#       col = ifelse(padj < 0.1, "red", "black"), cex = log(abs(stat))))
# rownames(res1)[with(res1, identify(baseMean, log2FoldChange))]

genes = c(A = "FBgn0053926", B = "FBgn0260933")
cols  = c(FBgn0053926 = "forestgreen", FBgn0260933 = "dodgerblue3", prior = "black")

df1 = tibble(
  k         = as.vector(counts(pasilla, normalized = TRUE)[genes, ]),
  condition = rep(colData(pasilla)[["condition"]], each = length(genes)),
  gene      = rep(genes, times = ncol(pasilla)))

beta    = seq(from = -1, to = 1.5, length = 500)
kounts  = counts(pasilla)[genes,, drop = FALSE]
disps   = dispersions(pasilla[genes, ]) |> `names<-`(genes)

sf      = sizeFactors(pasilla)
cond    = as.numeric(pasilla$condition)-1 

betaPriorVar = priorInfo(res2)$betaPriorVar
priorSigma = sqrt(betaPriorVar["conditiontreated"])
prior = dnorm(beta, mean = 0, sd = priorSigma)

likelihood = function(k, alpha, intercept) {
  z = vapply(beta, function(b) {
    prod(dnbinom(k, mu = sf * 2^(intercept + b * cond), size = 1/alpha))
  }, numeric(1))
  z / (sum(z) * diff(beta[1:2]))
}

posterior = function(k, alpha, intercept) {
  z = likelihood(k, alpha, intercept) * prior
  z / (sum(z) * diff(beta[1:2]))
}

intercepts = with(mcols(pasilla[genes,]), Intercept) |> `names<-`(genes)

df2 = bind_rows(
  tibble(beta = beta, y = prior, gene = "prior", what = "pre"),
  bind_rows(
  lapply(genes, function(i) bind_rows(
    tibble(beta = beta, gene = i, what = "pre",
           y = likelihood(k = kounts[i, ], alpha = disps[i],
                          intercept = intercepts[i])),
    tibble(beta = beta, gene = i, what = "post",
           y = posterior(k = kounts[i, ], alpha = disps[i],
                         intercept = intercepts[i]))))
  )
)

is_max = function(y)
  ifelse(seq(along = y) == which.max(y), y, NA_real_)

df2 %<>% group_by(gene, what) %>% mutate(py = is_max(y))

## some consistency checks:
deseqNoPrior = res1[genes, "log2FoldChange"]
deseqPrior   = res2[genes, "log2FoldChange"]
mleFromPlot  = c(beta[which.max(likelihood(kounts[1,], disps[1], intercepts[1]))],
                 beta[which.max(likelihood(kounts[2,], disps[2], intercepts[2]))])
mapFromPlot  = c(beta[which.max( posterior(kounts[1,], disps[1], intercepts[1]))],
                 beta[which.max( posterior(kounts[2,], disps[2], intercepts[2]))])
stopifnot(all(abs(deseqNoPrior - mleFromPlot) < 0.002))

## -----------------------------------------------------------------------------
#| label: fig-countdata-posterior
#| warning.known: !expr c("rows containing missing values")
#| echo: false
#| column: margin
#| layout-ncol: 1
#| results: hide
#| fig-width: 3
#| fig-height: 2
#| fig-cap: "Shrinkage estimation of logarithmic fold change estimates by use of an empirical prior in **[DESeq2](https://bioconductor.org/packages/DESeq2/)**. Two genes with similar mean count and MLE logarithmic fold change are highlighted in green and blue. The normalized counts for these genes (a) reveal low dispersion for the gene in blue and high dispersion for the gene in green. In (b), the density plots are shown of the normalized likelihoods (solid lines) and of the posteriors (dashed lines) for the green and blue gene. In addition, the solid black line shows the prior estimated from the MLEs of all genes. Due to the higher dispersion of the green gene, its likelihood is wider and less sharp (indicating less information), and the prior has more influence on its posterior than in the case of the blue gene."
#| fig-subcap:
#|   - ""
#|   - ""
library("ggbeeswarm")
ggplot(df1, aes(x = condition, y = k, col = gene)) + geom_beeswarm(cex = 5) +
      facet_grid(. ~ gene) + ylab("normalized counts") + scale_y_log10() +
      scale_color_manual(values = cols) + theme(legend.position = "none")
ggplot(df2, aes(x = beta, col = gene, linetype = what)) +
  geom_line(aes(y = y)) + geom_point(aes(y = py)) +
  scale_color_manual(values = cols) + theme(legend.position = "none") +
  scale_linetype_manual(values = c(pre = "solid", post = "dotted")) +
  xlab(expression(beta)) + ylab("density")

## -----------------------------------------------------------------------------
#| label: defvsp
vsp = varianceStabilizingTransformation(pasilla)

## -----------------------------------------------------------------------------
#| label: fig-countdata-plotvst
#| warning.known: !expr c("containing missing values")
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "Graph of variance-stabilizing transformation for the data of one of the samples, and for comparison also of the $\\log_2$ transformation. The variance-stabilizing transformation has finite values and finite slope even for counts close to zero, whereas the slope of $\\log_2$ becomes very steep for small counts and is undefined for counts of zero. For large counts, the two transformation are essentially the same."
j = 1
ggplot(
  tibble(
    counts = rep(assay(pasilla)[, j], 2),
    transformed = c(
      assay(vsp)[, j],
      log2(assay(pasilla)[, j])
      ),
    transformation = rep(c("VST", "log2"), each = nrow(pasilla))
  ),
  aes(x = counts, y = transformed, col = transformation)) +
  geom_line() + xlim(c(0, 600)) + ylim(c(0, 9))

## -----------------------------------------------------------------------------
#| label: fig-countdata-meansd
#| warning.known: !expr c("rows containing missing values", "rows containing non-finite values")
#| fig-width: 10
#| fig-height: 4
#| fig-margin: false
#| fig-cap: "Per-gene standard deviation (sd, taken across samples) against the rank of the mean, for the shifted logarithm $\\log_2(n+1)$, the variance-stabilizing transformation (vst) and the rlog. Note that for the leftmost $\\approx$ 2,500 genes, the counts are all zero, and hence their standard deviation is zero. The mean-sd dependence becomes more interesting for genes with non-zero counts. Note also the high value of the standard deviation for genes that are weakly detected (but not with all zero counts) when the shifted logarithm is used, and compare to the relatively flat shape of the mean-sd relationship for the variance-stabilizing transformation. "
library("vsn")
rlp = rlogTransformation(pasilla)

msd = function(x)
  meanSdPlot(x, plot = FALSE)$gg + ylim(c(0, 1)) +
     theme(legend.position = "none")

gridExtra::grid.arrange(
  msd(log2(counts(pasilla, normalized = TRUE) + 1)) +
    ylab("sd(log2)"),
  msd(assay(vsp)) + ylab("sd(vst)"),
  msd(assay(rlp)) + ylab("sd(rlog)"),
  ncol = 3
)

## -----------------------------------------------------------------------------
#| label: fig-countdata-lfcThresh
#| column: margin
#| fig-width: 2.7
#| fig-height: 9
#| fig-cap: "MA-plots of tests of $\\log_2$ fold change with respect to a threshold value. From top to bottom, the tests are for `altHypothesis = \"greaterAbs\"`, `\"lessAbs\"`, `\"greater\"`, and `\"less\"`."
par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))
myMA = function(h, v, theta = 0.5) {
  plotMA(pasilla, lfcThreshold = theta, altHypothesis = h,
         ylim = c(-2.5, 2.5))
  abline(h = v * theta, col = "dodgerblue", lwd = 2)
}
myMA("greaterAbs", c(-1, 1))
myMA("lessAbs",    c(-1, 1))
myMA("greater",          1)
myMA("less",         -1   )

## -----------------------------------------------------------------------------
#| label: fig-countdata-exbatch
#| column: margin
#| fig-width: 3
#| fig-height: 4
#| fig-cap: "p-values for the tests performed on `x1` and `x2` (see code)."
library("magrittr")
ng = 10000
ns = 12
x1 = x2 = matrix(rnorm(ns * ng), ncol = ns, nrow= ng)
group = factor(letters[1 + seq_len(ns) %% 2])  %T>% print
batch = factor(ifelse(seq_len(ns) <= ns/2, "B1", "B2")) %T>% print
table(group, batch)
x2[, batch=="B2"] = x2[, batch=="B2"] + 2 * rnorm(ng)
pvals = rbind(
  cbind(type = "x1", genefilter::rowttests(x1, fac = group)),
  cbind(type = "x2", genefilter::rowttests(x2, fac = group)))
ggplot(pvals, aes(x = p.value)) + 
  geom_histogram(binwidth = 0.02, boundary = 0) +
  facet_grid(type ~ .)

## -----------------------------------------------------------------------------
#| label: ui.R
#| eval: false
## library("shiny")
## shinyUI(fluidPage(
##   titlePanel("Breakdown"),
##   sidebarLayout(
##     sidebarPanel(     # select oulier shift
##       sliderInput("shift", "Outlier:", min = 0, max = 100, value = 0),
##       radioButtons("method", "Method:",
##                    c("Non-robust least squares" = "lm",
##                      "M-estimation" = "rlm"))
##     ),
##     mainPanel(       # show fit
##       plotOutput("regPlot")
##     )
##   )
## ))

## -----------------------------------------------------------------------------
#| label: server.R
#| eval: false
## library("shiny")
## library("ggplot2")
## library("MASS")
## shinyServer(function(input, output) {
##   output$regPlot = renderPlot({
##     whpt = 15
##     mtcars_new = mtcars
##     mtcars_new$mpg[whpt] = mtcars_new$mpg[whpt] + input$shift
##     reg = switch(input$method,
##       lm = lm(mpg ~ disp, data = mtcars_new),
##       rlm = rlm(mpg ~ disp, data = mtcars_new),
##       stop("Unimplemented method:", input$method)
##     )
##     ggplot(mtcars_new, aes(x = disp, y = mpg)) + geom_point() +
##       geom_abline(intercept = reg$coefficients["(Intercept)"],
##                   slope = reg$coefficients["disp"], col = "blue")
##   })
## })
