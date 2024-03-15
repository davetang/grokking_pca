
## -----------------------------------------------------------------------------
#| label: RAFisherSmoking
#| echo: false

## -----------------------------------------------------------------------------
#| label: dailies-icon
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Statisticians use the term **error** for any deviations of a measured value from the true value. This is different from the everyday use of the word. In statistics, error is an unavoidable aspect of life. It is not \"bad\", it is something to be cherished, reckoned with, tamed and controlled."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Maybe this is akin to the vision of \"personalized medicine\": better patient stratification that converts within group variation (incl. unsuccessful or unnecessary treatments) into between groups variation (where every group gets exactly what they need)."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-cointosser3-web
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: generatedata
#| cache: false
#| echo: false
library("dplyr")
bat6 = tibble(
  state  = factor(c("healthy", "disease")[rep(1:2, 6)]),
  time  = factor(rep(1:3, each = 4)),
  exprst = 0.5 * as.integer(time) + rep(c(0.5, 1), 6) + rnorm(12, 0, 0.3),
  exprs0 = rep(c(1.5, 2), 6) + rnorm(12,0,0.1),
  batch  = factor(c("Batch 1", "Batch 2")[rep(c(1, 2), 6)]))

ms0 = group_by(bat6, state) |> summarize(y = median(exprs0))

bat60 = tibble(
  state = factor(c("healthy", "disease")[rep(c(1, 2), 60)], levels = c("healthy", "disease")),
  exprs = rep(c(1.5, 2), 60) + rnorm(120, 0, 0.3))
## save(bat60, bat6, ms0, file = "../data/designI.rda")

## -----------------------------------------------------------------------------
#| label: fig-confounding-1
#| echo: false
#| column: margin
#| fig-width: 4.25
#| fig-height: 3
#| fig-cap: "Comparison of a (hypothetical) biomarker between samples from disease and healthy states. If we are only given the information shown in the left panel, we might conclude that this biomarker performs well in detecting the disease. If, in addition, we are told that the data were acquired in two separate batches (e.g., different labs, different machines, different time points) as indicated in the panel on the right hand side, the conclusion will be different."
library("ggplot2")
library("gridExtra")
library("ggbeeswarm")
## load("../data/designI.rda")
p0 = ggplot(bat6, aes(x = state, y = exprs0)) +
       geom_boxplot(alpha = 0.5, col="blue") + geom_beeswarm(size = 2, cex = 6) + # geom_point(size = 2) +
       ylab("biomarker level")
grid.arrange(p0, p0 + geom_beeswarm(aes(col = batch), size = 2, cex = 6),  # geom_point(aes(col = batch), size = 2),
  ncol = 2, widths = c(1.3, 2))

## -----------------------------------------------------------------------------
#| label: fig-Avicenna
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-effectsize
#| echo: false
#| column: margin
#| fig-width: 2
#| fig-height: 3
#| out-width: 50%
#| fig-cap: "The red arrow shows the effect size, as measured by the difference between the centers of the two groups. Here we locate the centers by the medians; sometimes the mean is used."
p0 + geom_segment(data = ms0, aes(y = y[1], yend = y[2]),
    x = 1.5, xend = 1.5, col = "red", arrow = arrow(length = unit(0.5, "cm"),
    ends = "both", type = "closed"))

## -----------------------------------------------------------------------------
#| label: fig-comparesamplesize
#| echo: false
#| column: margin
#| fig-width: 3.7
#| fig-height: 3
#| fig-cap: "On the left, the boxplot was created with samples of size 6. On the right the sample sizes are 60. The measurements have the same underlying error distribution in both cases."
p = ggplot(bat6, aes(x = state, y = exprst)) + geom_boxplot(alpha = 0.5, col = "blue") +
    ylim(c(0.5, 3)) + ylab("biomarker level")
p1 = p + geom_beeswarm(size = 2, cex = 6) # geom_point(size = 2)
p2 = p + geom_beeswarm(aes(col = time), size = 2, cex = 6) # geom_point(aes(col = time), size = 2)

mN = summarise(group_by(bat60, state), med = median(exprs))
pN = ggplot(bat60, aes(x = state, y = exprs)) + geom_boxplot(alpha = 0.5, col="blue") +
  ylab("biomarker level") + ylim(c(0.5,3)) + geom_beeswarm(size = 2, cex = 2)  +
  geom_segment(data = mN, aes(y = med[1],yend=med[2]), x = 1.5, xend = 1.5,
               col = "red", arrow = arrow(length = unit(0.5, "cm"), ends = "both", type = "closed"))
grid.arrange(p1, pN, ncol = 2, widths = c(1.6, 2.5))

## -----------------------------------------------------------------------------
#| label: fig-balance
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: HotellingsExpt
theta = round((2 * sample(8, 8) + rnorm(8)), 1)
theta

## -----------------------------------------------------------------------------
#| label: SimpleWeighing
X = theta + rnorm(length(theta), 0, 0.1)
X
errors1 = X - theta
errors1
sum(errors1^2)

## -----------------------------------------------------------------------------
#| label: HotellingsMethod
library("survey")
h8 = hadamard(6)
coef8 = 2*h8 - 1
coef8

## -----------------------------------------------------------------------------
Y = theta  %*% coef8 + rnorm(length(theta), 0, 0.1)

## -----------------------------------------------------------------------------
#| label: coef8
coef8 %*% t(coef8)
theta %*% coef8 %*% t(coef8) / ncol(coef8)

## -----------------------------------------------------------------------------
#| label: thetahat
thetahat = Y %*% t(coef8) / ncol(coef8)

## -----------------------------------------------------------------------------
#| label: Hoterrors
errors2 = as.vector(thetahat) - theta
errors2
sum(errors2^2)

## -----------------------------------------------------------------------------
#| label: bootstrapHotelling
B  = 10000
tc = t(coef8) / ncol(coef8)
sse = replicate(B, {
  theta = round((2 * sample(8, 8)) + rnorm(8), 1)
  X = theta + rnorm(length(theta), 0, 0.1)
  err1 = sum((X - theta)^2)
  Y = coef8 %*% theta + rnorm(length(theta), 0, 0.1)
  thetahat = tc %*% Y
  err2 = sum((thetahat - theta)^2)
  c(err1, err2)
})
rowMeans(sse)

## -----------------------------------------------------------------------------
#| label: fig-logsseratios
#| out-width: 50%
#| fig-width: 3.2
#| fig-height: 2
#| fig-cap: "Logarithm (base 2) of the ratios of sum of squared error for the two methods. The vertical orange line corresponds to 8."
ggplot(tibble(lr = log2(sse[1, ] / sse[2, ])), aes(x = lr)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = log2(8), col = "orange") +
  xlab("log2 ratio of SSE, Method 1 vs 2")

## -----------------------------------------------------------------------------
#| label: fig-blockbox
#| echo: false
#| column: margin
#| fig-width: 4.25
#| fig-height: 3
#| fig-cap: "On the left, two samples each of size 6 are being compared. On the right, the same data are shown, but colored by the time of data collection. We note a tendency of the data to fall into blocks according to these times. Because of this, comparison between the groups is diluted. This effect can be mitigated by comparing within times, i.,e., by blocking into three groups. Paired analysis, such as demonstrated in Questions [-@prp-design-paired]—[-@prp-design-powerPairedUnpaired], is a special case of blocking."
grid.arrange(p1, p2, ncol = 2, widths = c(1.3, 2))

## -----------------------------------------------------------------------------
#| label: fig-maizedarwin
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: zeamays
n = 15
effect = 0.2
pots   = rnorm(n, 0, 1)
noiseh = rnorm(n, 0, 0.25)
noisea = rnorm(n, 0, 0.25)
hybrid = pots + effect + noiseh
autoz  = pots + noisea

## -----------------------------------------------------------------------------
#| label: ttestpairedornot
t.test(hybrid, autoz, paired = FALSE)
t.test(hybrid, autoz, paired = TRUE)

## -----------------------------------------------------------------------------
#| label: bootstrapPower
B     = 1000
alpha = 0.05
what  = c(FALSE, TRUE)
pvs = replicate(B, {
  pots   = rnorm(n, 0, 1)
  noiseh = rnorm(n, 0, 0.25)
  noisea = rnorm(n, 0, 0.25)
  hybrid = pots + effect + noiseh
  autoz  = pots + noisea
  vapply(what,
    function(paired)
      t.test(hybrid, autoz, paired = paired)$p.value,
    double(1)) |> setNames(paste(what))
})
rowMeans(pvs <= alpha)

## -----------------------------------------------------------------------------
#| label: fig-pvaluescompare-1
#| out-width: 50%
#| fig-height: 2
#| fig-width: 3
#| fig-cap: "Results from the power calculation, comparing the p-value distributions from the ordinary unpaired and the paired $t$-test."
tidyr::pivot_longer(as.data.frame(t(pvs)), cols = everything(), names_to = "paired") |>
  ggplot(aes(x = value, fill = paired)) +
  geom_histogram(binwidth = 0.01, boundary = 0, alpha = 1/3)

## -----------------------------------------------------------------------------
#| label: powerPairedUnpaired
powercomparison = function(effect = 0.2, n = 15, alpha = 0.05,
                sdnoise, sdpots, B = 1000) {
  what = c(FALSE, TRUE)
  pvs = replicate(B, {
    pots   = rnorm(n, 0, sdpots)
    noiseh = rnorm(n, 0, sdnoise)
    noisea = rnorm(n, 0, sdnoise)
    hybrid = pots + effect + noiseh
    autoz  = pots + noisea
    vapply(what,
      function(paired)
        t.test(hybrid, autoz, paired = paired)$p.value,
      double(1)) |> setNames(paste(what))
  })
  rowMeans(pvs <= alpha)
}

## -----------------------------------------------------------------------------
powercomparison(sdpots = 0.5,  sdnoise = 0.25)
powercomparison(sdpots = 0.25, sdnoise = 0.25)
powercomparison(sdpots = 0.1,  sdnoise = 0.25)

## -----------------------------------------------------------------------------
powercomparison(sdpots = 0.5, sdnoise = 0.5, n = 100)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Beware of underpowered me-too studies."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-elephant
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: pwr.t.test
library("pwr")
str(pwr.t.test)

## -----------------------------------------------------------------------------
pwr.t.test(n = 15, d = 0.4, sig.level = 0.05, type = "two.sample")
pwr.t.test(n = 15, d = 0.4, sig.level = 0.05, type = "paired")

## -----------------------------------------------------------------------------
pwr.t.test(d = 0.4, sig.level = 0.05, type = "two.sample", power=0.8)
pwr.t.test(d = 0.4, sig.level = 0.05, type = "paired", power=0.8)

## -----------------------------------------------------------------------------
#| label: fig-effective-sample-size-sim-1
#| column: margin
#| fig-width: 4
#| fig-height: 2.5
#| fig-cap: "Density estimates for the polling result using the two sampling methods. The correlated method has higher spread. The truth is indicated by the vertical line."
doPoll = function(n = 100, numPeoplePolled = 12) {
  opinion = sort(rnorm(n))
  i1 = sample(n, numPeoplePolled)
  i2 = sample(seq(3, n, by = 3), numPeoplePolled / 3)
  i2 = c(i2, i2 - 1, i2 - 2)
  c(independent = mean(opinion[i1]), correlated = mean(opinion[i2]))
}
responses = replicate(5000, doPoll())

tidyr::pivot_longer(as.data.frame(t(responses)), 
        cols = everything(), names_to = "design") |>
ggplot(aes(x = value, col = design)) + geom_density() +
  geom_vline(xintercept = 0) + xlab("Opinion poll result")

## -----------------------------------------------------------------------------
#| label: fig-1896-Ford-Quadricycle
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Don't pretend you are dumb."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: dt2
#| cache: false
library("magrittr")
data("x", package = "Hiiragi2013")
xwdf = tibble(
  probe  = c("1420085_at", "1418863_at", "1425463_at", "1416967_at"),
  symbol = c(      "Fgf4",      "Gata4",      "Gata6",       "Sox2"))
xwdf %<>% bind_cols(as_tibble(Biobase::exprs(x)[xwdf$probe, ]))
dim(xwdf)
xwdf[, 1:5]

## -----------------------------------------------------------------------------
#| label: pivot_longer
library("tidyr")
xldf = pivot_longer(xwdf, cols = !all_of(c("probe", "symbol")),
                          names_to = "sample")
dim(xldf)
head(xldf)

## -----------------------------------------------------------------------------
#| label: fig-leakypipeline
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: vectorization1
#| cache: false
a = runif(1e6)
b = runif(length(a))
system.time({
  z1 = numeric(length(a))
  for (i in seq(along = a))
    z1[i] = a[i]^2 * b[i]
})
system.time({
  z2 = a^2 * b
})
identical(z1, z2)

## -----------------------------------------------------------------------------
#| label: vectorization2
#| echo: false
stopifnot(identical(z1, z2))

## -----------------------------------------------------------------------------
#| label: Rcpp1
#| cache: false
library("Rcpp")
cppFunction("
  NumericVector myfun(NumericVector x, NumericVector y) {
    int n = x.size();
    NumericVector out(n);
    for(int i = 0; i < n; ++i) {
      out[i] = pow(x[i], 2) * y[i];
    }
    return out;
  }")
z3 = myfun(a, b)
identical(z1, z3)

## -----------------------------------------------------------------------------
#| label: Rcpp2
#| echo: false
stopifnot(identical(z1, z3))
