
## -----------------------------------------------------------------------------
#| label: xkcdmulttest-newspapertitle
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-active-substance-discovery-robot-screening-robot
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: FDRsetup
#| cache: false
#| echo: false
library("tibble")
library("dplyr")
library("ggplot2")
library("gganimate")
library("magick")
require("transformr")

## -----------------------------------------------------------------------------
#| label: FDRcomputestuff
#| echo: false
makedata = function(px, f1, f2, pi0, xcut, what) {
  stopifnot(length(px)==length(f1), 
            length(px)==length(f2), 
            length(pi0)==1, 
            length(xcut)==1,
            length(what)==1,
            pi0>=0, pi0<=1)
  f1 = f1 * pi0
  f2 = f2 * (1-pi0)
  i1 = which(px >= xcut)
  i2 = seq(1, i1[1], by = 1L)
  maxf1 = max(f1)
  maxf2 = max(f2)
  bind_rows(  
  tibble(
    x = px[c(i1, rev(i1))],
    y = c(f1[i1], rep(0, length(i1))),
    outcome = "True Negative"),
  tibble(
    x = px[c(i2, rev(i2))],
    y = c(f1[i2], rep(0, length(i2))),
    outcome = "False Positive"),
  tibble(
    x = px[c(i1, rev(i1))],
    y = c(f2[i1], rep(0, length(i1))) + maxf1,
    outcome = "False Negative"),
  tibble(
    x = px[c(i2, rev(i2))],
    y = c(f2[i2], rep(0, length(i2))) + maxf1,
    outcome = "True Positive"),
  tibble(
    x = rep(xcut, 3L),
    y = c(0, maxf1+maxf2, 0),
    outcome = ""
  )) |>
  bind_cols(tibble(xcut = xcut, what = what))
} 

findclosest = function(x, x0) {x[which.min(abs(x-x0))]}

pi0 = 2/3
t_df = 4
pxa = seq(-4, 4, length.out = 500)
pxb = pt(pxa, df = t_df)
xcuta = findclosest(pxa, qt(0.05, df = t_df))
xcutb = findclosest(pxb,    0.05)
f1a = dt(pxa, df = t_df) 
f2a = dgamma(pxa + 4, shape = 2, rate = 0.8) 

chainrulefac = (diff(pxa)/diff(pxb)) |> {\(x) c(x, last(x))}()
f1b = f1a * chainrulefac |> {\(x) x/sum(x)}()
f2b = f2a * chainrulefac |> {\(x) x/sum(x)}()
f1b = f1b/sum(f1b)
f2b = f2b/sum(f2b)

df = bind_rows(
  makedata(pxa, f1a, f2a, pi0, xcuta, "x"),
  makedata(pxb, f1b, f2b, pi0, xcutb, "p-value")
) 

make_static_plot = function(df) {
  stopifnot(nrow(df)>=3)
  colpal = setNames(
      c(RColorBrewer::brewer.pal(12, "Paired")[c(6,5,1,2)], "black"),
      c("True Positive", "False Negative",
        "False Positive", "True Negative", ""))
  ggplot(df, aes(x = x, y = y, fill = outcome, col = outcome)) + 
    geom_polygon() +
    scale_fill_manual(values = colpal) +
    scale_colour_manual(values = colpal) +
    xlab("value") +
    theme(legend.position = "bottom",
          legend.title = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

## -----------------------------------------------------------------------------
#| label: fig-testing-FDRvspstatic1
#| echo: false
#| column: margin
#| fig-width: 5.6
#| fig-height: 3
#| fig-cap: "Making a binary (yes/no) decision. Here, we call the two possible decisions \"positive\" and \"negative\" based on some continuous-valued score $x$, shown along the $x$-axis. The curve shaded in blue shows the distribution density of $x$ for one of the classes (the negatives), the curve shaded in red, for the other class (the positives). The distributions are distinctive (the red values are generally lower), but have some overlap. The vertical black bar marks some choice of a decision boundary, which results in four possible outcomes highlighted by the color key."
make_static_plot(dplyr::filter(df, what == "x"))

## -----------------------------------------------------------------------------
#| label: fig-testing-FDRvspstatic2
#| echo: false
#| column: margin
#| fig-width: 5.6
#| fig-height: 3
#| fig-cap: "Analogous to @fig-testing-FDRvspstatic1, but now we have transformed $x$ from its original range to the range $[0,1]$ using a non-linear, strictly increasing transformation function $p=f(x)$, which we chose such that the resulting blue distribution is uniform. Such a function always exists: it is the cumulative distribution function of $x$ (we have seen it in @sec-graphics-ecdf). We call the result a **p-value**. The definition of the FDR in @eq-testing-simplefdr applies equally well in @fig-testing-FDRvspstatic1 and here."
make_static_plot(dplyr::filter(df, what == "p-value"))

## -----------------------------------------------------------------------------
#| label: fig-testing-FDRvspanim
#| echo: false
#| column: margin
#| fig-show: asis
#| fig-width: 5.6
#| fig-height: 3
#| fig-cap: "The animation highlights the analogies between using a generic score $x$ (as in @fig-testing-FDRvspstatic1) and a p-value from a formal hypothesis test (as in @fig-testing-FDRvspstatic2) for decision making. We will come back to these concepts in terms of the two-group model in @sec-testing-localfdr and @fig-testing-lfdr."
p1 <- make_static_plot(df) + 
  labs(title = "{closest_state}") +
  transition_states(what,
                    transition_length = 3,
                    state_length = 1) + 
  view_follow() +
  ease_aes("cubic-in-out")
animate(p1, renderer = magick_renderer(), width = 5.6, height = 3, units = "in", res = 150, device = "png")

## -----------------------------------------------------------------------------
#| label: whatprob1
set.seed(0xdada)
numFlips = 100
probHead = 0.6
coinFlips = sample(c("H", "T"), size = numFlips,
  replace = TRUE, prob = c(probHead, 1 - probHead))
head(coinFlips)

## -----------------------------------------------------------------------------
#| label: tableCoinFlips
table(coinFlips)

## -----------------------------------------------------------------------------
#| label: binomDens
library("dplyr")
k = 0:numFlips
numHeads = sum(coinFlips == "H")
binomDensity = tibble(k = k,
     p = dbinom(k, size = numFlips, prob = 0.5))

## -----------------------------------------------------------------------------
#| label: fig-testing-dbinom
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: !expr sprintf("The binomial distribution for the parameters $n=%s$ and $p=0.5$, according to @eq-testing-dbinom.", numFlips)
library("ggplot2")
ggplot(binomDensity) +
  geom_bar(aes(x = k, y = p), stat = "identity") +
  geom_vline(xintercept = numHeads, col = "blue")

## -----------------------------------------------------------------------------
#| label: fig-rbinom
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: !expr sprintf("An approximation of the binomial distribution from $%s$ simulations (same parameters as @fig-testing-dbinom).", numSimulations)
numSimulations = 10000
outcome = replicate(numSimulations, {
  coinFlips = sample(c("H", "T"), size = numFlips,
                     replace = TRUE, prob = c(0.5, 0.5))
  sum(coinFlips == "H")
})
ggplot(tibble(outcome)) + xlim(-0.5, 100.5) +
  geom_histogram(aes(x = outcome), binwidth = 1, center = 50) +
  geom_vline(xintercept = numHeads, col = "blue")

## -----------------------------------------------------------------------------
#| label: fig-testing-findrej
#| column: margin
#| fig-height: 3
#| fig-width: 3.5
#| fig-cap: !expr sprintf("As @fig-testing-dbinom, with rejection region (red) that has been chosen such that it contains the maximum number of bins whose total area is at most $\\alpha=%s$.", alpha)
library("dplyr")
alpha = 0.05
binomDensity = arrange(binomDensity, p) |>
        mutate(reject = (cumsum(p) <= alpha))

ggplot(binomDensity) +
  geom_bar(aes(x = k, y = p, col = reject), stat = "identity") +
  scale_colour_manual(
    values = c(`TRUE` = "red", `FALSE` = "darkgrey")) +
  geom_vline(xintercept = numHeads, col = "blue") +
  theme(legend.position = "none")

## -----------------------------------------------------------------------------
#| label: assertion
#| echo: false
stopifnot( .tmp1 + .tmp2 != numFlips )

## -----------------------------------------------------------------------------
#| label: binom.test
binom.test(x = numHeads, n = numFlips, p = 0.5)

## -----------------------------------------------------------------------------
#| label: tbl-typesoferror
#| echo: false
#| tbl-cap-location: margin
#| tbl-cap: "Types of error in a statistical test."
dat <- data.frame(
    c('**Reject null hypothesis**', '**Do not reject**'),
    c('Type I error (false positive)', 'True negative'),
    c('True positive', 'Type II error (false negative)'))
knitr::kable(dat, col.names = c('Test vs reality', 
                                'Null hypothesis is true', 
                                '$...$ is false')) 

## -----------------------------------------------------------------------------
#| label: checkbyexperimentalmaths
#| echo: false
#| results: hide
.myttest = function(x, y) {
  mx  = mean(x)
  my  = mean(y)
  s12 = sqrt((sum((x-mx)^2)+sum((y-my)^2)) / (length(x)+length(y)-2))
  (mx - my) / s12 * sqrt(length(x)*length(y)/(length(x)+length(y)))
}
replicate(100, {
  x = rnorm(ceiling(30 * runif(1)))
  y = rnorm(ceiling(30 * runif(1)))
  stopifnot(abs(.myttest(x, y) - t.test(x, y, var.equal=TRUE)$statistic) < 1e-9)
})

## -----------------------------------------------------------------------------
#| label: fig-testing-plantgrowth
#| column: margin
#| fig-width: 3
#| fig-height: 2.75
#| fig-cap: "The `PlantGrowth` data."
library("ggbeeswarm")
data("PlantGrowth")
ggplot(PlantGrowth, aes(y = weight, x = group, col = group)) +
  geom_beeswarm() + theme(legend.position = "none")
tt = with(PlantGrowth,
          t.test(weight[group =="ctrl"],
                 weight[group =="trt2"],
                 var.equal = TRUE))
tt

## -----------------------------------------------------------------------------
#| label: fig-testing-ttestperm
#| column: margin
#| fig-width: 3
#| fig-height: 2.75
#| fig-cap: "The null distribution of the (absolute) $t$-statistic determined by simulations -- namely, by random permutations of the group labels."
abs_t_null = with(
  dplyr::filter(PlantGrowth, group %in% c("ctrl", "trt2")),
    replicate(10000,
      abs(t.test(weight ~ sample(group))$statistic)))

ggplot(tibble(`|t|` = abs_t_null), aes(x = `|t|`)) +
  geom_histogram(binwidth = 0.1, boundary = 0) +
  geom_vline(xintercept = abs(tt$statistic), col = "red")

mean(abs(tt$statistic) <= abs_t_null)

## -----------------------------------------------------------------------------
#| label: tttestpermcheck
#| echo: false
stopifnot(abs(mean(abs(tt$statistic) <= abs_t_null) -  tt$p.value) < 0.0025)

## -----------------------------------------------------------------------------
#| label: ttdup
with(rbind(PlantGrowth, PlantGrowth),
       t.test(weight[group == "ctrl"],
              weight[group == "trt2"],
              var.equal = TRUE))

## -----------------------------------------------------------------------------
#| label: prohHead_assertion
#| echo: false
stopifnot(probHead != 0.5)

## -----------------------------------------------------------------------------
#| label: tbl-mterrors
#| echo: false
#| tbl-cap-location: margin
#| tbl-cap: "Types of error in multiple testing. The letters designate the number of times each type of error occurs."
dat <- data.frame(
    c('**Rejected**', '**Not rejected**', '**Total**'),
    c('$V$', '$U$', '$m_0$'),
    c('$S$', '$T$','$m-m_0$'),
    c('$R$', '$m-R$', '$m$'))
knitr::kable(dat, col.names = c('Test vs reality', 
                                'Null hypothesis is true', 
                                '$...$ is false', 'Total'))

## -----------------------------------------------------------------------------
#| label: typeerror3
1 - (1 - 1/1e6)^8e5

## -----------------------------------------------------------------------------
#| label: fig-testing-bonferroni
#| column: margin
#| fig-width: 3
#| fig-height: 2.75
#| fig-cap: !expr sprintf("Bonferroni method. The plot shows the graph of [-@eq-testing-bonferroni] for $m=%s$ as a function of $\\alpha$.", m)
m = 10000
ggplot(tibble(
  alpha = seq(0, 7e-6, length.out = 100),
  p     = 1 - (1 - alpha)^m),
  aes(x = alpha, y = p)) +  geom_line() +
  xlab(expression(alpha)) +
  ylab("Prob( no false rejection )") +
  geom_hline(yintercept = 0.05, col = "red")

## -----------------------------------------------------------------------------
#| label: mtdeseq2airway
library("DESeq2")
library("airway")
data("airway")
aw   = DESeqDataSet(se = airway, design = ~ cell + dex)
aw   = DESeq(aw)
awde = as.data.frame(results(aw)) |> dplyr::filter(!is.na(pvalue))

## -----------------------------------------------------------------------------
#| label: fig-testing-awpvhist
#| column: margin
#| fig-width: 3
#| fig-height: 2.75
#| fig-cap: "p-value histogram of for the `airway` data."
ggplot(awde, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.025, boundary = 0)

## -----------------------------------------------------------------------------
#| label: fig-testing-awpvvisfdr
#| column: margin
#| fig-width: 3
#| fig-height: 2.75
#| fig-cap: "Visual estimation of the FDR with the p-value histogram."
alpha = binw = 0.025
pi0 = 2 * mean(awde$pvalue > 0.5)
ggplot(awde,
  aes(x = pvalue)) + geom_histogram(binwidth = binw, boundary = 0) +
  geom_hline(yintercept = pi0 * binw * nrow(awde), col = "blue") +
  geom_vline(xintercept = alpha, col = "red")

## -----------------------------------------------------------------------------
#| label: fdrvis
pi0 * alpha / mean(awde$pvalue <= alpha)

## -----------------------------------------------------------------------------
#| label: fig-testing-BH
#| column: margin
#| fig-width: 3
#| fig-height: 2.75
#| fig-cap: "Visualization of the Benjamini-Hochberg procedure. Shown is a zoom-in to the 7000 lowest p-values."
phi  = 0.10
awde = mutate(awde, rank = rank(pvalue))
m    = nrow(awde)

ggplot(dplyr::filter(awde, rank <= 7000), aes(x = rank, y = pvalue)) +
  geom_line() + geom_abline(slope = phi / m, col = "red")

## -----------------------------------------------------------------------------
#| label: kmax
kmax = with(arrange(awde, rank),
         last(which(pvalue <= phi * rank / m)))
kmax

## -----------------------------------------------------------------------------
#| label: fig-testing-SchwederSpjotvoll
#| out-width: 50%
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "Schweder and Spj\\o{}tvoll plot, as described in the answer to @prp-testing-SchwederSpjotvoll."
awdef = awde |>
  dplyr::filter(baseMean >=1) |> 
  arrange(pvalue) |>
  mutate(oneminusp = 1 - pvalue,
         N = n() - row_number())
jj = round(nrow(awdef) * c(1, 0.5))
slope = with(awdef, diff(N[jj]) / diff(oneminusp[jj]))
ggplot(awdef) +
  geom_point(aes(x = oneminusp, y = N), size = 0.15) + 
  xlab(expression(1-p[i])) +
  ylab(expression(N(p[i]))) +
  geom_abline(intercept = 0, slope = slope, col = "red3") +
  geom_hline(yintercept = slope, linetype = "dotted") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_text(x = 0, y = slope, label = paste(round(slope)), 
            hjust = -0.1, vjust = -0.25) 

## -----------------------------------------------------------------------------
#| label: fig-testing-sunexplode
#| echo: false
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: fig-testing-lfdr
#| echo: false
#| column: margin
#| fig-width: 3
#| fig-height: 6
#| fig-cap: !expr sprintf("Local false discovery rate and the two-group model, with some choice of $f_{\\text{alt}}(p)$, and $\\pi_0=%s$; densities (top) and distribution functions (bottom).", pi0)
pi0 = 0.6
f1 = function(t, shape2 = 7) {
   rv = dbeta(t, 1, shape2)
   rv / sum(rv) * (length(rv)-1) * (1-pi0)
}

t = seq(0, 1, length.out = 101)
t0 = 0.1

f0  = rep(pi0, length(t))
f   = f0 + f1(t)
F0  = cumsum(f0) / (length(t)-1)
F   = cumsum(f)  / (length(t)-1)
stopifnot(abs(F[length(F)] - 1) < 1e-2)

myplot = function(y, y0, ylim, yat, havepi0, colo = RColorBrewer::brewer.pal(12, "Paired")) {
  plot(x = t, y = y, type = "l", xlim = c(0, 1), ylim = ylim,
    xaxs = "i", yaxs = "i", ylab = "", yaxt = "n", xaxt = "n", xlab = "", main = deparse(substitute(y)))
  axis(side = 1, at = c(0, 1))
  axis(side = 2, at = yat)
  xa  =  t[t<=t0]
  xb  =  t[t>=t0]
  y0a = y0[t<=t0]
  y0b = y0[t>=t0]
  ya  =  y[t<=t0]
  polygon(x = c(xa, rev(xa)), y = c(y[t<=t0], rev(y0a)), col = colo[2])
  polygon(x = c(xb, rev(xb)), y = c(y[t>=t0], rev(y0b)), col = colo[1])
  polygon(x = c(xa, rev(xa)), y = c(rep(0, length(xa)), rev(y0a)), col = "#c0c0c0")
  polygon(x = c(xb, rev(xb)), y = c(rep(0, length(xb)), rev(y0b)), col = "#f0f0f0")
  segments(x0 = rep(t0, 2), x1 = rep(t0, 2), y0 = c(0, last(y0a)), y1 = c(last(y0a), last(ya)),
           col = colo[5:6], lwd = 3)
  text(t0, 0, adj = c(0, 1.8), labels = expression(p), cex = 1, xpd = NA)
  if (havepi0)
      text(0, pi0, adj = c(1.5, 0.5), labels = expression(pi[0]), cex = 1, xpd = NA)
}

par(mai = c(1, 0.6, 0.4, 0.3), mfcol = c(2,1))
myplot(f, f0, ylim = c(0, f[1]), yat = c(0:3),       havepi0 = TRUE)
myplot(F, F0, ylim = c(0, 1),    yat = c(0, 0.5, 1), havepi0 = FALSE)

## -----------------------------------------------------------------------------
#| label: fdrtool
#| results: hide
#| fig.show: hide
library("fdrtool")
ft = fdrtool(awde$pvalue, statistic = "pvalue")

## -----------------------------------------------------------------------------
#| label: qvalue31
ft$param[,"eta0"]

## -----------------------------------------------------------------------------
#| label: awde_basemean_counts
awde$baseMean[1]
cts = counts(aw, normalized = TRUE)[1, ]
cts
mean(cts)

## -----------------------------------------------------------------------------
#| label: makesure
#| echo: false
stopifnot(abs(mean(cts)-awde$baseMean[1])<1e-9)

## -----------------------------------------------------------------------------
#| label: fig-testing-basemean-hist
#| column: margin
#| fig-width: 3
#| fig-height: 2.4
#| fig-cap: !expr sprintf("Histogram of `baseMean`. We see that it covers a large dynamic range, from close to 0 to around %s.", round(max(awde[["baseMean"]]), -4))
ggplot(awde, aes(x = asinh(baseMean))) +
  geom_histogram(bins = 60)

## -----------------------------------------------------------------------------
#| label: fig-testing-basemean-scp
#| column: margin
#| fig-width: 3
#| fig-height: 2.4
#| fig-cap: "Scatterplot of the rank of `baseMean` versus the negative logarithm of the p-value. For small values of `baseMean`, no small p-values occur. Only for genes whose read counts across all observations have a certain size, the test for differential expression has power to come out with a small p-value."
ggplot(awde, aes(x = rank(baseMean), y = -log10(pvalue))) +
  geom_hex(bins = 60) +
  theme(legend.position = "none")

## -----------------------------------------------------------------------------
#| label: awde_stratify
awde = mutate(awde, stratum = cut(baseMean, include.lowest = TRUE,
  breaks = signif(quantile(baseMean,probs=seq(0,1,length.out=7)),2)))

## -----------------------------------------------------------------------------
#| label: fig-testing-awde-stratified-hist
#| column: margin
#| fig-width: 3
#| fig-height: 4.5
#| fig-cap: "p-value histograms of the airway data, stratified into equally sized groups defined by increasing value of `baseMean`."
ggplot(awde, aes(x = pvalue)) + facet_wrap( ~ stratum, nrow = 4) +
  geom_histogram(binwidth = 0.025, boundary = 0)

## -----------------------------------------------------------------------------
#| label: fig-testing-awde-stratified-ecdf
#| column: margin
#| fig-width: 4
#| fig-height: 3.4
#| fig-cap: "Same data as in @fig-testing-awde-stratified-hist, shown with ECDFs."
ggplot(awde, aes(x = pvalue, col = stratum)) +
  stat_ecdf(geom = "step") + theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
#| label: ihw_do
#| cache: false
library("IHW")
ihw_res = ihw(awde$pvalue, awde$baseMean, alpha = 0.1)
rejections(ihw_res)

## -----------------------------------------------------------------------------
#| label: ihwcompare
padj_BH = p.adjust(awde$pvalue, method = "BH")
sum(padj_BH < 0.1)

## -----------------------------------------------------------------------------
#| label: fig-testing-ihwplot
#| column: margin
#| fig-width: 3.5
#| fig-height: 2
#| fig-cap: !expr sprintf("Hypothesis weights determined by the `ihw` function. Here the function's default settings chose %s strata, while in our manual exploration above (Figures [-@fig-testing-awde-stratified-hist], [-@fig-testing-awde-stratified-ecdf]) we had used %s; in practice, this is a minor detail.", nrow(weights(ihw_res, levels_only=TRUE)), nlevels(awde$stratum))
plot(ihw_res)
