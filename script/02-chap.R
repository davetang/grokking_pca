
## -----------------------------------------------------------------------------
#| label: StatDiagram
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "In a statistical setting, we start with the data $X$ and use them to *estimate* the parameters. These estimates are denoted by Greek letters with what we call hats on them, as in $\\widehat{\\theta}$."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: Parameters
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-probabilitydiag
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: e99
load("../data/e100.RData")
e99 = e100[-which.max(e100)]

## -----------------------------------------------------------------------------
#| label: fig-twopoisson
#| column: margin
#| fig-height: 4
#| fig-cap: "The observed distribution of the epitope data without the outlier."
barplot(table(e99), space = 0.8, col = "chartreuse4")

## -----------------------------------------------------------------------------
#| label: fig-stat-rooto
#| column: margin
#| fig-height: 4
#| fig-cap: "Rootogram showing the square root of the theoretical values as red dots and the square root of the observed frequencies as drop down rectangles. (We'll see a bit below how the `goodfit` function decided which $\\lambda$ to use.)"
library("vcd")
gf1 = goodfit( e99, "poisson")
rootogram(gf1, xlab = "", rect_gp = gpar(fill = "chartreuse4"))

## -----------------------------------------------------------------------------
#| label: RootogramPoisson
#| fig-show: hide
simp = rpois(100, lambda = 0.05)
gf2 = goodfit(simp, "poisson")
rootogram(gf2, xlab = "")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "The parameter is called the Poisson mean because it is the mean of the theoretical distribution *and*, as it turns out, is estimated by the sample mean. This overloading of the word is confusing to everyone."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: table100
table(e100)

## -----------------------------------------------------------------------------
#| label: table3
table(rpois(100, 3))

## -----------------------------------------------------------------------------
#| label: table101
#| echo: false
counts  =  table(e100)
stopifnot(identical(names(counts), c("0", "1", "2", "7")), all(counts==c(58, 34, 7, 1)))

## -----------------------------------------------------------------------------
#| label: poism3
prod(dpois(c(0, 1, 2, 7), lambda = 3) ^ (c(58, 34, 7, 1)))

## -----------------------------------------------------------------------------
#| label: anspois
prod(dpois(c(0, 1, 2, 7), lambda = 0.4) ^ (c(58, 34, 7, 1)))

## -----------------------------------------------------------------------------
#| label: functionll
loglikelihood  =  function(lambda, data = e100) {
  sum(log(dpois(data, lambda)))
}

## -----------------------------------------------------------------------------
#| label: fig-poislikel-1
#| column: margin
#| fig-width: 3.5
#| fig-cap: "The red curve is the log-likelihood function. The vertical line shows the value of `m` (the mean) and the horizontal line the log-likelihood of `m`. It looks like `m` maximizes the likelihood."
lambdas = seq(0.05, 0.95, length = 100)
loglik = vapply(lambdas, loglikelihood, numeric(1))
plot(lambdas, loglik, type = "l", col = "red", ylab = "", lwd = 2,
     xlab = expression(lambda))
m0 = mean(e100)
abline(v = m0, col = "blue", lwd = 2)
abline(h = loglikelihood(m0), col = "purple", lwd = 2)
m0

## -----------------------------------------------------------------------------
#| label: gfpoisson
gf  =  goodfit(e100, "poisson")
names(gf)
gf$par

## -----------------------------------------------------------------------------
#| label: colorblind
#| echo: false
cb  =  c(rep(0, 110), rep(1, 10))

## -----------------------------------------------------------------------------
#| label: cb
table(cb)

## -----------------------------------------------------------------------------
#| label: meancb
mean(cb)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "However, be careful: sometimes, maximum likelihood estimates are harder to guess and to compute, as well as being much less intuitive (see @exr-models-mlmax)."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-likely1-1
#| column: margin
#| fig-width: 4
#| fig-cap: !expr sprintf("Plot of the likelihood as a function of the probabilities. The likelihood is a function on $[0, 1]$. Here we have zoomed into the range of $[%s, %s]$, as the likelihood is practically zero for larger values of $p$.", probs[1], probs[length(probs)])
probs  =  seq(0, 0.3, by = 0.005)
likelihood = dbinom(sum(cb), prob = probs, size = length(cb))
plot(probs, likelihood, pch = 16, xlab = "probability of success",
       ylab = "likelihood", cex=0.6)
probs[which.max(likelihood)]

## -----------------------------------------------------------------------------
#| label: check
#| echo: false
stopifnot(abs(probs[which.max(likelihood)]-1/12) < diff(probs[1:2]))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "One can come up with different criteria than maximum likelihood, which lead to other estimators. They all carry hats. We'll see other examples in @sec-mixtures."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: loglike1
loglikelihood = function(p, n = 300, y = 40) {
  log(choose(n, y)) + y * log(p) + (n - y) * log(1 - p)
}

## -----------------------------------------------------------------------------
#| label: fig-loglikelihood-1
#| column: margin
#| fig-width: 3
#| fig-height: 3.2
#| fig-cap: "Plot of the log likelihood function for $n=300$ and $y=40$."
p_seq = seq(0, 1, by = 0.001)
plot(p_seq, loglikelihood(p_seq), xlab = "p", ylab = "log f(p|y)", type = "l")

## -----------------------------------------------------------------------------
#| label: staph
#| cache: false
library("Biostrings")
staph = readDNAStringSet("../data/staphsequence.ffn.txt", "fasta")

## -----------------------------------------------------------------------------
#| label: firstgenestaph
staph[1]
letterFrequency(staph[[1]], letters = "ACGT", OR = 0)

## -----------------------------------------------------------------------------
#| label: compareprop
letterFrq = vapply(staph, letterFrequency, FUN.VALUE = numeric(4),
         letters = "ACGT", OR = 0)
colnames(letterFrq) = paste0("gene", seq(along = staph))
tab10 = letterFrq[, 1:10]
computeProportions = function(x) { x/sum(x) }
prop10 = apply(tab10, 2, computeProportions)
round(prop10, digits = 2)
p0 = rowMeans(prop10)
p0

## -----------------------------------------------------------------------------
#| label: outerex
cs = colSums(tab10)
cs
expectedtab10 = outer(p0, cs, FUN = "*")
round(expectedtab10)

## -----------------------------------------------------------------------------
#| label: genrandomtabs
randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) } )
all(colSums(randomtab10) == cs)

## -----------------------------------------------------------------------------
#| label: assertgenrandomtabs
#| echo: false
stopifnot(all(colSums(randomtab10) == cs))

## -----------------------------------------------------------------------------
#| label: fig-quant12-1
#| column: margin
#| fig-width: 4
#| fig-height: 3.5
#| fig-cap: "Histogram of `simulstat`. The value of `S1` is marked by the vertical red line, those of the 0.95 and 0.99 quantiles (see next section) by the dotted lines."
stat = function(obsvd, exptd) {
   sum((obsvd - exptd)^2 / exptd)
}
B = 1000
simulstat = replicate(B, {
  randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) })
  stat(randomtab10, expectedtab10)
})
S1 = stat(tab10, expectedtab10)
sum(simulstat >= S1)

hist(simulstat, col = "lavender", breaks = seq(0, 75, length.out=50))
abline(v = S1, col = "red")
abline(v = quantile(simulstat, probs = c(0.95, 0.99)),
       col = c("darkgreen", "blue"), lty = 2)

## -----------------------------------------------------------------------------
#| label: checksimulstat
#| echo: false
stopifnot(max(simulstat)<75, S1<75)

## -----------------------------------------------------------------------------
#| label: quantiles3
#| results: hide
qs = ppoints(100)
quantile(simulstat, qs)
quantile(qchisq(qs, df = 30), qs)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "A name collision occurs here. Statisticians call the summary statistic we just computed as `simulstat` (sum of squares of weighted differences), the **chi-squared** or $\\chi^2$ *statistic*. The theoretical *distribution* $\\chi^2_\\nu$ is a distribution in its own right, with a parameter $\\nu$ called the degrees of freedom. When reading about the chi-squared or $\\chi^2$, you will need to pay attention to the context to see which meaning is appropriate."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-qqplot3-1
#| column: margin
#| fig-width: 3.4
#| fig-height: 4
#| fig-cap: "Our simulated statistic's distribution compared to $\\chi_{30}^2$ using a quantile-quantile (QQ) plot, which shows the theoretical **quantiles** for the $\\chi^2_{30}$ distribution on the horizontal axis and the sampled ones on the vertical axis."
qqplot(qchisq(ppoints(B), df = 30), simulstat, main = "",
  xlab = expression(chi[nu==30]^2), asp = 1, cex = 0.5, pch = 16)
abline(a = 0, b = 1, col = "red")

## -----------------------------------------------------------------------------
#| label: pvalueBias
1 - pchisq(S1, df = 30)

## -----------------------------------------------------------------------------
#| label: ChargaffColdSpring-web
#| echo: false

## -----------------------------------------------------------------------------
#| label: Chargaff
#| cache: false
load("../data/ChargaffTable.RData")
ChargaffTable

## -----------------------------------------------------------------------------
#| label: fig-ChargaffBars
#| echo: false
#| fig-width: 12
#| fig-height: 6
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Barplots for the different rows in `ChargaffTable`. Can you spot the pattern?"
stopifnot(nrow(ChargaffTable) == 8)
mycolors = c("chocolate", "aquamarine4", "cadetblue4", "coral3",
            "chartreuse4","darkgoldenrod4","darkcyan","brown4")
par(mfrow=c(2, 4), mai = c(0, 0.7, 0.7, 0))
for (i in 1:8) {
  cbp = barplot(ChargaffTable[i, ], horiz = TRUE, axes = FALSE, axisnames = FALSE, col = mycolors[i])
  ax = axis(3, las = 2, labels = FALSE, col = mycolors[i], cex = 0.5, at = c(0, 10, 20))
  mtext(side = 3, at = ax,  text = paste(ax), col = mycolors[i], line = 0, las = 1, cex = 0.9)
  mtext(side = 2, at = cbp, text = colnames(ChargaffTable), col = mycolors[i], line = 0, las = 2, cex = 1)
  title(paste(rownames(ChargaffTable)[i]), col = mycolors[i], cex = 1.1)
}

## -----------------------------------------------------------------------------
#| label: fig-permstatChf-1
#| out-width: 50%
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "Histogram of our statistic `statChf` computed from simulations using per-row permutations of the columns. The value it yields for the observed data is shown by the red line."
statChf = function(x){
  sum((x[, "C"] - x[, "G"])^2 + (x[, "A"] - x[, "T"])^2)
}
chfstat = statChf(ChargaffTable)
permstat = replicate(100000, {
     permuted = t(apply(ChargaffTable, 1, sample))
     colnames(permuted) = colnames(ChargaffTable)
     statChf(permuted)
})
pChf = mean(permstat <= chfstat)
pChf
hist(permstat, breaks = 100, main = "", col = "lavender")
abline(v = chfstat, lwd = 2, col = "red")

## -----------------------------------------------------------------------------
#| label: vcdHC
HairEyeColor[,, "Female"]

## -----------------------------------------------------------------------------
#| label: answerHC
#| eval: !expr c(1)
str(HairEyeColor)
?HairEyeColor

## -----------------------------------------------------------------------------
#| label: Deuto
load("../data/Deuteranopia.RData")
Deuteranopia

## -----------------------------------------------------------------------------
#| label: chisq.test.Deuteranopia
chisq.test(Deuteranopia)

## -----------------------------------------------------------------------------
#| label: fig-HardyWeinberg-1
#| column: margin
#| fig-cap: !expr paste("Plot of the log-likelihood for the", Mourant$Country[216], "data.")
library("HardyWeinberg")
data("Mourant")
Mourant[214:216,]
nMM = Mourant$MM[216]
nMN = Mourant$MN[216]
nNN = Mourant$NN[216]
loglik = function(p, q = 1 - p) {
  2 * nMM * log(p) + nMN * log(2*p*q) + 2 * nNN * log(q)
}
xv = seq(0.01, 0.99, by = 0.01)
yv = loglik(xv)
plot(x = xv, y = yv, type = "l", lwd = 2,
     xlab = "p", ylab = "log-likelihood")
imax = which.max(yv)
abline(v = xv[imax], h = yv[imax], lwd = 1.5, col = "blue")
abline(h = yv[imax], lwd = 1.5, col = "purple")

## -----------------------------------------------------------------------------
#| label: phat
phat  =  af(c(nMM, nMN, nNN))
phat
pMM   =  phat^2
qhat  =  1 - phat

## -----------------------------------------------------------------------------
#| label: hweq
pHW = c(MM = phat^2, MN = 2*phat*qhat, NN = qhat^2)
sum(c(nMM, nMN, nNN)) * pHW

## -----------------------------------------------------------------------------
#| label: fig-HWtern
#| echo: !expr -c(1)
#| results: false
#| out-width: 50%
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "This **de Finetti plot** shows the points as barycenters of the three genotypes using the frequencies as weights on each of the corners of the triangle. The Hardy-Weinberg model is the red curve, the acceptance region is between the two purple lines. We see that the US is the furthest from being in HW equilibrium."
par(mai = rep(0.1, 4))
pops = c(1, 69, 128, 148, 192)
genotypeFrequencies = as.matrix(Mourant[, c("MM", "MN", "NN")])
HWTernaryPlot(genotypeFrequencies[pops, ],
        markerlab = Mourant$Country[pops],
        alpha = 0.0001, curvecols = c("red", rep("purple", 4)),
        mcex = 0.75, vertex.cex = 1)

## -----------------------------------------------------------------------------
#| label: quesTern
#| echo: !expr -c(1)
#| fig-show: hide
HWTernaryPlot(genotypeFrequencies[pops, ],
              markerlab = Mourant$Country[pops],
              curvecols = c("red", rep("purple", 4)),
              alpha = 0.0001, mcex = 0.75, vertex.cex = 1)
HWTernaryPlot(genotypeFrequencies[-pops, ], 
              newframe = FALSE, alpha = 0.0001, cex = 0.5)

## -----------------------------------------------------------------------------
#| label: newMNdata
#| fig-show: hide
newgf = round(genotypeFrequencies / 50)
HWTernaryPlot(newgf[pops, ],
              markerlab = Mourant$Country[pops],
              curvecols = c("red", rep("purple", 4)),
              alpha = 0.0001, mcex = 0.75, vertex.cex = 1)

## -----------------------------------------------------------------------------
#| label: fig-seqlogo-1
#| out-width: 75%
#| fig-height: 5
#| fig-width: 5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Here is a diagram called a sequence logo for the position dependent multinomial used to model the Kozak motif. It codifies the amount of variation in each of the positions on a log scale. The large letters represent positions where there is no uncertainty about which nucleotide occurs."
library("seqLogo")
load("../data/kozak.RData")
kozak
pwm = makePWM(kozak)
seqLogo(pwm, ic.scale = FALSE)

## -----------------------------------------------------------------------------
#| label: 4stateMC
#| echo: false
library("markovchain")
library("igraph")
sequence = toupper(c("a", "c", "a", "c", "g", "t", "t", "t", "t", "c", "c",
                     "a", "c", "g", "t", "a", "c","c","c","a","a","a","t","a",
                     "c","g","g","c","a","t","g","t","g","t","g","a","g","c","t","g"))
mcFit   =  markovchainFit(data = sequence)
MCgraph =  markovchain:::.getNet(mcFit$estimate, round = TRUE)
edgelab =  round(E(MCgraph)$weight / 100, 2)

## -----------------------------------------------------------------------------
#| label: fig-statsfourstateMC
#| echo: false
#| column: margin
#| fig-width: 4
#| fig-height: 3.5
#| fig-cap: "Visualisation of a 4-state Markov chain. The probability of each possible digram (e.,g., CA) is given by the weight of the edge between the corresponding nodes. So for instance, the probability of CA is given by the edge C$\to$ A. We'll see in @sec-images how to use R packages to draw these type of network graphs."
par(mai=c(0,0,0,0))
plot.igraph(MCgraph, edge.label = edgelab,
       vertex.size = 40, xlim = c(-1, 1.25))

## -----------------------------------------------------------------------------
#| label: fig-FreqBayes-turtles
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-STRDefinition
#| echo: false
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: fig-YSTRPositions
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-USY-STR
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: haplo6
haplo6 = read.table("../data/haplotype6.txt", header = TRUE)
haplo6

## -----------------------------------------------------------------------------
#| label: haplo6check
#| echo: false
with(haplo6, stopifnot(Individual[1] == "H1", DYS19[1] == 14, DXYS156Y[1] == 12))

## -----------------------------------------------------------------------------
#| label: fig-histobeta2
#| echo: false
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "Beta distributions with $\\alpha=10,20,50$ and $\\beta=30,60,150$. We can use these as a **prior** for probability of success in a binomial experiment. These three distributions have the same mean ($\\frac{\\alpha}{\\alpha +\\beta}$), but different concentrations around the mean."
dfbetas = data.frame(
  p = rep(p_seq, 3),
  dbeta = c(dbeta(p_seq,  10,  30),
            dbeta(p_seq,  20,  60), 
            dbeta(p_seq,  50, 150)),
  pars = rep(c("Beta(10,30)", "Beta(20,60)", "Beta(50,150)"), each = length(p_seq)))
library("ggplot2")
ggplot(dfbetas) +
  geom_line(aes(x = p, y = dbeta, colour = pars)) +
  theme(legend.title = element_blank()) +
  geom_vline(aes(xintercept = 0.25), colour = "#990000", linetype = "dashed")

## -----------------------------------------------------------------------------
#| label: fig-histmarginal-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "Marginal Distribution of $Y$."
rp = rbeta(100000, 50, 350)
y = vapply(rp, 
           function(x) rbinom(1, prob = x, size = 300), 
           integer(1))
hist(y, breaks = 50, col = "orange", main = "", xlab = "")

## -----------------------------------------------------------------------------
#| label: freqquesvectorize
set.seed(0xbebe)
y1 = vapply(rp, 
            function(x) rbinom(1, prob = x, size = 300), 
            integer(1))
set.seed(0xbebe)
y2 = rbinom(length(rp), rp, size = 300)
stopifnot(identical(y1, y2))

## -----------------------------------------------------------------------------
#| label: fig-densityposterior-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "Only choosing the values of the distribution with $Y=40$ gives the posterior distribution of $p$. The histogram (green) shows the simulated values for the posterior distribution, the line the density of a Beta distribution with the theoretical parameters."
pPostEmp = rp[ y == 40 ]
hist(pPostEmp, breaks = 40, col = "chartreuse4", main = "",
  probability = TRUE, xlab = "posterior p")

p_seq = seq(0, 1, by = 0.001)
densPostTheory = dbeta(p_seq, 50 + 40, 350 + 260)
lines(p_seq, densPostTheory, type = "l", lwd = 3)

## -----------------------------------------------------------------------------
#| label: comparetheory1
mean(pPostEmp)
dp = p_seq[2] - p_seq[1]
sum(p_seq * densPostTheory * dp)

## -----------------------------------------------------------------------------
#| label: comparetheory2
#| echo: false
stopifnot(abs(mean(pPostEmp) - sum(p_seq * densPostTheory * dp)) < 1e-3)

## -----------------------------------------------------------------------------
#| label: mcint
pPostMC = rbeta(n = 100000, 90, 610)
mean(pPostMC)

## -----------------------------------------------------------------------------
#| label: fig-qqplotbeta-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "Quantile-quantile (QQ) plot of our Monte Carlo sample `pPostMC` from the theoretical distribution and our simulation sample `pPostEmp`. We could also similarly compare either of these two distributions to the theoretical distribution function `pbeta(., 90, 610)`. If the curve lies on the line $y=x$, this indicates a good agreement. There are some random differences at the tails."
qqplot(pPostMC, pPostEmp, type = "l", asp = 1)
abline(a = 0, b = 1, col = "blue")

## -----------------------------------------------------------------------------
#| label: postbeta
densPost2 = dbeta(p_seq, 115, 735)
mcPost2   = rbeta(1e6, 115, 735)
sum(p_seq * densPost2 * dp)   # mean, by numeric integration
mean(mcPost2)                 # mean by MC
p_seq[which.max(densPost2)]   # MAP estimate

## -----------------------------------------------------------------------------
#| label: quantilespost
quantile(mcPost2, c(0.025, 0.975))

## -----------------------------------------------------------------------------
#| label: fig-DESeq2-Prediction-Interval
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: callBios
library("Biostrings")

## -----------------------------------------------------------------------------
#| label: BiostringExplore
#| eval: false
#| results: hide
## GENETIC_CODE
## IUPAC_CODE_MAP
## vignette(package = "Biostrings")
## vignette("BiostringsQuickOverview", package = "Biostrings")

## -----------------------------------------------------------------------------
#| label: BiostringCheck
#| echo: false
#| results: hide
GENETIC_CODE
IUPAC_CODE_MAP

## -----------------------------------------------------------------------------
#| label: BSgenomes
library("BSgenome")
ag = available.genomes()
length(ag)
ag[1:2]

## -----------------------------------------------------------------------------
#| label: BSGenomeEcoli
#| results: hide
library("BSgenome.Ecoli.NCBI.20080805")
Ecoli
shineDalgarno = "AGGAGGT"
ecoli = Ecoli$NC_010473

## -----------------------------------------------------------------------------
#| label: window
window = 50000
starts = seq(1, length(ecoli) - window, by = window)
ends   = starts + window - 1
numMatches = vapply(seq_along(starts), function(i) {
  countPattern(shineDalgarno, ecoli[starts[i]:ends[i]],
               max.mismatch = 0)
  }, numeric(1))
table(numMatches)

## -----------------------------------------------------------------------------
#| label: fig-poissonness
#| out-width: 50%
#| fig-height: 5
#| fig-cap: "Evaluation of a Poisson model for motif counts along the sequence `Ecoli$NC_010473`."
library("vcd")
gf = goodfit(numMatches, "poisson")
summary(gf)
distplot(numMatches, type = "poisson")

## -----------------------------------------------------------------------------
#| label: pattlocIranges1
#| results: hide
sdMatches = matchPattern(shineDalgarno, ecoli, max.mismatch = 0)

## -----------------------------------------------------------------------------
#| label: pattlocIranges2
betweenmotifs = gaps(sdMatches)

## -----------------------------------------------------------------------------
#| label: fig-expplotdata-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 4
#| fig-cap: "Evaluation of fit to the exponential distribution (black line) of the gaps between the motifs."
library("Renext")
expplot(width(betweenmotifs), rate = 1/mean(width(betweenmotifs)),
        labels = "fit")

## -----------------------------------------------------------------------------
#| label: gof
#| echo: false
#| eval: false
## gofExp.test(width(betweenmotifs))

## -----------------------------------------------------------------------------
#| label: chr8HS
#| cache: false
library("BSgenome.Hsapiens.UCSC.hg19")
chr8  =  Hsapiens$chr8
CpGtab = read.table("../data/model-based-cpg-islands-hg19.txt",
                    header = TRUE)
nrow(CpGtab)
head(CpGtab)
irCpG = with(dplyr::filter(CpGtab, chr == "chr8"),
         IRanges(start = start, end = end))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "We use the  `::` operator to call the  `filter` function specifically from the  `dplyr` package—and not from any other packages that may happen to be loaded and defining functions of the same name. This precaution is particularly advisable in the case of the  `filter` function, since this name is used by quite a few other packages. You can think of the normal (without  `::`) way of calling R functions like calling people by their first (given) names; whereas the fully qualified version with  `::` corresponds to calling someone by their full name. At least within the reach of the CRAN and Bioconductor repositories, such fully qualified names are guaranteed to be unique."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: grCpG
grCpG = GRanges(ranges = irCpG, seqnames = "chr8", strand = "+")
genome(grCpG) = "hg19"

## -----------------------------------------------------------------------------
#| label: fig-freqandbayes-ideo
#| column: margin
#| fig-height: 4
#| fig-cap: "**[Gviz](https://bioconductor.org/packages/Gviz/)** plot of CpG locations in a selected region of chromosome 8."
library("Gviz")
ideo = IdeogramTrack(genome = "hg19", chromosome = "chr8")
plotTracks(
  list(GenomeAxisTrack(),
    AnnotationTrack(grCpG, name = "CpG"), ideo),
    from = 2200000, to = 5800000,
    shape = "box", fill = "#006400", stacking = "dense")

## -----------------------------------------------------------------------------
#| label: CGIview
CGIview    = Views(unmasked(Hsapiens$chr8), irCpG)
NonCGIview = Views(unmasked(Hsapiens$chr8), gaps(irCpG))

## -----------------------------------------------------------------------------
#| label: CGIview2
seqCGI      = as(CGIview, "DNAStringSet")
seqNonCGI   = as(NonCGIview, "DNAStringSet")
dinucCpG    = sapply(seqCGI, dinucleotideFrequency)
dinucNonCpG = sapply(seqNonCGI, dinucleotideFrequency)
dinucNonCpG[, 1]
NonICounts = rowSums(dinucNonCpG)
IslCounts  = rowSums(dinucCpG)

## -----------------------------------------------------------------------------
#| label: transitions
TI  = matrix( IslCounts, ncol = 4, byrow = TRUE)
TnI = matrix(NonICounts, ncol = 4, byrow = TRUE)
dimnames(TI) = dimnames(TnI) =
  list(c("A", "C", "G", "T"), c("A", "C", "G", "T"))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "The transition probabilities are probabilities so the rows need to sum to 1."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: MI
MI = TI /rowSums(TI)
MI
MN = TnI / rowSums(TnI)
MN

## -----------------------------------------------------------------------------
#| label: STATI
freqIsl = alphabetFrequency(seqCGI, baseOnly = TRUE, collapse = TRUE)[1:4]
freqIsl / sum(freqIsl)
freqNon = alphabetFrequency(seqNonCGI, baseOnly = TRUE, collapse = TRUE)[1:4]
freqNon / sum(freqNon)

## -----------------------------------------------------------------------------
#| label: book-chunk-1
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: alphabeta
alpha = log((freqIsl/sum(freqIsl)) / (freqNon/sum(freqNon)))
beta  = log(MI / MN)

## -----------------------------------------------------------------------------
#| label: scorepatt
x = "ACGTTATACTACG"
scorefun = function(x) {
  s = unlist(strsplit(x, ""))
  score = alpha[s[1]]
  if (length(s) >= 2)
    for (j in 2:length(s))
      score = score + beta[s[j-1], s[j]]
  score
}
scorefun(x)

## -----------------------------------------------------------------------------
#| label: scorefun1
generateRandomScores = function(s, len = 100, B = 1000) {
  alphFreq = alphabetFrequency(s)
  isGoodSeq = rowSums(alphFreq[, 5:ncol(alphFreq)]) == 0
  s = s[isGoodSeq]
  slen = sapply(s, length)
  prob = pmax(slen - len, 0)
  prob = prob / sum(prob)
  idx  = sample(length(s), B, replace = TRUE, prob = prob)
  ssmp = s[idx]
  start = sapply(ssmp, function(x) sample(length(x) - len, 1))
  scores = sapply(seq_len(B), function(i)
    scorefun(as.character(ssmp[[i]][start[i]+(1:len)]))
  )
  scores / len
}
scoresCGI    = generateRandomScores(seqCGI)
scoresNonCGI = generateRandomScores(seqNonCGI)

## -----------------------------------------------------------------------------
#| label: fig-ScoreMixture-1
#| column: margin
#| fig-height: 5
#| fig-cap: "Island and non-island scores as generated by the function `generateRandomScores`. This is the first instance of a **mixture** we encounter. We will revisit them in @sec-mixtures."
rgs = range(c(scoresCGI, scoresNonCGI))
br = seq(rgs[1], rgs[2], length.out = 50)
h1 = hist(scoresCGI,    breaks = br, plot = FALSE)
h2 = hist(scoresNonCGI, breaks = br, plot = FALSE)
plot(h1, col = rgb(0, 0, 1, 1/4), xlim = c(-0.5, 0.5), ylim=c(0,120))
plot(h2, col = rgb(1, 0, 0, 1/4), add = TRUE)

## -----------------------------------------------------------------------------
#| label: savescoresforChap4
#| echo: false
#| eval: false
## ###This is for provenance reasons, keep track of how the data
## ###were generated for the EM exercise in Chapter 4.
## Mdata=c(scoresCGI,scoresNonCGI)
## MM1=sample(Mdata[1:1000],800)
## MM2=sample(Mdata[1001:2000],1000)
## Myst=c(MM1,MM2);names(Myst)=NULL
## saveRDS(c(MM1,MM2),"../data/Myst.rds")
## ###True value of m1,m2,s1 and s2
## ###

## -----------------------------------------------------------------------------
#| label: checkhists
#| echo: false
stopifnot(max(h1$counts) < 120, max(h2$counts) < 120,
          h1$breaks[1] >= br[1], h1$breaks[length(h1$breaks)] <= br[length(br)],
          h2$breaks[1] >= br[1], h2$breaks[length(h2$breaks)] <= br[length(br)])

## -----------------------------------------------------------------------------
#| label: mtbotub
mtb = read.table("../data/M_tuberculosis.txt", header = TRUE)
head(mtb, n = 4)

## -----------------------------------------------------------------------------
#| label: ProlMyc
pro  =  mtb[ mtb$AmAcid == "Pro", "Number"]
pro/sum(pro)

## -----------------------------------------------------------------------------
#| label: staphread
staph = readDNAStringSet("../data/staphsequence.ffn.txt", "fasta")

## -----------------------------------------------------------------------------
#| label: staphex
staph[1:3, ]
staph

## -----------------------------------------------------------------------------
#| label: GCfreq
letterFrequency(staph[[1]], letters = "ACGT", OR = 0)
GCstaph = data.frame(
  ID = names(staph),
  GC = rowSums(alphabetFrequency(staph)[, 2:3] / width(staph)) * 100
)

## -----------------------------------------------------------------------------
#| label: fig-SlidingGC-1
#| out-width: 75%
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "GC content along sequence 364 of the *Staphylococcus Aureus* genome."
window = 100
gc = rowSums( letterFrequencyInSlidingView(staph[[364]], window,
      c("G","C")))/window
plot(x = seq(along = gc), y = gc, type = "l")

## -----------------------------------------------------------------------------
#| label: fig-SmoothSlidingGC-1
#| out-width: 75%
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "Similar to @fig-SlidingGC-1, with smoothing."
plot(x = seq(along = gc), y = gc, type = "l")
lines(lowess(x = seq(along = gc), y = gc, f = 0.2), col = 2)

## -----------------------------------------------------------------------------
#| label: fig-histobeta4-1
#| out-width: 50%
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "Beta densities for different parameter choices."
dfbetas = data.frame(
  p = rep(p_seq, 5),
  dbeta = c(dbeta(p_seq, 0.5, 0.5), 
            dbeta(p_seq,   1,   1), 
            dbeta(p_seq,  10,  30),
            dbeta(p_seq,  20,  60), 
            dbeta(p_seq,  50, 150)),
  pars = rep(c("Beta(0.5,0.5)", "U(0,1)=Beta(1,1)", 
               "Beta(10,30)", "Beta(20,60)", 
               "Beta(50,150)"), each = length(p_seq)))
ggplot(dfbetas) +
  geom_line(aes(x = p, y = dbeta, colour = pars)) +
  theme(legend.title = element_blank()) +
  geom_vline(aes(xintercept = 0.25), colour = "#990000", linetype = "dashed")
