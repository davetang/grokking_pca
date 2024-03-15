
## -----------------------------------------------------------------------------
#| label: Pile-ou-face
#| echo: false

## -----------------------------------------------------------------------------
#| label: dpois1
dpois(x = 3, lambda = 5)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Note how the output from R is formatted: the first line begins with the first item in the vector, hence the [1], and the second line begins with the 9th item, hence the [9]. This helps you keep track of elements in long vectors. The term *vector* is R parlance for an ordered list of elements of the same type (in this case, numbers)."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-Poisson5
#| echo: !expr -c(1,5)
#| column: margin
#| fig-height: 4
#| fig-cap: "Probabilities of seeing 0,1,2,...,12 mutations, as modeled by the Poisson(5) distribution. The plot shows that we will often see 4 or 5 mutations but rarely as many as 12. The distribution continues to higher numbers ($13,...$), but the probabilities will be successively smaller, and here we don't visualize them."
.oldopt = options(digits = 2)
0:12
dpois(x = 0:12, lambda = 5)
barplot(dpois(0:12, 5), names.arg = 0:12, col = "red")
options(.oldopt)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Think of a categorical variable as having different alternative values. These are the levels, similar to the different alternatives at a gene locus: *alleles*."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "`c()` is one of the most basic functions. It collates elements of the same type into a vector. In the code shown here, the elements of `genotype` are character strings."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: genotype1
genotype = c("AA","AO","BB","AO","OO","AO","AA","BO","BO",
             "AO","BB","AO","BO","AB","OO","AB","BB","AO","AO")
table(genotype)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "It is not obvious from the output of the `table` function that the input was a factor; however if there had been another level with no instances, the table would also have contained that level, with a zero count."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: genotype2
genotypeF = factor(genotype)
levels(genotypeF)
table(genotypeF)

## -----------------------------------------------------------------------------
#| label: fig-twoballs
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: rbinom1
rbinom(15, prob = 0.5, size = 1)

## -----------------------------------------------------------------------------
#| label: rbinom2
rbinom(12, prob = 2/3, size = 1)

## -----------------------------------------------------------------------------
#| label: rbinom3
rbinom(1, prob = 2/3, size = 12)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "What does `set.seed` do here?"

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: rbinom4
set.seed(235569515)
rbinom(1, prob = 0.3, size = 15)

## -----------------------------------------------------------------------------
#| label: dbinom
probabilities = dbinom(0:15, prob = 0.3, size = 15)
round(probabilities, 2)

## -----------------------------------------------------------------------------
#| label: fig-binombarplot
#| column: margin
#| fig-width:  7
#| fig-height: 6
#| fig-cap: !expr sprintf("Theoretical distribution of $B(15,0.3)$ . The highest bar is at $x=%s$. We have chosen to represent theoretical values in red throughout.", .freqoutcome)
barplot(probabilities, names.arg = 0:15, col = "red")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Instead of $\\frac{n!}{(n-k)!k!}$ we can use the special notation ${n \\choose k}$ as a shortcut."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-Simeon-Poisson
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: testbinomvspois
#| echo: false
#| eval: false
## plot(dbinom(0:12, prob = 5e-4, size = 1e4),
##      dpois(0:12, lambda = 5), asp = 1)
## abline(a = 0, b = 1, col = "blue")

## -----------------------------------------------------------------------------
#| label: dpois2
5^3 * exp(-5) / factorial(3)

## -----------------------------------------------------------------------------
#| label: fig-gen-simpoisson
#| column: margin
#| fig-width:  7
#| fig-height: 5.5
#| fig-cap: !expr paste("Simulated distribution of B(10000, $10^{-4}$) for", length(simulations), "simulations.")
rbinom(1, prob = 5e-4, size = 10000)
simulations = rbinom(n = 300000, prob = 5e-4, size = 10000)
barplot(table(simulations), col = "lavender")

## -----------------------------------------------------------------------------
#| label: fig-antibody
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: brack100
#| echo: false
`[<-`(rep(0, 100), 22, 1)

## -----------------------------------------------------------------------------
#| label: fig-typicalP
#| echo: false
#| column: margin
#| fig-width:  6
#| fig-height: 4
#| fig-cap: "Plot of typical data from our generative model for the background, i.\\,e., for the false positive hits: 100 positions along the protein, at each position the count is drawn from a Poisson(0.5) random variable."
s100 = rpois(100, lambda=0.5)
barplot(s100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5,100.5),
  names.arg = seq(along = s100), col="lavender")

## -----------------------------------------------------------------------------
#| label: For_the_record
#| eval: false
#| echo: false
## set.seed(8969311)
## e100 = rpois(100,lambda = 0.5)
## e100[42] = 7
## save(e100, file = "../data/e100.RData")

## -----------------------------------------------------------------------------
#| label: epitopeData
#| fig-show: hide
load("../data/e100.RData")
barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
  names.arg = seq(along = e100), col = "darkolivegreen")

## -----------------------------------------------------------------------------
#| label: fig-epitopedata
#| echo: false
#| column: margin
#| fig-width:  6
#| fig-height: 4
#| fig-cap: "Output of the ELISA array results for 50 patients in the 100 positions."
barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
  names.arg = seq(along = e100), col = "darkolivegreen")
text(35, 7, adj = c(-0.05, 0.5), labels = "?", xpd = NA, col = "red",
  cex = 1.25, font = 2)

## -----------------------------------------------------------------------------
#| label: ppois
1 - ppois(6, 0.5)
ppois(6, 0.5, lower.tail = FALSE)

## -----------------------------------------------------------------------------
#| label: montecarlomethod
maxes = replicate(100000, {
  max(rpois(100, 0.5))
})
table(maxes)

## -----------------------------------------------------------------------------
#| label: meanmaxes
mean( maxes >= 7 )

## -----------------------------------------------------------------------------
#| label: fig-Boxes4
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "You are secretly meeting a continuous distribution here, the uniform distribution: `runif`."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: dmultinom
dmultinom(c(4, 2, 0, 0), prob = rep(1/4, 4))

## -----------------------------------------------------------------------------
#| label: pvec
pvec = rep(1/4, 4)
t(rmultinom(1, prob = pvec, size = 8))

## -----------------------------------------------------------------------------
#| label: SampleSize
#| echo: false

## -----------------------------------------------------------------------------
#| label: obsunder0start
obsunder0 = rmultinom(1000, prob = pvec, size = 20)
dim(obsunder0)
obsunder0[, 1:11]

## -----------------------------------------------------------------------------
#| label: assertpvec
#| echo: false
thep = unique(pvec); stopifnot(length(thep)==1, thep == 0.25)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Notice that the top of every column there is an index of the form `[,1][,2]...` These are the column indices. The rows are labeled `[1,][2,]...`. The object `obsunder0` is not a simple vector as those we have seen before, but an array of numbers in a matrix."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: obsunder0look
expected0 = pvec * 20
sum((obsunder0[, 1] - expected0)^2 / expected0)
sum((obsunder0[, 2] - expected0)^2 / expected0)
sum((obsunder0[, 3] - expected0)^2 / expected0)

## -----------------------------------------------------------------------------
#| label: stat
stat = function(obsvd, exptd = 20 * pvec) {
  sum((obsvd - exptd)^2 / exptd)
}
stat(obsunder0[, 1])

## -----------------------------------------------------------------------------
#| label: fig-histS0
#| column: margin
#| fig-cap: "The histogram of simulated values `S0` of the statistic `stat` under the null (fair) distribution provides an approximation of the **sampling distribution** of the statistic `stat`."
S0 = apply(obsunder0, 2, stat)
summary(S0)
hist(S0, breaks = 25, col = "lavender", main = "")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "The `apply` function is shorthand for a loop over the rows or columns of an array. Here the second argument, 2, indicates looping over the columns."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: quantile
q95 = quantile(S0, probs = 0.95)
q95

## -----------------------------------------------------------------------------
#| label: saveS0
#| echo: false
#| eval: false
## ## This was done to save this object for its reuse in Chapter 2.
## save(S0, file = "../data/S0.RData")

## -----------------------------------------------------------------------------
#| label: roulette-chunk-1
#| echo: false

## -----------------------------------------------------------------------------
#| label: alternativeA
pvecA = c(3/8, 1/4, 1/4, 1/8)
observed = rmultinom(1000, prob = pvecA, size = 20)
dim(observed)
observed[, 1:7]
apply(observed, 1, mean)
expectedA = pvecA * 20
expectedA

## -----------------------------------------------------------------------------
#| label: computestat1
stat(observed[, 1])
S1 = apply(observed, 2, stat)
q95
sum(S1 > q95)
power = mean(S1 > q95)
power

## -----------------------------------------------------------------------------
#| label: assertS1
#| echo: false
stopifnot(stat(observed[, 1]) < q95)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "We read the vertical line as **given** or **conditional on**."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-ProbaDiagram
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "$P(H_0\\;|\\;\\text{data})$ is not the same as a p-value $P(\\text{data}\\;|\\;H_0)$."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
dbinom(2, size = 10, prob = 0.3)
pbinom(2, size = 10, prob = 0.3)
sum(dbinom(0:2, size = 10, prob = 0.3)) 

## -----------------------------------------------------------------------------
#| label: poismax1
poismax = function(lambda, n, m) {
  epsilon = 1 - ppois(m - 1, lambda)
  1 - exp( -n * epsilon)
}
poismax(lambda = 0.5, n = 100, m = 7)
poismax(lambda = mean(e100), n = 100, m = 7)

## -----------------------------------------------------------------------------
#| label: poismax2
poismax = function(lambda, n = 100, m = 7) {
  1 - exp( -n * (1 - ppois(m - 1, lambda)))
}
poismax(0.5)
poismax(0.5, m = 9)

## -----------------------------------------------------------------------------
#| label: biocex1
#| eval: false
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install(c("Biostrings", "BSgenome.Celegans.UCSC.ce2"))

## -----------------------------------------------------------------------------
#| label: BSgenome.Celegans.UCSC.ce2
library("BSgenome.Celegans.UCSC.ce2")
Celegans
seqnames(Celegans)
Celegans$chrM
class(Celegans$chrM)
length(Celegans$chrM)

## -----------------------------------------------------------------------------
#| label: Biostr
library("Biostrings")
lfM = letterFrequency(Celegans$chrM, letters=c("A", "C", "G", "T"))
lfM
sum(lfM)
lfM / sum(lfM)

## -----------------------------------------------------------------------------
#| label: rf
t(rmultinom(1, length(Celegans$chrM), p = rep(1/4, 4)))

## -----------------------------------------------------------------------------
#| label: lengthM
length(Celegans$chrM) / 4

## -----------------------------------------------------------------------------
#| label: oestat
oestat = function(o, e) {
  sum((o-e)^2 / e)
}
oe = oestat(o = lfM, e = length(Celegans$chrM) / 4)
oe

## -----------------------------------------------------------------------------
#| label: oesim
B = 10000
n = length(Celegans$chrM)
expected = rep(n / 4, 4)
oenull = replicate(B,
  oestat(e = expected, o = rmultinom(1, n, p = rep(1/4, 4))))

## -----------------------------------------------------------------------------
#| label: resundernull
#| eval: false
#| echo: false
## hist(oenull, breaks = 100, col = "skyblue", main = "")

## -----------------------------------------------------------------------------
#| label: assert
#| echo: false
stopifnot( oe/10 > max(oenull) )
