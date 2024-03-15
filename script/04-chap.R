
## -----------------------------------------------------------------------------
#| label: t-distribution
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-twocoins
#| column: margin
#| fig-asp: 1
#| fig-height: 3.5
#| fig-cap: "Histogram of 10,000 random draws from a fair mixture of two normals. The left hand part of the histogram is dominated by numbers generated from (A), on the right from (B)."
coinflips = (runif(10000) > 0.5)
table(coinflips)
oneFlip = function(fl, mean1 = 1, mean2 = 3, sd1 = 0.5, sd2 = 0.5) {
  if (fl) {
   rnorm(1, mean1, sd1)
  } else {
   rnorm(1, mean2, sd2)
  }
}
fairmix = vapply(coinflips, oneFlip, numeric(1))
library("ggplot2")
library("dplyr")
ggplot(tibble(value = fairmix), aes(x = value)) +
     geom_histogram(fill = "purple", binwidth = 0.1)

## -----------------------------------------------------------------------------
#| label: efficient
means = c(1, 3)
sds   = c(0.5, 0.5)
values = rnorm(length(coinflips),
               mean = ifelse(coinflips, means[1], means[2]),
               sd   = ifelse(coinflips, sds[1],   sds[2]))

## -----------------------------------------------------------------------------
#| label: fig-limitinghistogram
#| fig-asp: 1
#| fig-height: 3
#| out-width: 33%
#| fig-cap: "Similar to @fig-twocoins, but with one million observations."
fair = tibble(
  coinflips = (runif(1e6) > 0.5),
  values = rnorm(length(coinflips),
                 mean = ifelse(coinflips, means[1], means[2]),
                 sd   = ifelse(coinflips, sds[1],   sds[2])))
ggplot(fair, aes(x = values)) +
     geom_histogram(fill = "purple", bins = 200)

## -----------------------------------------------------------------------------
#| label: fig-overlaydensity
#| fig-asp: 1
#| fig-height: 3
#| out-width: 33%
#| fig-cap: !expr sprintf("Histogram of half a million draws from the normal distribution $N(\\mu=%s,\\sigma^2=%s^2)$. The curve is the theoretical density $\\phi(x)$ calculated using the function `dnorm`.", means[1], sds[1])
ggplot(dplyr::filter(fair, coinflips), aes(x = values)) +
  geom_histogram(aes(y = after_stat(density)), fill = "purple", binwidth = 0.01) +
  stat_function(fun = dnorm, color = "red",
                args = list(mean = means[1], sd = sds[1]))

## -----------------------------------------------------------------------------
#| label: fig-twodensity
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "The theoretical density of the mixture."
fairtheory = tibble(
  x = seq(-1, 5, length.out = 1000),
  f = 0.5 * dnorm(x, mean = means[1], sd = sds[1]) +
      0.5 * dnorm(x, mean = means[2], sd = sds[2]))
ggplot(fairtheory, aes(x = x, y = f)) +
  geom_line(color = "red", linewidth = 1.5) + ylab("mixture density")

## -----------------------------------------------------------------------------
#| label: fig-histmystery
#| echo: false
#| results: hide
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "A mixture of two normals that is harder to recognize."
mystery = tibble(
  coinflips = (runif(1e3) > 0.5),
  values = rnorm(length(coinflips),
                 mean = ifelse(coinflips, 1, 2),
                 sd   = ifelse(coinflips, sqrt(.5), sqrt(.5))))
br2 = with(mystery, seq(min(values), max(values), length.out = 30))
ggplot(mystery, aes(x = values)) +
geom_histogram(fill = "purple", breaks = br2)

## -----------------------------------------------------------------------------
#| label: fig-betterhistogram-1
#| fig-asp: 1
#| fig-height: 3
#| out-width: 33%
#| fig-cap: "The mixture from @fig-histmystery, but with the two components colored in red and blue."
head(mystery, 3)
br = with(mystery, seq(min(values), max(values), length.out = 30))
ggplot(mystery, aes(x = values)) +
  geom_histogram(data = dplyr::filter(mystery, coinflips),
     fill = "red", alpha = 0.2, breaks = br) +
  geom_histogram(data = dplyr::filter(mystery, !coinflips),
     fill = "darkblue", alpha = 0.2, breaks = br) 

## -----------------------------------------------------------------------------
#| label: checkbr
#| echo: false
stopifnot(identical(br2, br))

## -----------------------------------------------------------------------------
#| label: fig-comparecomponents-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "As @fig-betterhistogram-1, with stacked bars for the two mixture components."
ggplot(mystery, aes(x = values, fill = coinflips)) +
  geom_histogram(data = dplyr::filter(mystery, coinflips),
     fill = "red", alpha = 0.2, breaks = br) +
  geom_histogram(data = dplyr::filter(mystery, !coinflips),
     fill = "darkblue", alpha = 0.2, breaks = br) +
  geom_histogram(fill = "purple", breaks = br, alpha = 0.2)

## -----------------------------------------------------------------------------
#| label: book-chunk-1
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: mixnorm1
mus = c(-0.5, 1.5)
lambda = 0.5
u = sample(2, size = 100, replace = TRUE, prob = c(lambda, 1-lambda))
x = rnorm(length(u), mean = mus[u])
dux = tibble(u, x)
head(dux)

## -----------------------------------------------------------------------------
#| label: mixnorm2
group_by(dux, u) |> summarize(mu = mean(x), sigma = sd(x))
table(dux$u) / nrow(dux)

## -----------------------------------------------------------------------------
#| label: mixtools
#| echo: !expr c(2:5)
.o = options(digits = 3)
library("mixtools")
y = c(rnorm(100, mean = -0.2, sd = 0.5),
      rnorm( 50, mean =  0.5, sd =   1))
gm = normalmixEM(y, k = 2, 
                    lambda = c(0.5, 0.5),
                    mu = c(-0.01, 0.01), 
                    sigma = c(3, 3))
with(gm, c(lambda, mu, sigma, loglik))
options(.o)

## -----------------------------------------------------------------------------
#| label: mosaics
#| echo: false
#| results: hide
## PROVENANCE: here's a record of how the data were created
library("mosaics")
library("mosaicsExample")
for(f in c("wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam",
           "wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam"))
  constructBins(infile = system.file(file.path("extdata", f), package="mosaicsExample"),
    fileFormat = "bam", outfileLoc = "../data/",
    byChr = FALSE, useChrfile = FALSE, chrfile = NULL, excludeChr = NULL,
    PET = FALSE, fragLen = 200, binSize = 200, capping = 0)

datafiles = c("../data/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt",
              "../data/wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt")
binTFBS = readBins(type = c("chip", "input"), fileName = datafiles)

## -----------------------------------------------------------------------------
#| label: binTFBS
binTFBS

## -----------------------------------------------------------------------------
#| label: fig-chipseqzeros
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "The number of binding sites found in 200nt windows along chromosome 22 in a ChIP-Seq dataset."
bincts = print(binTFBS)
ggplot(bincts, aes(x = tagCount)) +
  geom_histogram(binwidth = 1, fill = "forestgreen")

## -----------------------------------------------------------------------------
#| label: fig-ChipseqHistlogY
#| warning.known: !expr c("Transformation introduced infinite values in continuous y-axis", "Removed .. rows containing missing values")
#| fig-asp: 1
#| fig-height: 3
#| out-width: 33%
#| fig-cap: "As @fig-chipseqzeros, but using a logarithm base 10 scale on the $y$-axis. The fraction of zeros seems elevated compared to that of ones, twos, ..."
ggplot(bincts, aes(x = tagCount)) + scale_y_log10() +
   geom_histogram(binwidth = 1, fill = "forestgreen")

## -----------------------------------------------------------------------------
#| label: fig-nucleotideweights-1
#| column: margin
#| fig-height: 3
#| fig-cap: "Simulation of 7,000 nucleotide mass measurements."
masses = c(A =  331, C =  307, G =  347, T =  322)
probs  = c(A = 0.12, C = 0.38, G = 0.36, T = 0.14)
N  = 7000
sd = 3
nuclt   = sample(length(probs), N, replace = TRUE, prob = probs)
quadwts = rnorm(length(nuclt),
                mean = masses[nuclt],
                sd   = sd)
ggplot(tibble(quadwts = quadwts), aes(x = quadwts)) +
  geom_histogram(bins = 100, fill = "purple")

## -----------------------------------------------------------------------------
#| label: fig-ecdfZea
#| column: margin
#| fig-height: 2
#| fig-width: 4
#| fig-cap: "The observed sample can be seen as a mixture of point masses at each of the values (real point masses would be bars without any width whatsoever)."
library("HistData")
ZeaMays$diff
ggplot(ZeaMays, aes(x = diff, ymax = 1/15, ymin = 0)) +
  geom_linerange(linewidth = 1, col = "forestgreen") + ylim(0, 0.1)

## -----------------------------------------------------------------------------
#| label: checkZeaMays
#| echo: false
stopifnot(nrow(ZeaMays) == 15)

## -----------------------------------------------------------------------------
#| label: fig-samplingdist
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-bootpple
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-bootmedian
#| column: margin
#| fig-height: 2.5
#| fig-width: 4
#| fig-cap: "The bootstrap sampling distribution of the median of the Zea Mays differences."
B = 1000
meds = replicate(B, {
  i = sample(15, 15, replace = TRUE)
  median(ZeaMays$diff[i])
})
ggplot(tibble(medians = meds), aes(x = medians)) +
  geom_histogram(bins = 30, fill = "purple")

## -----------------------------------------------------------------------------
#| label: bootpkg
#| results: hide
library("bootstrap")
bootstrap(ZeaMays$diff, B, mean)
bootstrap(ZeaMays$diff, B, median)

## -----------------------------------------------------------------------------
#| child: devil.qmd

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
c(N3 = choose(5, 3), N15 = choose(29, 15))

## -----------------------------------------------------------------------------
#| label: fig-LaplacePortrait
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: Laplace1
w = rexp(10000, rate = 1)

## -----------------------------------------------------------------------------
#| label: fig-Laplacedistribution
#| column: margin
#| fig-width: 3
#| fig-height: 2.4
#| fig-cap: "Data sampled from a Laplace distribution."
mu  = 0.3
lps = rnorm(length(w), mean = mu, sd = sqrt(w))
ggplot(data.frame(lps), aes(x = lps)) +
  geom_histogram(fill = "purple", binwidth = 0.1)

## -----------------------------------------------------------------------------
#| label: fig-ALaplacedistribution
#| column: margin
#| fig-width: 3
#| fig-height: 2.4
#| fig-cap: "Histogram of data generated from an asymmetric Laplace distribution -- a scale mixture of many normals whose means and variances are dependent. We write $X \\sim AL(\\theta, \\mu, \\sigma)$."
mu = 0.3; sigma = 0.4; theta = -1
w  = rexp(10000, 1)
alps = rnorm(length(w), theta + mu * w, sigma * sqrt(w))
ggplot(tibble(alps), aes(x = alps)) +
  geom_histogram(fill = "purple", binwidth = 0.1)

## -----------------------------------------------------------------------------
#| label: fig-promoterlengthsandexpression
#| echo: false
#| column: margin
#| layout-ncol: 1
#| fig-width: 3
#| fig-height: 3
#| fig-margin: false
#| fig-cap: "Histogram of real data. Both distributions can be modeled by asymmetric Laplace distributions."
#| fig-subcap:
#|   - "The lengths of the promoters shorter than 2000bp from Saccharomyces cerevisiae as studied by  @Kristiansson2009."
    c('LaplaceMixturePromoterLengths.png', 'tcellhist.png')))

## -----------------------------------------------------------------------------
#| label: fig-three-worlds-web
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-gammahist1
#| column: margin
#| layout-ncol: 1
#| fig-height: 2.5
#| fig-subcap:
#|   - "gamma$(2,\\frac{1}{3})$"
#|   - "gamma$(10,\\frac{3}{2})$"
#| fig-cap: "Histograms of random samples of gamma distributions. The gamma is a flexible two parameter distribution: [see Wikipedia](http://en.wikipedia.org/wiki/Gamma_distribution)."
ggplot(tibble(x = rgamma(10000, shape = 2, rate = 1/3)),
   aes(x = x)) + geom_histogram(bins = 100, fill= "purple")
ggplot(tibble(x = rgamma(10000, shape = 10, rate = 3/2)),
   aes(x = x)) + geom_histogram(bins = 100, fill= "purple")

## -----------------------------------------------------------------------------
#| label: fig-generatepoissongamma
#| column: margin
#| fig-height: 2.5
#| fig-cap: "Histogram of `gp`, generated via a gamma-Poisson hierachical model."
lambda = rgamma(10000, shape = 10, rate = 3/2)
gp = rpois(length(lambda), lambda = lambda)
ggplot(tibble(x = gp), aes(x = x)) +
  geom_histogram(bins = 100, fill= "purple")

## -----------------------------------------------------------------------------
#| label: fig-goofy
#| out-width: 75%
#| fig-width: 8
#| fig-height: 4
#| fig-cap: "Goodness of fit plot. The **rootogram** shows the theoretical probabilities of the gamma-Poisson distribution (a.k.a. negative binomial) as red dots and the square roots of the observed frequencies as the height of the rectangular bars. The bars all end close to the horizontal axis, which indicates a good fit to the negative binomial distribution."
library("vcd")
ofit = goodfit(gp, "nbinomial")
plot(ofit, xlab = "")
ofit$par

## -----------------------------------------------------------------------------
#| label: setupgammapois
#| echo: false
x    = 0:95
mu   = 50
vtot = 80
v1   = vtot - mu
scale = v1/mu    # 0.6
shape = mu^2/v1  # 83.3
p1   = dgamma(x = x, scale = 0.6, shape = 80)
p2   = dpois(x = x, lambda = mu*1.2)
p3   = dnbinom(x = x, mu = mu, size = mu^2/vtot)

## -----------------------------------------------------------------------------
#| label: fig-mixtures-dgammapois
#| echo: false
#| out-width: 75%
#| fig-height: 6
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: !expr sprintf("Visualization of the hierarchical model that generates the gamma-Poisson distribution. The top panel shows the density of a gamma distribution with mean %s (vertical black line) and variance %s. Assume that in one particular experimental replicate, the value %s is realized. This is our latent variable. The observable outcome is distributed according to the Poisson distribution with that rate parameter, shown in the middle panel. In one particular experiment the outcome may be, say, %s, indicated by the dashed green line. Overall, if we repeat these two subsequent random process many times, the outcomes will be distributed as shown in the bottom panel -- the gamma-Poisson distribution.", mu, v1, mu*1.2, mu*1.1)
library("RColorBrewer")
cols = brewer.pal(8, "Paired")
par(mfrow=c(3,1), mai=c(0.5, 0.5, 0.01, 0.01))
xlim = x[c(1, length(x))]
plot(NA, NA, xlim=xlim, ylim=c(0,0.07), type="n", ylab="", xlab="")
polygon(x, p1, col=cols[1])
abline(v=mu, col="black", lwd=3)
abline(v=mu*1.2, col=cols[2], lty=2, lwd=3)
plot(x, p2, col=cols[3], xlim=xlim, ylab="", xlab="", type="h", lwd=2)
abline(v=mu*1.2, col=cols[2], lwd=2)
abline(v=mu*1.1, col=cols[4], lty=2, lwd=3)
plot(x, p3, col=cols[4], xlim=xlim, type="h", lwd=2, ylab="", xlab="")

## -----------------------------------------------------------------------------
#| label: fig-seriesofpoisson
#| column: margin
#| fig-width: 3
#| fig-height: 6
#| fig-cap: "Poisson distributed measurement data, for eight different choices of the mean `lambda`. In the upper panel, the $y$-axis is proportional to the data; in the lower panel, it is on a square-root scale. Note how the distribution widths change in the first case, but less so in the second."
simdat = lapply(seq(10, 100, by = 10), function(lam)
    tibble(n = rpois(200, lambda = lam),
           `sqrt(n)` = sqrt(n),
	   lambda = lam)) |>
  bind_rows() |>
  tidyr::pivot_longer(cols = !lambda)
ggplot(simdat, aes(x = factor(lambda), y = value)) + xlab(expression(lambda)) +
  geom_violin() + facet_grid(rows = vars(name), scales = "free")

## -----------------------------------------------------------------------------
#| label: vstpois
#| echo: !expr -c(1,3)
.o = options(digits = 3)
summarise(group_by(simdat, name, lambda), sd(value)) |> tidyr::pivot_wider(values_from = `sd(value)`)
options(.o)

## -----------------------------------------------------------------------------
#| label: fig-seriesofnb
#| column: margin
#| fig-width: 3
#| fig-height: 3
#| fig-cap: !expr sprintf("gamma-Poisson distributed measurement data, for a range of $\\mu$ from %s to %s.", dplyr::first(muvalues), dplyr::last(muvalues))
muvalues = 2^seq(0, 10, by = 1)
simgp = lapply(muvalues, function(mu) {
  u = rnbinom(n = 1e4, mu = mu, size = 4)
  tibble(mean = mean(u), sd = sd(u),
         lower = quantile(u, 0.025),
         upper = quantile(u, 0.975),
         mu = mu)
  } ) |> bind_rows()
head(as.data.frame(simgp), 2)
ggplot(simgp, aes(x = mu, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + geom_errorbar()

## -----------------------------------------------------------------------------
#| label: fig-pcwlin-1
#| column: margin
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "Piecewise linear function that stabilizes the variance of the data in @fig-seriesofnb."
simgp = mutate(simgp,
  slopes = 1 / sd,
  trsf   = cumsum(slopes * mean))
ggplot(simgp, aes(x = mean, y = trsf)) +
  geom_point() + geom_line() + xlab("")

## -----------------------------------------------------------------------------
#| label: fig-plotvstgammapoisson-1
#| out-width: 50%
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "Graph of the function @eq-mixtures-vstgammapoisson for different choices of $\\alpha$."
f = function(x, a) 
  ifelse (a==0, 
    sqrt(x), 
    log(2*sqrt(a) * sqrt(x*(a*x+1)) + 2*a*x+1) / (2*sqrt(a)))
x  = seq(0, 24, by = 0.1)
df = lapply(c(0, 0.05*2^(0:5)), function(a) 
  tibble(x = x, a = a, y = f(x, a))) %>% bind_rows()
ggplot(df, aes(x = x, y = y, col = factor(a))) + 
  geom_line() + labs(col = expression(alpha))

## -----------------------------------------------------------------------------
#| label: checkmaths
#| echo: !expr -c(3)
f2 = function(x, a) ifelse (a==0, sqrt(x), acosh(2*a*x + 1) / (2*sqrt(a)))  
with(df, max(abs(f2(x,a) - y)))
stopifnot(with(df, max(abs(f2(x,a) - y))) < 1e-10)

## -----------------------------------------------------------------------------
#| label: testasymptotics
  a = c(0.2, 0.5, 1)
  f(1e6, a) 
  1/(2*sqrt(a)) * (log(1e6) + log(4*a))

## -----------------------------------------------------------------------------
#| label: fig-EMillustrate-1
#| column: margin
#| fig-width: 2.5
#| fig-height: 2.5
#| fig-cap: "Histogram of `mx`, our example data for the EM algorithm."
mx = readRDS("../data/Myst.rds")$yvar
str(mx)
ggplot(tibble(mx), aes(x = mx)) + geom_histogram(binwidth = 0.025)

## -----------------------------------------------------------------------------
#| label: EM1
wA = runif(length(mx))
wB = 1 - wA

## -----------------------------------------------------------------------------
#| label: EM2
iter      = 0
loglik    = -Inf
delta     = +Inf
tolerance = 1e-12
miniter   = 50
maxiter   = 1000

## -----------------------------------------------------------------------------
#| label: EM3
while((delta > tolerance) && (iter <= maxiter) || (iter < miniter)) {
  lambda = mean(wA)
  muA = weighted.mean(mx, wA)
  muB = weighted.mean(mx, wB)
  sdA = sqrt(weighted.mean((mx - muA)^2, wA))
  sdB = sqrt(weighted.mean((mx - muB)^2, wB))

  pA   =    lambda    * dnorm(mx, mean = muA, sd = sdA)
  pB   = (1 - lambda) * dnorm(mx, mean = muB, sd = sdB)
  ptot = pA + pB
  wA   = pA / ptot
  wB   = pB / ptot

  loglikOld = loglik
  loglik = sum(log(pA + pB))
  delta = abs(loglikOld - loglik)
  iter = iter + 1
}
iter

## -----------------------------------------------------------------------------
#| label: EM4
#| echo: !expr c(2)
.o = options(digits = 3)
c(lambda, muA, muB, sdA, sdB)
options(.o)

## -----------------------------------------------------------------------------
#| label: compareEM
#| echo: !expr c(2:3)
.o = options(digits = 3)
gm = mixtools::normalmixEM(mx, k = 2)
with(gm, c(lambda[1], mu, sigma))
options(.o)

## -----------------------------------------------------------------------------
#| label: flexmix
library("flexmix")
data("NPreg")

## -----------------------------------------------------------------------------
#| label: flexmix2
m1 = flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)

## -----------------------------------------------------------------------------
#| label: fig-npreg
#| fig-asp: 1
#| fig-width: 3
#| out-width: 33%
#| fig-cap: "The points seem to come from two different generative processes, one is linear; the other quadratic."
ggplot(NPreg, aes(x = x, y = yn)) + geom_point()

## -----------------------------------------------------------------------------
#| label: parm1c1
modeltools::parameters(m1, component = 1)
modeltools::parameters(m1, component = 2)

## -----------------------------------------------------------------------------
#| label: tableNP
table(NPreg$class, modeltools::clusters(m1))

## -----------------------------------------------------------------------------
#| label: summ1
#| results: false
summary(m1)

## -----------------------------------------------------------------------------
#| label: fig-npregC
#| out-width: 50%
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "Regression example using `flexmix` with the points colored according to their estimated class. You can see that at the intersection we have an `identifiability' problem: we cannot distinguish points that belong to the straight line from ones that belong to the parabole."
NPreg = mutate(NPreg, gr = factor(class))
ggplot(NPreg, aes(x = x, y = yn, group = gr)) +
   geom_point(aes(colour = gr, shape = gr)) +
   scale_colour_hue(l = 40, c = 180)
