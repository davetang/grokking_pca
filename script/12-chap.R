
## -----------------------------------------------------------------------------
#| label: BuildWall
#| echo: false
#| column: margin
#| out-width: 50%

## -----------------------------------------------------------------------------
#| label: fig-overfitting-1
#| echo: false
#| column: margin
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "An example for **overfitting**: two regression lines are fit to data in the $(x, y)$-plane (black points). We can think of such a line as a rule that predicts the $y$-value, given an $x$-value. Both lines are smooth, but the fits differ in what is called their **bandwidth**, which intuitively can be interpreted their stiffness. The blue line seems overly keen to follow minor wiggles in the data, while the orange line captures the general trend but is less detailed. The effective number of parameters needed to describe the blue line is much higher than for the orange line. Also, if we were to obtain additional data, it is likely that the blue line would do a **worse** job than the orange line in modeling the new data. We'll formalize these concepts --training error and test set error-- later in this chapter. Although exemplified here with line fitting, the concept applies more generally to prediction models."
library("tibble")
library("ggplot2")
ov = tibble(
  x = seq(0, 30, by = 1),
  y = 2 + 0.01 * x^2 + 0.1 * x + 2 * rnorm(length(x)))
ggplot(ov, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(span = 0.2, col = "dodgerblue3", se = FALSE) +
  geom_smooth(span = 0.8, col = "darkorange1", se = FALSE)

## -----------------------------------------------------------------------------
#| label: fig-fourtypes
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: diabetes
data("diabetes", package = "rrcov")
head(diabetes)

## -----------------------------------------------------------------------------
#| label: fig-ldagroups-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 7.5
#| fig-cap: "We see already from the one-dimensional distributions that some of the individual variables could potentially predict which group a patient is more likely to belong to. Our goal is to combine variables to improve over such one-dimensional prediction models."
library("reshape2")
ggplot(melt(diabetes, id.vars = "group"), aes(x = value, col = group)) +
 geom_density() + facet_wrap( ~variable, ncol = 1, scales = "free") +
 theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
#| label: fig-cellshape
#| echo: false
#| out-width: 75%
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: fig-scatterdiabetes-1
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "Scatterplot of two of the variables in the `diabetes` data. Each point is a sample, and the color indicates the diabetes type as encoded in the `group` variable."
ggdb = ggplot(mapping = aes(x = sspg, y = glucose)) +
  geom_point(aes(colour = group), data = diabetes)
ggdb

## -----------------------------------------------------------------------------
#| label: ldaresults
library("MASS")
diabetes_lda = lda(group ~ sspg + glucose, data = diabetes)
diabetes_lda
ghat = predict(diabetes_lda)$class
table(ghat, diabetes$group)
mean(ghat != diabetes$group)

## -----------------------------------------------------------------------------
#| label: make1Dgrid
make1Dgrid = function(x) {
  rg = grDevices::extendrange(x)
  seq(from = rg[1], to = rg[2], length.out = 100)
}

## -----------------------------------------------------------------------------
#| label: diabetes_grid_1
diabetes_grid = with(diabetes,
  expand.grid(sspg = make1Dgrid(sspg),
              glucose = make1Dgrid(glucose)))

## -----------------------------------------------------------------------------
#| label: diabetes_grid_2
diabetes_grid$ghat =
  predict(diabetes_lda, newdata = diabetes_grid)$class

## -----------------------------------------------------------------------------
#| label: centers
centers = diabetes_lda$means

## -----------------------------------------------------------------------------
#| label: unitcircle
unitcircle = exp(1i * seq(0, 2*pi, length.out = 360)) |>
          (\(z) cbind(Re(z), Im(z)))() 
ellipse = unitcircle %*% solve(diabetes_lda$scaling) |> as_tibble()

## -----------------------------------------------------------------------------
#| label: ellipses
#| cache: false
library("dplyr")
ellipses = lapply(rownames(centers), function(gr) {
  mutate(ellipse,
     sspg    = sspg    + centers[gr, "sspg"],
     glucose = glucose + centers[gr, "glucose"],
     group   = gr)
}) |> bind_rows()

## -----------------------------------------------------------------------------
#| label: fig-modeldiabetes-1
#| column: margin
#| fig-height: 4
#| fig-cap: "As @fig-scatterdiabetes-1, with the classification regions from the LDA model shown. The three ellipses represent the class centers and the covariance matrix of the LDA model; note that there is only one covariance matrix, which is the same for all three classes. Therefore also the sizes and orientations of the ellipses are the same for the three classes, only their centers differ. They represent contours of equal class membership probability."
ggdb + geom_raster(aes(fill = ghat),
            data = diabetes_grid, alpha = 0.25, interpolate = TRUE) +
    geom_point(data = as_tibble(centers), pch = "+", size = 8) +
    geom_path(aes(colour = group), data = ellipses) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

## -----------------------------------------------------------------------------
#| label: fig-diabetes-lda-uniform-prior-1
#| fig-width: 5
#| fig-height: 4
#| out-width: 50%
#| fig-cap: "As @fig-modeldiabetes-1, but with uniform class priors."
diabetes_up = lda(group ~ sspg + glucose, data = diabetes,
  prior = (\(n) rep(1/n, n)) (nlevels(diabetes$group)))

diabetes_grid$ghat_up =
  predict(diabetes_up, newdata = diabetes_grid)$class

stopifnot(all.equal(diabetes_up$means, diabetes_lda$means))

ellipse_up  = unitcircle %*% solve(diabetes_up$scaling) |> as_tibble()
ellipses_up = lapply(rownames(centers), function(gr) {
  mutate(ellipse_up,
     sspg    = sspg    + centers[gr, "sspg"],
     glucose = glucose + centers[gr, "glucose"],
     group   = gr)
}) |> bind_rows()

ggdb + geom_raster(aes(fill = ghat_up),
            data = diabetes_grid, alpha = 0.4, interpolate = TRUE) +
    geom_point(data = data.frame(centers), pch = "+", size = 8) +
    geom_path(aes(colour = group), data = ellipses_up) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

## -----------------------------------------------------------------------------
#| label: all5diab
diabetes_lda5 = lda(group ~ rw + fpg + glucose + sspg + insulin, data = diabetes)
diabetes_lda5
ghat5 = predict(diabetes_lda5)$class
table(ghat5, diabetes$group)
mean(ghat5 != diabetes$group)

## -----------------------------------------------------------------------------
#| label: loadHiiragi2
#| cache: false
library("Hiiragi2013")
data("x")
probes = c("1426642_at", "1418765_at", "1418864_at", "1416564_at")
embryoCells = t(Biobase::exprs(x)[probes, ]) |> as_tibble() |>
  mutate(Embryonic.day = x$Embryonic.day) |>
  dplyr::filter(x$genotype == "WT")

## -----------------------------------------------------------------------------
#| label: annoHiiragi
annotation(x)
library("mouse4302.db")
anno = AnnotationDbi::select(mouse4302.db, keys = probes,
                             columns = c("SYMBOL", "GENENAME"))
anno
mt = match(anno$PROBEID, colnames(embryoCells))
colnames(embryoCells)[mt] = anno$SYMBOL

## -----------------------------------------------------------------------------
#| label: assertprobeid
#| echo: false
stopifnot(!any(is.na(mt)))

## -----------------------------------------------------------------------------
#| label: fig-HiiragiFourGenesPairs-1
#| column: margin
#| fig-width: 6
#| fig-height: 6
#| fig-margin: false
#| fig-cap: "Expression values of the discriminating genes, with the prediction target Embryonic.day shown by color."
library("GGally")
ggpairs(embryoCells, mapping = aes(col = Embryonic.day),
  columns = anno$SYMBOL, upper = list(continuous = "points"))

## -----------------------------------------------------------------------------
#| label: ldacells
ec_lda = lda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)
round(ec_lda$scaling, 1)

## -----------------------------------------------------------------------------
#| label: fig-edcontour-1
#| column: margin
#| fig-width: 4.5
#| fig-height: 3.5
#| fig-cap: "LDA classification regions for Embryonic.day."
ec_rot = predict(ec_lda)$x |> as_tibble() |>
           mutate(ed = embryoCells$Embryonic.day)
ec_lda2 = lda(ec_rot[, 1:2], predict(ec_lda)$class)
ec_grid = with(ec_rot, expand.grid(
  LD1 = make1Dgrid(LD1),
  LD2 = make1Dgrid(LD2)))
ec_grid$edhat = predict(ec_lda2, newdata = ec_grid)$class
ggplot() +
  geom_point(aes(x = LD1, y = LD2, colour = ed), data = ec_rot) +
  geom_raster(aes(x = LD1, y = LD2, fill = edhat),
            data = ec_grid, alpha = 0.4, interpolate = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed()

## -----------------------------------------------------------------------------
#| label: fig-qdamouse-1
#| fig-width: 9
#| fig-height: 9
#| fig-cap: "QDA for the mouse cell data. Shown are all pairwise plots of the four features. In each plot, the other two features are set to the median."
library("gridExtra")

ec_qda = qda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)

variables = colnames(ec_qda$means)
pairs = combn(variables, 2)
lapply(seq_len(ncol(pairs)), function(i) {
  grid = with(embryoCells,
    expand.grid(x = make1Dgrid(get(pairs[1, i])),
                y = make1Dgrid(get(pairs[2, i])))) |>
    `colnames<-`(pairs[, i])

  for (v in setdiff(variables, pairs[, i]))
    grid[[v]] = median(embryoCells[[v]])

  grid$edhat = predict(ec_qda, newdata = grid)$class

  x <- pairs[1,i]
  y <- pairs[2,i]
  ggplot() + 
    geom_point(
      data = embryoCells,
      aes(x = .data[[x]], y = .data[[y]], colour = Embryonic.day)
    ) +
    geom_raster(
      aes(x = .data[[x]], y = .data[[y]], fill = edhat),
      data = grid, alpha = 0.4, interpolate = TRUE
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    if (i != ncol(pairs)) theme(legend.position = "none")
}) |> (\(g) grid.arrange(grobs = g, ncol = 2))()

## -----------------------------------------------------------------------------
#| label: ladallvariables
#| warning.known: !expr c("variables are collinear")
#| results: hide
#| error: true
lda(t(Biobase::exprs(x))[, 1:1000], x$Embryonic.day)
warnings()
qda(t(Biobase::exprs(x))[, 1:1000], x$Embryonic.day)

## -----------------------------------------------------------------------------
#| label: fig-learnbyheart-1
#| warning.known: !expr c("variables are collinear")
#| column: margin
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "Misclassification rate of LDA applied to random data. While the number of observations `n` is held constant (at 20), we are increasing the number of features `p` starting from 2 up to 21. The misclassification rate becomes almost zero as `p` approaches 20. The LDA model becomes so elaborate and over-parameterized that it manages to learn the random labels \"by heart\". (As `p` becomes even larger, the \"performance\" degrades again somewhat, apparently due to numerical properties of the `lda` implementation used here.)"
p = 2:21
n = 20

mcl = lapply(p, function(pp) {
  replicate(100, {
    xmat = matrix(rnorm(n * pp), nrow = n)
    resp = sample(c("apple", "orange"), n, replace = TRUE)
    fit  = lda(xmat, resp)
    pred = predict(fit)$class
    mean(pred != resp)
  }) |> mean() |> (\(x) tibble(mcl = x, p = pp))()
}) |> bind_rows()

ggplot(mcl, aes(x = p, y = mcl)) + 
  geom_line() + geom_point() +
  ylab("Misclassification rate")

## -----------------------------------------------------------------------------
#| label: book-chunk-1
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-mclcv-1
#| warning.known: !expr c("variables are collinear")
#| column: margin
#| fig-cap: "Cross-validation: the misclassification rate of LDA applied to random data, when evaluated on test data that were not used for learning, hovers around 0.5 independent of `p`. The misclassification rate on the training data is also shown. It behaves similar to what we already saw in @fig-learnbyheart-1."
estimate_mcl_loocv = function(x, resp) {
  vapply(seq_len(nrow(x)), function(i) {
    fit  = lda(x[-i, ], resp[-i])
    ptrn = predict(fit, newdata = x[-i,, drop = FALSE])$class
    ptst = predict(fit, newdata = x[ i,, drop = FALSE])$class
    c(train = mean(ptrn != resp[-i]), test = (ptst != resp[i]))
  }, FUN.VALUE = numeric(2)) |> rowMeans() |> t() |> as_tibble()
}

xmat = matrix(rnorm(n * last(p)), nrow = n)
resp = sample(c("apple", "orange"), n, replace = TRUE)

mcl = lapply(p, function(k) {
  estimate_mcl_loocv(xmat[, 1:k], resp)
}) |> bind_rows() |> data.frame(p) |> melt(id.var = "p")

ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
  geom_point() + ylab("Misclassification rate")

## -----------------------------------------------------------------------------
#| label: fig-curseofdim
#| warning.known: !expr c("variables are collinear")
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "As we increase the number of features included in the model, the misclassification rate initially improves; as we start including more and more irrelevant features, it increases again, as we are fitting noise."
p   = 2:20
mcl = replicate(100, {
  xmat = matrix(rnorm(n * last(p)), nrow = n)
  resp = sample(c("apple", "orange"), n, replace = TRUE)
  xmat[, 1:6] = xmat[, 1:6] + as.integer(factor(resp))

  lapply(p, function(k) {
    estimate_mcl_loocv(xmat[, 1:k], resp)
  }) |> bind_rows() |> cbind(p = p) |> melt(id.var = "p")
}, simplify = FALSE) |> bind_rows()

mcl = group_by(mcl, p, variable) |> summarise(value = mean(value))

ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
   geom_point() + ylab("Misclassification rate")

## -----------------------------------------------------------------------------
#| label: fig-BiasVarianceTradeoff
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-cursedimans1-1
#| fig-width: 3.5
#| fig-height: 2.5
#| out-width: 50%
#| fig-cap: "Side length of a $p$-dimensional hybercube expected to contain 10 points out of 1 million uniformly distributed ones, as a function of the $p$. While for $p=1$, this length is conveniently small, namely $10/10^6=10^{-5}$, for larger $p$ it approaches 1, i.,e., becomes the same as the range of each the features. This means that a \"local neighborhood\" of 10 points encompasses almost the same data range as the whole dataset."
sideLength = function(p, pointDensity = 1e6, pointsNeeded = 10)
  (pointsNeeded / pointDensity) ^ (1 / p)
ggplot(tibble(p = 1:400, sideLength = sideLength(p)),
       aes(x = p, y = sideLength)) + geom_line(col = "red") +
  geom_hline(aes(yintercept = 1), linetype = 2)

## -----------------------------------------------------------------------------
#| label: fig-cursedimans2-1
#| fig-width: 3.5
#| fig-height: 2.5
#| out-width: 50%
#| fig-cap: "Fraction of a unit cube's total volume that is in its \"shell\" (here operationalised as those points that are closer than 0.01 to its surface) as a function of the dimension $p$."
tibble(
  p = 1:400,
  volOuterCube = 1 ^ p,
  volInnerCube = 0.98 ^ p,  # 0.98 = 1 - 2 * 0.01
  `V(shell)` = volOuterCube - volInnerCube) |>
ggplot(aes(x = p, y =`V(shell)`)) + geom_line(col = "blue")

## -----------------------------------------------------------------------------
#| label: fig-cursedimans3-1
#| fig-width: 3.5
#| fig-height: 2.5
#| out-width: 50%
#| fig-cap: "Coefficient of variation (CV) of the distance between randomly picked points in the unit hypercube, as a function of the dimension. As the dimension increases, everybody is equally far away from everyone else: there is almost no variation in the distances any more."
n = 1000
df = tibble(
  p = round(10 ^ seq(0, 4, by = 0.25)),
  cv = vapply(p, function(k) {
    x1 = matrix(runif(k * n), nrow = n)
    x2 = matrix(runif(k * n), nrow = n)
    d = sqrt(rowSums((x1 - x2)^2))
    sd(d) / mean(d)
  }, FUN.VALUE = numeric(1)))
ggplot(df, aes(x = log10(p), y = cv)) + geom_line(col = "orange") +
  geom_point()

## -----------------------------------------------------------------------------
#| label: confusiontable
#| eval: false
## table(truth, response)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Note that the $\\log\\text{loss}$ will be infinite if a prediction is totally confident ($\\hat{p}_i$ is exactly $0$ or $1$) but wrong."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-supervised-bullseye
#| echo: false
#| column: margin
#| layout-nrow: 1
#| fig-subcap:
#|   - ""
#|   - ""
    IMGS, c('TargetBias.png','TargetVariance.png')))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Here, $|\\beta|^\\nu=\\sum_i\\beta_i^\\nu$ is the $L_\\nu$-norm of the vector $\\beta$. Variations are possible, for instead we could include in this summation only some but not all of the elements of $\\beta$; or we could scale different elements differently, for instance based on some prior belief of their scale and importance."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: colon1
#| cache: false
#| results: hide
library("ExperimentHub")
eh = ExperimentHub()
zeller = eh[["EH361"]]

## -----------------------------------------------------------------------------
#| label: colon1b
table(zeller$disease)

## -----------------------------------------------------------------------------
#| label: colon2
zellerNC = zeller[, zeller$disease %in% c("n", "cancer")]

## -----------------------------------------------------------------------------
#| label: ehzellertest
#| echo: false
stopifnot(is.numeric(Biobase::exprs(zellerNC)), !any(is.na(Biobase::exprs(zellerNC))))

## -----------------------------------------------------------------------------
#| label: zellerpData-1
pData(zellerNC)[ sample(ncol(zellerNC), 3), ]

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "We define the helper function `formatfn` to line wrap these long character strings for the available space here."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: zellerrownames
formatfn = function(x)
   gsub("|", "| ", x, fixed = TRUE) |> lapply(strwrap)

rownames(zellerNC)[1:4]
rownames(zellerNC)[nrow(zellerNC) + (-2:0)] |> formatfn()

## -----------------------------------------------------------------------------
#| label: fig-zellerHist-1
#| column: margin
#| fig-width: 3
#| fig-height: 4
#| fig-cap: "Histograms of the distributions for two randomly selected features. The distributions are highly skewed, with many zero values and a thin, long tail of non-zero values."
ggplot(melt(Biobase::exprs(zellerNC)[c(510, 527), ]), aes(x = value)) +
    geom_histogram(bins = 25) +
    facet_wrap( ~ Var1, ncol = 1, scales = "free")

## -----------------------------------------------------------------------------
#| label: glmnet
library("glmnet")
glmfit = glmnet(x = t(Biobase::exprs(zellerNC)),
                y = factor(zellerNC$disease),
                family = "binomial")

## -----------------------------------------------------------------------------
#| label: colonPred
predTrsf = predict(glmfit, newx = t(Biobase::exprs(zellerNC)),
                   type = "class", s = 0.04)
table(predTrsf, zellerNC$disease)

## -----------------------------------------------------------------------------
#| label: fig-plotglmfit-1
#| echo: !expr -c(1)
#| column: margin
#| fig-width: 3.6
#| fig-height: 3.2
#| fig-cap: "Regularization paths for `glmfit`."
par(mai = c(0.5, 0.5, 0.575, 0.05))
plot(glmfit, col = brewer.pal(8, "Dark2"), lwd = sqrt(3), ylab = "")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "You'll already realize from the description of this strategy that if we optimize $\\lambda$ in this way, the resulting apparent classification performance will likely be exaggerated. We need a truly independent dataset, or at least another, outer cross-validation loop to get a more realistic impression of the generalizability. We will get back to this question at the end of the chapter."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-colonCV-1
#| column: margin
#| fig-width: 4
#| fig-height: 4
#| fig-cap: "Diagnostic plot for `cv.glmnet`: shown is a measure of cross-validated prediction performance, the deviance, as a function of $\\lambda$. The dashed vertical lines show `lambda.min` and `lambda.1se`."
cvglmfit = cv.glmnet(x = t(Biobase::exprs(zellerNC)),
                     y = factor(zellerNC$disease),
                     family = "binomial")
plot(cvglmfit)

## -----------------------------------------------------------------------------
#| label: lambda.min
cvglmfit$lambda.min

## -----------------------------------------------------------------------------
#| label: lambda.1se
cvglmfit$lambda.1se

## -----------------------------------------------------------------------------
#| label: predictwithlambda1se
s0 = cvglmfit$lambda.1se
predict(glmfit, newx = t(Biobase::exprs(zellerNC)),type = "class", s = s0) |>
    table(zellerNC$disease)

## -----------------------------------------------------------------------------
#| label: zellercoef
coefs = coef(glmfit)[, which.min(abs(glmfit$lambda - s0))]
topthree = order(abs(coefs), decreasing = TRUE)[1:3]
as.vector(coefs[topthree])
formatfn(names(coefs)[topthree])

## -----------------------------------------------------------------------------
#| label: fig-colonCVTrsf-1
#| fig-width: 4
#| fig-height: 4
#| out-width: 50%
#| fig-cap: "like @fig-colonCV-1, but using an $\\text{asinh}$ transformation of the data."
cv.glmnet(x = t(asinh(Biobase::exprs(zellerNC))),
          y = factor(zellerNC$disease),
          family = "binomial") |> plot()

## -----------------------------------------------------------------------------
#| label: fig-mousecvglmfit-1
#| column: margin
#| fig-width: 4
#| fig-height: 4
#| fig-cap: "Cross-validated misclassification error versus penalty parameter for the mouse cells data."
sx = x[, x$Embryonic.day == "E3.25"]
embryoCellsClassifier = cv.glmnet(t(Biobase::exprs(sx)), sx$genotype,
                family = "binomial", type.measure = "class")
plot(embryoCellsClassifier)

## -----------------------------------------------------------------------------
#| label: checkclaimMouseCellsClassifier
#| echo: false
stopifnot(sum((diff(embryoCellsClassifier$cvm) * diff(embryoCellsClassifier$lambda)) < 0) <= 2)

## -----------------------------------------------------------------------------
#| label: fig-mousecellsrowttst-1
#| fig-width: 4
#| fig-height: 2.5
#| out-width: 50%
#| fig-cap: "Histogram of p-values for the per-feature $t$-tests between genotypes in the E3.25 cells."
mouse_de = rowttests(sx, "genotype")
ggplot(mouse_de, aes(x = p.value)) +
  geom_histogram(boundary = 0, breaks = seq(0, 1, by = 0.01))

## -----------------------------------------------------------------------------
#| label: mousecellsnn1
dists = as.matrix(dist(scale(t(Biobase::exprs(x)))))
diag(dists) = +Inf

## -----------------------------------------------------------------------------
#| label: mousecellsnn2
nn = sapply(seq_len(ncol(dists)), function(i) which.min(dists[, i]))
table(x$sampleGroup, x$sampleGroup[nn]) |> `colnames<-`(NULL)

## -----------------------------------------------------------------------------
#| label: caret1
library("caret")
caretMethods = names(getModelInfo())
head(caretMethods, 8)
length(caretMethods)

## -----------------------------------------------------------------------------
#| label: caret2
getModelInfo("nnet", regex = FALSE)[[1]]$parameter

## -----------------------------------------------------------------------------
#| label: caret3
#| results: hide
trnCtrl = trainControl(
  method = "repeatedcv",
  repeats = 3,
  classProbs = TRUE)
tuneGrid = expand.grid(
  size = c(2, 4, 8),
  decay = c(0, 1e-2, 1e-1))
nnfit = train(
  Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
  data = embryoCells,
  method = "nnet",
  tuneGrid  = tuneGrid,
  trControl = trnCtrl,
  metric = "Accuracy")

## -----------------------------------------------------------------------------
#| label: fig-ML-nnfit
#| column: margin
#| fig-width: 3.75
#| fig-height: 4.25
#| fig-cap: "Parameter tuning of the neural net by cross-validation."
nnfit
plot(nnfit)
predict(nnfit) |> head(10)

## -----------------------------------------------------------------------------
#| label: kernelsvm
#| eval: false
#| echo: false
## library("kernlab")
## kfunction= function(linear =0, quadratic=0)
## {  k = function (v,w){ linear*sum((v)*(w)) + quadratic*sum((v^2)*(w^2))}
##   class(k) = "kernel"
##   return(k) }
## subx=subx[,2:3]
## svp = ksvm(subx,dftxy$tg,type="C-svc",C = 100, kernel=kfunction(1,0),scaled=c())
## plot(c(min(subx[,1]), max(subx[,1])),c(min(subx[,2]), max(subx[,2])),
##             type='n',xlab='x1',ylab='x2')
## ymat = ymatrix(svp)
## points(subx[-SVindex(svp),1], subx[-SVindex(svp),2],
##          pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
## points(subx[SVindex(svp),1], subx[SVindex(svp),2],
##          pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
## 
## # Extract w and b from the model
## w = colSums(coef(svp)[[1]] * subx[SVindex(svp),])
## b = b(svp)
## # Draw the lines
## abline(b/w[2],-w[1]/w[2])
## abline((b+1)/w[2],-w[1]/w[2],lty=2)
## abline((b-1)/w[2],-w[1]/w[2],lty=2)
