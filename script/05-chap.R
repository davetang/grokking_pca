
## -----------------------------------------------------------------------------
#| label: starlings-copyrightfree
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-SnowMapSmallest-web
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: book-chunk-1
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-bombings
#| echo: false
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: fig-cancerHC
#| echo: false
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: fig-ClusteringA
#| echo: false
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "**Of a feather**: how the distances are measured and similarities between observations defined has a strong impact on the clustering result."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-fourdistances
#| echo: false
#| column: margin
#| layout-nrow: 2
#| fig-cap: "Equal-distance contour plots according to four different distances: points on any one curve are all the same distance from the center point."
#| fig-subcap:
#|   - ""
#|   - ""
#|   - ""
    'FourDistances_a.png', 'FourDistances_b.png',
    'FourDistances_c.png', 'FourDistances_d.png')))

## -----------------------------------------------------------------------------
#| label: fig-Mahalanobis
#| echo: false
#| column: margin
#| fig-width: 4.7
#| fig-height: 4.5
#| fig-cap: "An example for the use of Mahalanobis distances to measure the distance of a new data point (red) from two cluster centers."
library("MASS")
library("RColorBrewer")
set.seed(101)
n <- 60000
S1=matrix(c(1,.72,.72,1), ncol=2)
S2=matrix(c(1.5,-0.6,-0.6,1.5),ncol=2)
mu1=c(.5,2.5)
mu2=c(6.5,4)
X1 = mvrnorm(n, mu=c(.5,2.5), Sigma=matrix(c(1,.72,.72,1), ncol=2))
X2 = mvrnorm(n,mu=c(6.5,4), Sigma=matrix(c(1.5,-0.6,-0.6,1.5),ncol=2))
# A color palette from blue to yellow to red
k = 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
plot(X1, xlim=c(-4,12),ylim=c(-2,9), xlab="Orange", ylab="Red", pch='.', cex=1)
points(X2, pch='.', cex=1)
# Draw the colored contour lines
# compute 2D kernel density, see MASS book, pp. 130-131
z1 = kde2d(X1[,1], X1[,2], n=50)
z2 = kde2d(X2[,1], X2[,2], n=50)
contour(z1, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
contour(z2, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
points(3.2,2,pch=20,cex=2.2,col="red")
lines(c(3.2,6.5),c(2,4),col="red",lwd=3)
lines(c(3.2,.5),c(2,2.5),col="red",lwd=3)

## -----------------------------------------------------------------------------
#| label: fig-DistanceTriangle
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: xz
mx  = c(0, 0, 0, 1, 1, 1)
my  = c(1, 0, 1, 1, 0, 1)
mz  = c(1, 1, 1, 0, 1, 1)
mat = rbind(mx, my, mz)
dist(mat)
dist(mat, method = "binary")

## -----------------------------------------------------------------------------
#| label: morder
load("../data/Morder.RData")
sqrt(sum((Morder[1, ] - Morder[2, ])^2))
as.matrix(dist(Morder))[2, 1]

## -----------------------------------------------------------------------------
#| label: HIVmut
mut = read.csv("../data/HIVmutations.csv")
mut[1:3, 10:16]

## -----------------------------------------------------------------------------
#| label: answer-hiv
#| echo: !expr c(2:6)
.o = options(digits = 3)
library("vegan")
mutJ = vegdist(mut, "jaccard")
mutC = sqrt(2 * (1 - cor(t(mut))))
mutJ
as.dist(mutC)
options(.o)

## -----------------------------------------------------------------------------
#| label: fig-xkcd-birds
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "The centers of the groups are sometimes called medoids, thus the name PAM (partitioning around medoids)."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-clust-kmeansastep1
#| warning.known: !expr c("did not converge in 1 iteration")
#| echo: false
#| layout-nrow: 1
#| fig-height: 3
#| fig-cap: "An example run of the $k$-means algorithm. The initial, randomly chosen centers (black circles) and groups (colors) are shown in (a). The group memberships are assigned based on their distance to centers. At each iteration (b) and (c), the group centers are redefined, and the points reassigned to the cluster centers."
#| fig-subcap:
#|   - ""
#|   - ""
#|   - ""
#|   - ""
#
# This code is from Wolfgang. Should it be moved elsewhere, e.g. convert into an exercise?
# The code could be simplied using a loop
#
set.seed(248811)
Xmat = matrix(runif(100), ncol = 2)
nk = 3
cents = Xmat[sample(nrow(Xmat), nk, replace = FALSE), ]
# default distance: Euclidean
dist1 = function(vec){dist(rbind(vec, cents[1,]))}
dist2 = function(vec){dist(rbind(vec, cents[2,]))}
dist3 = function(vec){dist(rbind(vec, cents[3,]))}
dists123 = cbind(apply(Xmat, 1, dist1),
                 apply(Xmat, 1, dist2),
                 apply(Xmat, 1, dist3))
clust0 = apply(dists123, 1, which.min)
out1 = kmeans(Xmat, cents, iter.max=1)
out2 = kmeans(Xmat, cents, iter.max=3)
data0 = data.frame(x = Xmat[,1],
                   y = Xmat[,2],
                   cluster = as.factor(clust0))
data1 = data.frame(x = Xmat[,1],
                   y = Xmat[,2],
                   cluster = as.factor(out1$cluster))
data2 = data.frame(x = Xmat[,1],
                   y = Xmat[,2],
                   cluster = as.factor(out2$cluster))
library("ggplot2")
.mp = function(v, cdg) {
  ggplot(data = v, aes(x = x, y = y)) +
  geom_point(aes(col = cluster, shape = cluster), size = 5) + 
  geom_point(data = cdg, fill = "black", size = 7, shape = 1) +
  scale_shape_discrete(solid = TRUE, guide = "none") + 
  xlab("") + ylab("") + guides(col = "none") + coord_fixed()
}

# centers of clusters:
cdg = data.frame(x = cents[,1],y = cents[,2])
.mp(data0, cdg)

cents = out1$centers
cdg1 = data.frame(x=cents[,1],y=cents[,2])
.mp(data1, cdg1)

cents = out2$centers
cdg2 = data.frame(x=cents[,1],y=cents[,2])
.mp(data2, cdg2)

## -----------------------------------------------------------------------------
#| label: fig-quiltclust-1
#| echo: !expr c(3:9)
#| out-width: 75%
#| fig-height: 4
#| fig-width: 5
#| fig-cap: "Comparison of clustering results (rows), for different numbers of included genes and for varying numbers of clusters, $k$. Each column of the heatmap corresponds to a cell, and the colors represent the cluster assignments."
.oldMar = par("mar")
par(mar = c(1.1, 6, 4.1, 1.1))
library("clusterExperiment")
fluidigm = scRNAseq::ReprocessedFluidigmData()
se = fluidigm[, fluidigm$Coverage_Type == "High"]
assays(se) = list(normalized_counts = 
   round(limma::normalizeQuantiles(assay(se))))
ce = clusterMany(se, clusterFunction = "pam", ks = 5:10, run = TRUE,
  isCount = TRUE, reduceMethod = "var", nFilterDims = c(60, 100, 150))
clusterLabels(ce) = sub("FilterDims", "", clusterLabels(ce))
plotClusters(ce, whichClusters = "workflow", axisLine = -1)
par(mar = .oldMar)

## -----------------------------------------------------------------------------
#| label: quiltclust2-check
#| echo: false
stopifnot(length(clusterLabels(ce)) == 18)

## -----------------------------------------------------------------------------
#| label: book-chunk-2
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: flowCore
library("flowCore")
library("flowViz")
fcsB = read.FCS("../data/Bendall_2011.fcs", truncate_max_range = FALSE)
slotNames(fcsB)

## -----------------------------------------------------------------------------
#| label: RenameCols
markersB = readr::read_csv("../data/Bendall_2011_markers.csv")
mt = match(markersB$isotope, colnames(fcsB))
stopifnot(!any(is.na(mt)))
colnames(fcsB)[mt] = markersB$marker

## -----------------------------------------------------------------------------
#| label: fig-ObviousClusters
#| column: margin
#| fig-cap: "Cell measurements that show clear clustering in two dimensions."
flowPlot(fcsB, plotParameters = colnames(fcsB)[2:3], logy = TRUE)

## -----------------------------------------------------------------------------
#| label: v1v3
#| fig-show: hide
v1 = seq(0, 1, length.out = 100)
plot(log(v1), asinh(v1), type = 'l')
plot(v1, asinh(v1), type = 'l')
v3 = seq(30, 3000, length = 100)
plot(log(v3), asinh(v3), type= 'l')

## -----------------------------------------------------------------------------
#| label: fig-plotTransformations
#| layout-nrow: 1
#| out-width: 75%
#| fig-height: 3
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Panel (a) shows the histogram of the CD3all variable: the cells are clustered around 0 with a few large values. In (b), we see that after an asinh transformation, the cells cluster and fall into two groups or types."
#| fig-subcap:
#|   - ""
#|   - ""
asinhtrsf = arcsinhTransform(a = 0.1, b = 1)
fcsBT = transform(fcsB, transformList(colnames(fcsB)[-c(1, 2, 41)], asinhtrsf))
densityplot(~`CD3all`, fcsB)
densityplot(~`CD3all`, fcsBT)

## -----------------------------------------------------------------------------
#| label: fcsBT
kf = kmeansFilter("CD3all" = c("Pop1","Pop2"), filterId="myKmFilter")
fres = flowCore::filter(fcsBT, kf)
summary(fres)
fcsBT1 = flowCore::split(fcsBT, fres, population = "Pop1")
fcsBT2 = flowCore::split(fcsBT, fres, population = "Pop2")

## -----------------------------------------------------------------------------
#| label: fig-flowCD3CD56-1
#| results: hide
#| column: margin
#| fig-cap: "After transformation these cells were clustered using `kmeans`."
library("flowPeaks")
fp = flowPeaks(Biobase::exprs(fcsBT)[, c("CD3all", "CD56")])
plot(fp)

## -----------------------------------------------------------------------------
#| label: fig-groupcontourCD3CD56
#| column: margin
#| fig-width: 5
#| fig-height: 5
#| fig-cap: "Like @fig-flowCD3CD56-1, using contours."
flowPlot(fcsBT, plotParameters = c("CD3all", "CD56"), logy = FALSE)
contour(fcsBT[, c(40, 19)], add = TRUE)

## -----------------------------------------------------------------------------
#| label: ggcytoCD4CD8
library("ggcyto")
library("labeling")

p1 = ggcyto(fcsB, aes(x = CD4)) + geom_histogram(bins = 60)
p2 = ggcyto(fcsB, aes(x = CD8)) + geom_histogram(bins = 60)
p3 = ggcyto(fcsB, aes(x = CD4, y = CD8)) + geom_density2d(colour = "black")

fcsBT = transform(fcsB, transformList(colnames(fcsB)[-c(1, 2, 41)], 
                                      arcsinhTransform(a = 0, b = 1)))
                                      
p1t = ggcyto(fcsBT, aes(x = CD4))            + geom_histogram(bins = 90)
p2t = ggcyto(fcsBT, aes(x = CD4,y = CD8))    + geom_density2d(colour = "black")
p3t = ggcyto(fcsBT, aes(x = CD45RA,y = CD20))+ geom_density2d(colour = "black")

## -----------------------------------------------------------------------------
#| label: dbscanfcs5
library("dbscan")
mc5 = Biobase::exprs(fcsBT)[, c(15,16,19,40,33)]
res5 = dbscan::dbscan(mc5, eps = 0.65, minPts = 30)
mc5df = data.frame(mc5, cluster = as.factor(res5$cluster))
table(mc5df$cluster)

## -----------------------------------------------------------------------------
#| label: fig-dbscanfcs5
#| warning.known: !expr c("zero contours were generated", "no non-missing arguments to")
#| layout-nrow: 1
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "These two plots show the results of clustering with `dbscan` using five markers. Here we only show the projections of the data into the CD4-CD8 and C3all-CD20 planes."
#| fig-subcap:
#|   - ""
#|   - ""
ggplot(mc5df, aes(x=CD4,    y=CD8,  col=cluster))+geom_density2d()
ggplot(mc5df, aes(x=CD3all, y=CD20, col=cluster))+geom_density2d()

## -----------------------------------------------------------------------------
#| label: mc6
mc6 = Biobase::exprs(fcsBT)[, c(15, 16, 19, 33, 25, 40)]
res = dbscan::dbscan(mc6, eps = 0.65, minPts = 20)
mc6df = data.frame(mc6, cluster = as.factor(res$cluster))
table(mc6df$cluster)

## -----------------------------------------------------------------------------
#| label: mc7
mc7 = Biobase::exprs(fcsBT)[, c(11, 15, 16, 19, 25, 33, 40)]
res = dbscan::dbscan(mc7, eps = 0.95, minPts = 20)
mc7df = data.frame(mc7, cluster = as.factor(res$cluster))
table(mc7df$cluster)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "It is important that the method looks for high density of points in a neighborhood. Other methods exist that try to define clusters by a void, or \"missing points\" between clusters. But these are vulnerable to the curse of dimensionality; these can create spurious \"voids\"."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-LinnaeusClass
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-Cluster-dendrogramorder
#| echo: false
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: fig-single
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-complete
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-betweenwithin
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-mobile
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: threetreeshapes
#| eval: false
#| echo: false
## ####This is for documentation purposes, had to paste
## ####trees together par(mfrow)) not working for pheatmap
## library("pheatmap")
## load("../data/d14.RData")
## pheatmap(d14,clustering_distance_rows=d14,treeheight_col =200,
## cellwidth=20,cellheight=10,lwd=5,treeheight_row=0,clustering_method = "single",
## labels_col=1:11,main="single")
## pheatmap(d14,clustering_distance_rows=d14,treeheight_col =200,cellwidth=20,
## cellheight=10,lwd=5,treeheight_row=0,clustering_method = "average",
## labels_col=1:11,main="average")
## pheatmap(d14,clustering_distance_rows=d14,treeheight_col =200,cellwidth=20,
## cellheight=10,lwd=5,treeheight_row=0,clustering_method = "complete",labels_col=1:11,
## main="complete")

## -----------------------------------------------------------------------------
#| label: fortherecord-2
#| eval: false
#| echo: false
## ####### For the eecord: this is what was done to the data
## ####### Melanoma/Tcell Data: Peter Lee, Susan Holmes, PNAS.
## load("../data/Msig3transp.RData")
## celltypes=factor(substr(rownames(Msig3transp),7,9))
## status=factor(substr(rownames(Msig3transp),1,3))
## Msig2=as.matrix(Msig3transp)
## rownames(Msig2)=substr(rownames(Msig2),1,9)
## hm1=heatmap(as.matrix(dist(Msig2)))
## Morder=Msig2[hm1$rowInd,]
## save(Morder,file="../data/Morder.RData")
## write.table(Morder,"../data/Morder.txt")

## -----------------------------------------------------------------------------
#| label: hclust30Tcells
#| eval: false
#| echo: false
## library("gplots")
## library("pheatmap")
## library("RColorBrewer")
## load("../data/Morder.RData")
## celltypes=factor(substr(rownames(Morder),7,9))
## status=factor(substr(rownames(Morder),1,3))
## ##Just the Euclidean distance
## pheatmap(as.matrix(dist(Morder)),cluster_rows=FALSE,
##         cluster_cols=FALSE,cellwidth=10,cellheight=10)
## ###Manhattan
## pheatmap(as.matrix(dist(Morder,"manhattan")),cluster_rows=FALSE,
##         cluster_cols=FALSE,cellwidth=10,cellheight=10)

## -----------------------------------------------------------------------------
#| label: corT
#| eval: false
#| echo: false
## pheatmap(corT,clustering_distance_rows=distcor,
##        annotation_row=samplesdata[,c("celltypes","status")],
##        show_rownames = FALSE, show_colnames = FALSE)
## pheatmap(corT,clustering_distance_rows=distcor,treeheight_row =150,
##        annotation_row=samplesdata[,c("celltypes","status")],
##        show_rownames = FALSE, show_colnames = FALSE)
## pheatmap(corT,clustering_distance_rows=distcor,treeheight_row =150,
##        annotation_row=samplesdata[,c("celltypes","status")],
##        treeheight_col =150,
##        show_rownames = FALSE, show_colnames = FALSE)
## pheatmap(corT,clustering_distance_rows=distcor,treeheight_row =150,
##           annotation_col=samplesdata[,c("celltypes","status")],
##        annotation_row=samplesdata[,c("celltypes","status")],
##        treeheight_col =150,
##        show_rownames = FALSE, show_colnames = FALSE)

## -----------------------------------------------------------------------------
#| label: fig-treeshapes
#| echo: false
#| layout-nrow: 1
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Three hierarchical clustering plots made with different agglomeration choices. Note the comb-like structure for single linkage in (a). The average (b) and complete linkage (c) trees only differ by the lengths of their inner branches."
#| fig-subcap:
#|   - ""
#|   - ""
    'single14heatmap.png', 'average14heatmap.png', 'complete14heatmap.png')))

## -----------------------------------------------------------------------------
#| label: fig-apeclust14
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: Fake4
library("dplyr")
simdat = lapply(c(0, 8), function(mx) {
  lapply(c(0,8), function(my) {
    tibble(x = rnorm(100, mean = mx, sd = 2),
           y = rnorm(100, mean = my, sd = 2),
           class = paste(mx, my, sep = ":"))
   }) %>% bind_rows
}) %>% bind_rows
simdat
simdatxy = simdat[, c("x", "y")] # without class label

## -----------------------------------------------------------------------------
#| label: fig-simdat-1
#| column: margin
#| fig-width: 3
#| fig-height: 2.3
#| fig-cap: "The `simdat` data colored by the class labels. Here, we know the labels since we generated the data -- usually we do not know them."
ggplot(simdat, aes(x = x, y = y, col = class)) + geom_point() +
  coord_fixed()

## -----------------------------------------------------------------------------
#| label: fig-WSS
#| column: margin
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "The barchart of the WSS statistic as a function of $k$ shows that the last substantial jump is just before $k=4$. This indicates that the best choice for these data is $k=4$."
wss = tibble(k = 1:8, value = NA_real_)
wss$value[1] = sum(scale(simdatxy, scale = FALSE)^2)
for (i in 2:nrow(wss)) {
  km  = kmeans(simdatxy, centers = wss$k[i])
  wss$value[i] = sum(km$withinss)
}
ggplot(wss, aes(x = k, y = value)) + geom_col()

## -----------------------------------------------------------------------------
#| label: fig-CHIndex-1
#| out-width: 33%
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "The Calinski-Harabasz index, i.,e., the ratio of the between and within group variances for different choices of $k$, computed on the `simdat` data."
library("fpc")
library("cluster")
CH = tibble(
  k = 2:8,
  value = sapply(k, function(i) {
    p = pam(simdatxy, i)
    calinhara(simdatxy, p$cluster)
  })
)
ggplot(CH, aes(x = k, y = value)) + geom_line() + geom_point() +
  ylab("CH index")

## -----------------------------------------------------------------------------
#| label: roulette-chunk-3
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-GapStat-1
#| out-width: 33%
#| fig-width: 3
#| fig-height: 4
#| fig-cap: "The gap statistic, see @prp-clustering-gapstat."
library("cluster")
library("ggplot2")
pamfun = function(x, k)
  list(cluster = pam(x, k, cluster.only = TRUE))

gss = clusGap(simdatxy, FUN = pamfun, K.max = 8, B = 50,
              verbose = FALSE)
plot_gap = function(x) {
  gstab = data.frame(x$Tab, k = seq_len(nrow(x$Tab)))
  ggplot(gstab, aes(k, gap)) + geom_line() +
    geom_errorbar(aes(ymax = gap + SE.sim,
                      ymin = gap - SE.sim), width=0.1) +
    geom_point(size = 3, col=  "red")
}
plot_gap(gss)

## -----------------------------------------------------------------------------
#| label: Hiiragi
#| cache: false
library("Hiiragi2013")
data("x")

## -----------------------------------------------------------------------------
#| label: submatH
selFeats = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:50]
embmat = t(Biobase::exprs(x)[selFeats, ])
embgap = clusGap(embmat, FUN = pamfun, K.max = 24, verbose = FALSE)
k1 = maxSE(embgap$Tab[, "gap"], embgap$Tab[, "SE.sim"])
k2 = maxSE(embgap$Tab[, "gap"], embgap$Tab[, "SE.sim"],
           method = "Tibs2001SEmax")
c(k1, k2)

## -----------------------------------------------------------------------------
#| label: checkassertion
#| echo: false
stopifnot("firstSEmax" == eval(formals(maxSE)$method)[1])

## -----------------------------------------------------------------------------
#| label: fig-gapHiiragi-1
#| column: margin
#| fig-width: 4.5
#| fig-height: 4.5
#| fig-cap: "The gap statistic for the **[Hiiragi2013](https://bioconductor.org/packages/Hiiragi2013/)** data."
plot(embgap, main = "")
cl = pamfun(embmat, k = k1)$cluster
table(pData(x)[names(cl), "sampleGroup"], cl)

## -----------------------------------------------------------------------------
#| label: fig-BootstrapClusterNew
#| echo: false
#| layout-nrow: 1
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Different samples from the same distribution $F$ lead to different clusterings. In (a), we see the true sampling variability. The bootstrap simulates this sampling variability by drawing subsamples using the empirical distribution function $\\hat{F}_n$ as shown in (b)."
#| fig-subcap:
#|   - ""
    'BootstrapClusterNew.png', 'BootstrapCluster2New.png')))

## -----------------------------------------------------------------------------
#| label: BootstrapCluster
clusterResampling = function(x, ngenes = 50, k = 2, B = 250,
                             prob = 0.67) {
  mat = Biobase::exprs(x)
  ce = cl_ensemble(list = lapply(seq_len(B), function(b) {
    selSamps = sample(ncol(mat), size = round(prob * ncol(mat)),
                      replace = FALSE)
    submat = mat[, selSamps, drop = FALSE]
    sel = order(rowVars(submat), decreasing = TRUE)[seq_len(ngenes)]
    submat = submat[sel,, drop = FALSE]
    pamres = pam(t(submat), k = k)
    pred = cl_predict(pamres, t(mat[sel, ]), "memberships")
    as.cl_partition(pred)
  }))
  cons = cl_consensus(ce)
  ag = sapply(ce, cl_agreement, y = cons)
  list(agreements = ag, consensus = cons)
}

## -----------------------------------------------------------------------------
#| label: ce1
iswt = (x$genotype == "WT")
cr1 = clusterResampling(x[, x$Embryonic.day == "E3.25" & iswt])
cr2 = clusterResampling(x[, x$Embryonic.day == "E3.5"  & iswt])

## -----------------------------------------------------------------------------
#| label: fig-figClue1
#| echo: !expr -c(12)
#| fig-width: 8
#| fig-height: 4
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Cluster stability analysis with E3.25 and E3.5 samples. Left: beeswarm plots of the cluster agreements with the consensus, for the `B` clusterings; $1$ indicates perfect agreement, lower values indicate lower degrees of agreement. Right: membership probabilities of the consensus clustering. For E3.25, the probabilities are diffuse, indicating that the individual clusterings often disagree, whereas for E3.5, the distribution is bimodal, with only one ambiguous sample."
ag1 = tibble(agreements = cr1$agreements, day = "E3.25")
ag2 = tibble(agreements = cr2$agreements, day = "E3.5")
p1 <- ggplot(bind_rows(ag1, ag2), aes(x = day, y = agreements)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(cex = 1.5, col = "#0000ff40")
mem1 = tibble(y = sort(cl_membership(cr1$consensus)[, 1]),
              x = seq(along = y), day = "E3.25")
mem2 = tibble(y = sort(cl_membership(cr2$consensus)[, 1]),
              x = seq(along = y), day = "E3.5")
p2 <- ggplot(bind_rows(mem1, mem2), aes(x = x, y = y, col = day)) +
  geom_point() + facet_grid(~ day, scales = "free_x")
gridExtra::grid.arrange(p1, p2, widths = c(2.4,4.0))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "**Computational complexity**. An algorithm is said to be $O(n^k)$, if, as $n$ gets larger, the resource consumption (CPU time or memory) grows proportionally to $n^k$. There may be other (sometimes considerable) baseline costs, or costs that grow proportionally to lower powers of $n$, but these always become negligible compared to the leading term as $n\\to\\infty$."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-seqradius
#| column: margin
#| fig-width: 3.1
#| fig-height: 2.6
#| fig-cap: "Although both groups have noise distributions with the same variances, the apparent radii of the groups are very different. The $10^{5}$ instances in `seq2` have many more opportunities for errors than what we see in `seq1`, of which there are only $10^{3}$. Thus we see that frequencies are important in clustering the data."
library("mixtools")
seq1 = rmvnorm(n = 1e3, mu = -c(1, 1), sigma = 0.5 * diag(c(1, 1)))
seq2 = rmvnorm(n = 1e5, mu =  c(1, 1), sigma = 0.5 * diag(c(1, 1)))
twogr = data.frame(
  rbind(seq1, seq2),
  seq = factor(c(rep(1, nrow(seq1)),
                 rep(2, nrow(seq2))))
)
colnames(twogr)[1:2] = c("x", "y")
library("ggplot2")
ggplot(twogr, aes(x = x, y = y, colour = seq,fill = seq)) +
  geom_hex(alpha = 0.5, bins = 50) + coord_fixed()

## -----------------------------------------------------------------------------
#| label: book-chunk-3
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: seqradius2
n    = 2000
len  = 200
perr = 0.001
seqs = matrix(runif(n * len) >= perr, nrow = n, ncol = len)

## -----------------------------------------------------------------------------
#| label: dists
dists = as.matrix(dist(seqs, method = "manhattan"))

## -----------------------------------------------------------------------------
#| label: fig-diameter
#| out-width: 50%
#| fig-width: 4
#| fig-height: 3
#| fig-cap: "The diameter of a set of sequences as a function of the number of sequences."
library("tibble")
dfseqs = tibble(
  k = 10 ^ seq(log10(2), log10(n), length.out = 20),
  diameter = vapply(k, function(i) {
    s = sample(n, i)
    max(dists[s, s])
    }, numeric(1)))
ggplot(dfseqs, aes(x = k, y = diameter)) + geom_point()+geom_smooth()

## -----------------------------------------------------------------------------
#| label: fig-seqradiusex
#| out-width: 50%
#| fig-cap: "`distplot` for the `simseq10K` data."
simseq10K = replicate(1e5, sum(rpois(200, 0.0005)))
mean(simseq10K)
vcd::distplot(simseq10K, "poisson")

## -----------------------------------------------------------------------------
#| label: rates
#| results: false
derepFs = readRDS(file="../data/derepFs.rds")
derepRs = readRDS(file="../data/derepRs.rds")
library("dada2")
ddF = dada(derepFs, err = NULL, selfConsist = TRUE)
ddR = dada(derepRs, err = NULL, selfConsist = TRUE)

## -----------------------------------------------------------------------------
#| label: fig-rerrorprofile1
#| warning.known: !expr c("transformation introduced infinite values")
#| out-width: 75%
#| fig-width: 8
#| fig-height: 6
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Forward transition error rates as provided by `plotErrors(ddF)`. This shows the frequencies of each type of nucleotide transition as a function of quality."
plotErrors(ddF)

## -----------------------------------------------------------------------------
#| label: dada
#| results: false
dadaFs = dada(derepFs, err=ddF[[1]]$err_out, pool = TRUE)
dadaRs = dada(derepRs, err=ddR[[1]]$err_out, pool = TRUE)

## -----------------------------------------------------------------------------
#| label: merge
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs)

## -----------------------------------------------------------------------------
#| label: seqtab
seqtab.all = makeSequenceTable(mergers[!grepl("Mock",names(mergers))])

## -----------------------------------------------------------------------------
#| label: checkdada
#| echo: false
dadada = unique(vapply(dadaRs, class, character(1)))
stopifnot(is.list(dadaRs), identical("dada", dadada))

## -----------------------------------------------------------------------------
#| label: answer-mergers
#| echo: false
length(dadaRs)
length(dadaFs)
class(dadaRs)
names(dadaRs)
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs)
class(mergers)
length(mergers)

## -----------------------------------------------------------------------------
#| label: chimeras
seqtab = removeBimeraDenovo(seqtab.all)

## -----------------------------------------------------------------------------
#| label: Silhouette4
#| fig-show: hide
library("cluster")
pam4 = pam(simdatxy, 4)
sil = silhouette(pam4, 4)
plot(sil, col=c("red","green","blue","purple"), main="Silhouette")

## -----------------------------------------------------------------------------
#| label: ggplotdistheatmap
#| column: margin
#| include: false
#| eval: false
## jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
## paletteSize <- 256
## jBuPuPalette <- jBuPuFun(paletteSize)
## dd=as.matrix(dist.dune)
## prdune <- data.frame(sample = colnames(dd),
##                         probe = rownames(dd),
##                         dist = dd)
## ggplot(prdune, aes(x = probe, y = sample, fill = dist)) +
##   geom_tile() +
##   scale_fill_gradient2(low = jBuPuPalette[1],
##                        mid = jBuPuPalette[paletteSize/2],
##                        high = jBuPuPalette[paletteSize],
##                        midpoint = (max(prdune$dist) + min(prdune$dist)) / 2,
##                        name = "Distance")

## -----------------------------------------------------------------------------
#| label: heatmapDist
#| include: false
#| eval: false
## ## To Do: use pheatmap
## library("graphics")
## library("gplots")
## rc=heat.colors(21, alpha = 1)
## dr=round(as.matrix(dist.dune),1)
## heatmap.2(1-as.matrix(dist.dune),symm = TRUE, margins = c(3,3),Rowv = NA, Colv = NA,col=rc,
## distfun=function(c) as.dist(c), cellnote=dr,key=FALSE)

## -----------------------------------------------------------------------------
#| label: fig-kmeanspital1
#| echo: false
#| layout-nrow: 1
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "An example of non-convex clusters. In (a), we show the result of $k$-means clustering with $k=2$. In (b), we have the output from `dbscan`. The colors represent the three clusters found by the algorithm for the settings \\texttt{eps = 0.16, minPts = 3}."
#| fig-subcap:
#|   - ""
#|   - ""
library("kernlab")
data("spirals")
clusts = kmeans(spirals,2)$cluster
plot(spirals, col = c("blue", "red")[clusts])
data("spirals", package = "kernlab")
res.dbscan = dbscan::dbscan(spirals, eps = 0.16, minPts = 3)
plot(spirals,col=c("blue","red","forestgreen")[res.dbscan$cluster])

## -----------------------------------------------------------------------------
#| label: checkdbscan
#| echo: false
stopifnot(identical(range(res.dbscan$cluster), c(1L, 3L)))

## -----------------------------------------------------------------------------
#| label: specc
#| eval: false
#| echo: false
## sc = specc(spirals, centers=2)
## plot(spirals, col=sc)

## -----------------------------------------------------------------------------
#| label: dada2setup
base_dir = "../data"
miseq_path = file.path(base_dir, "MiSeq_SOP")
filt_path = file.path(miseq_path, "filtered")
fnFs = sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs = sort(list.files(miseq_path, pattern="_R2_001.fastq"))
sampleNames = sapply(strsplit(fnFs, "_"), `[`, 1)
if (!file_test("-d", filt_path)) dir.create(filt_path)
filtFs = file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
fnFs = file.path(miseq_path, fnFs)
fnRs = file.path(miseq_path, fnRs)
print(length(fnFs))

## -----------------------------------------------------------------------------
#| label: fig-profile-1
#| layout-nrow: 1
#| fig-height: 6
#| fid-width: 4
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Quality scores. The lines show positional summary statistics: green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles."
#| fig-subcap:
#|   - ""
#|   - ""
plotQualityProfile(fnFs[1:2]) + ggtitle("Forward")
plotQualityProfile(fnRs[1:2]) + ggtitle("Reverse")

## -----------------------------------------------------------------------------
#| label: profilerev
#| warning.known: !expr c("is deprecated. please use")
#| fig-show: hide
#| results: false
ii = sample(length(fnFs), 4)
plotQualityProfile(fnFs[ii]) + ggtitle("Forward")
plotQualityProfile(fnRs[ii]) + ggtitle("Reverse")

## -----------------------------------------------------------------------------
#| label: filter
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
        maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,  trimLeft=10,
        compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

## -----------------------------------------------------------------------------
#| label: derep
derepFs = derepFastq(filtFs, verbose = FALSE)
derepRs = derepFastq(filtRs, verbose = FALSE)
names(derepFs) = sampleNames
names(derepRs) = sampleNames
