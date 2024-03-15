
## -----------------------------------------------------------------------------
#| label: Brighton-West-Pier-20090214-sunset
#| echo: false

## -----------------------------------------------------------------------------
#| label: UkrCities
library("pheatmap")
data("ukraine_dists", package = "MSMB")
as.matrix(ukraine_dists)[1:4, 1:4]

## -----------------------------------------------------------------------------
#| label: fig-HeatDists
#| column: margin
#| fig-height: 6
#| fig-width: 6.9
#| fig-cap: "A heatmap of the `ukraine_dists` distance matrix. Distances are measured in metres. The function has re-arranged the order of the cities, and grouped the closest ones."
pheatmap(as.matrix(ukraine_dists), 
  color = colorRampPalette(c("#0057b7", "#ffd700"))(50),
  breaks = seq(0, max(ukraine_dists)^(1/2), length.out = 51)^2,
  treeheight_row = 10, treeheight_col = 10)

## -----------------------------------------------------------------------------
#| label: cmdscaleDistE
ukraine_mds = cmdscale(ukraine_dists, eig = TRUE)

## -----------------------------------------------------------------------------
#| label: defplotscree
#| cache: false
library("dplyr")
library("ggplot2")
plotscree = function(x, m = length(x$eig)) {
  ggplot(tibble(eig = x$eig[seq_len(m)], k = seq(along = eig)),
    aes(x = k, y = eig)) + theme_minimal() +
    scale_x_discrete("k", limits = as.factor(seq_len(m))) + 
    geom_bar(stat = "identity", width = 0.5, fill = "#ffd700", col = "#0057b7")
}

## -----------------------------------------------------------------------------
#| label: fig-plotscreeeig
#| column: margin
#| fig-height: 4
#| fig-width: 4
#| fig-cap: "Screeplot of the first four eigenvalues. There is a pronounced drop after the first two eigenvalues, which indicates that the data are well described by a two-dimensional embedding."
plotscree(ukraine_mds, m = 4)

## -----------------------------------------------------------------------------
#| label: fig-plotscreeeigall
#| fig-width: 6
#| fig-height: 4
#| out-width: 50%
#| fig-cap: "Screeplot of all the eigenvalues."
ukraine_mds$eig |> signif(3)
plotscree(ukraine_mds)

## -----------------------------------------------------------------------------
#| label: checkassertionbarploteigenvaluescmdscale
#| echo: false
stopifnot(any(ukraine_mds$eig < 0))

## -----------------------------------------------------------------------------
#| label: fig-ukrainemds1
#| column: margin
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "MDS map based on the distances."
#| echo: !expr -(2:4)
#| eval: true
ukraine_mds_df = tibble(
  PCo1 = ukraine_mds$points[, 1],
  PCo2 = ukraine_mds$points[, 2],
  labs = rownames(ukraine_mds$points)
)
if (with(ukraine_mds_df, labs[which.max(PCo1)] != "Luhansk"))
  ukraine_mds_df$PCo1 = -ukraine_mds_df$PCo1
if (with(ukraine_mds_df, labs[which.max(PCo2)] != "Sevastopol"))
  ukraine_mds_df$PCo2 = -ukraine_mds_df$PCo2
if(with(ukraine_mds_df,
     labs[which.max(PCo1)] != "Luhansk" ||
     labs[which.max(PCo2)] != "Sevastopol"))
  stop("There is an error with 'ukraine_mds_df'.")
library("ggrepel")
g = ggplot(ukraine_mds_df, aes(x = PCo1, y = PCo2, label = labs)) +
  geom_point() + geom_text_repel(col = "#0057b7") + coord_fixed() 
g

## -----------------------------------------------------------------------------
#| label: fig-ukrainemds2
#| column: margin
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "Same as @fig-ukrainemds1, but with y-axis flipped."
g %+% mutate(ukraine_mds_df, PCo1 = PCo1, PCo2 = -PCo2)

## -----------------------------------------------------------------------------
#| label: fig-ukrainecoord1
#| column: margin
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "True latitudes and longitudes, taken from the `ukraine_coords` dataframe."
data("ukraine_coords", package = "MSMB")
print.data.frame(ukraine_coords[1:4,  c("city", "lat", "lon")])
ggplot(ukraine_coords, aes(x = lon, y = lat, label = city)) +
  geom_point() + geom_text_repel(col = "#0057b7")

## -----------------------------------------------------------------------------
#| label: EarthRadius
#| echo: false
earthradius = 6371

## -----------------------------------------------------------------------------
#| label: fig-ukrainecoord3
#| eval: false
## library("sf")
## library("rnaturalearth")
## library("rnaturalearthdata")
## world = ne_countries(scale = "medium", returnclass = "sf")
## ggplot() +
##   geom_sf(data = world,
##           fill = ifelse(world$geounit == "Ukraine", "#ffd700", "#f0f0f0")) +
##   coord_sf(xlim = range(ukraine_coords$lon) + c(-0.5, + 2),
##            ylim = range(ukraine_coords$lat) + c(-0.5, + 1)) +
##   geom_point(aes(x = lon, y = lat), data = ukraine_coords) +
##   geom_text_repel(aes(x = lon, y = lat, label = city),
##                   color = "#0057b7", data = ukraine_coords)

## -----------------------------------------------------------------------------
#| label: fig-ukrainecoord2
#| echo: false
#| fig-width: 6
#| fig-height: 4
#| out-width: 75%
#| fig-cap: "International borders added to @fig-ukrainecoord1."
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
world = ne_countries(scale = "medium", returnclass = "sf")
ggplot() + 
  geom_sf(data = world, 
          fill = ifelse(world$geounit == "Ukraine", "#ffd700", "#f0f0f0")) +
  coord_sf(xlim = range(ukraine_coords$lon) + c(-0.5, + 2), 
           ylim = range(ukraine_coords$lat) + c(-0.5, + 1)) +
  geom_point(aes(x = lon, y = lat), data = ukraine_coords) + 
  geom_text_repel(aes(x = lon, y = lat, label = city), 
                  color = "#0057b7", data = ukraine_coords)

## -----------------------------------------------------------------------------
#| label: cos48
#| echo: false
stopifnot(48 == round(mean(range(ukraine_coords$lat))))

## -----------------------------------------------------------------------------
#| label: D2
X = with(ukraine_coords, cbind(lon, lat * cos(48)))
DdotD = as.matrix(dist(X)^2)

## -----------------------------------------------------------------------------
#| label: CheckMDS
n = nrow(X)
H = diag(rep(1,n))-(1/n) * matrix(1, nrow = n, ncol = n)
Xc = sweep(X,2,apply(X,2,mean))
Xc[1:2, ]
HX = H %*% X
HX[1:2, ]
apply(HX, 2, mean)

## -----------------------------------------------------------------------------
B0 = H  %*% DdotD %*% H
B2 = HX %*% t(HX)
B2[1:3, 1:3] / B0[1:3, 1:3]
max(abs(-0.5 * B0 - B2))

## -----------------------------------------------------------------------------
#| label: ekmandis
ekm = read.table("../data/ekman.txt", header=TRUE)
rownames(ekm) = colnames(ekm)
disekm = 1 - ekm - diag(1, ncol(ekm))
disekm[1:5, 1:5]
disekm = as.dist(disekm)

## -----------------------------------------------------------------------------
#| label: fig-ekmanMDSeig
#| column: margin
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "The screeplot shows us that the phenomenon is largely two dimensional."
mdsekm = cmdscale(disekm, eig = TRUE)
plotscree(mdsekm)

## -----------------------------------------------------------------------------
#| label: fig-ekmanMDS-1
#| column: margin
#| fig-width: 4.6
#| fig-height: 4.5
#| fig-cap: "The layout of the scatterpoints in the first two dimensions has a horseshoe shape. The labels and colors show that the arch corresponds to the wavelengths."
dfekm = mdsekm$points[, 1:2] |>
  `colnames<-`(paste0("MDS", 1:2)) |>
  as_tibble() |>
  mutate(
    name = rownames(ekm),
    rgb = photobiology::w_length2rgb(as.numeric(sub("w", "", name))))
ggplot(dfekm, aes(x = MDS1, y = MDS2)) +
  geom_point(col = dfekm$rgb, size = 4) +
  geom_text_repel(aes(label = name)) + coord_fixed()

## -----------------------------------------------------------------------------
#| label: simstress
#| output: false
#| results: hide
library("vegan")
nmds.stress = function(x, sim = 100, kmax = 4) {
  sapply(seq_len(kmax), function(k)
    replicate(sim, metaMDS(x, k = k, autotransform = FALSE)$stress))
}
stress = nmds.stress(disekm, sim = 100)
dim(stress)

## -----------------------------------------------------------------------------
#| label: fig-NMDSscreeplot-1
#| column: margin
#| fig-height: 4
#| fig-width: 4
#| fig-cap: "Several replicates at each dimension were run to evaluate the stability of the <span style=\"font-variant:small-caps;\">stress</span>. We see that the <span style=\"font-variant:small-caps;\">stress</span> drops dramatically with two or more dimensions, thus indicating that a two dimensional solution is appropriate here."
dfstr = reshape2::melt(stress, varnames = c("replicate","dimensions"))
ggplot(dfstr, aes(y = value, x = dimensions, group = dimensions)) +
  geom_boxplot()

## -----------------------------------------------------------------------------
#| label: fig-Shepardsplot-1
#| results: hide
#| column: margin
#| fig-cap: "The Shepard's plot compares the original distances or dissimilarities (along the horizonal axis) to the reconstructed distances, in this case for $k=2$ (vertical axis)."
nmdsk2 = metaMDS(disekm, k = 2, autotransform = FALSE)
stressplot(nmdsk2, pch = 20)

## -----------------------------------------------------------------------------
#| label: fig-ekmannonMDS-1
#| echo: !expr -c(1:2)
#| layout-nrow: 1
#| fig-width: 4.6
#| fig-height: 4.5
#| fig-margin: false
#| fig-cap: "Comparison of the output from (a) the classical multidimensional scaling (same as @fig-ekmanMDS-1) and (b) the nonmetric version."
#| fig-subcap:
#|   - ""
#|   - ""
ggplot(dfekm, aes(x = MDS1, y = MDS2)) +
  geom_point(col = dfekm$rgb, size = 4) +
  geom_text_repel(aes(label = name)) + coord_fixed()
  
nmdsk2$points[, 1:2] |> 
  `colnames<-`(paste0("NmMDS", 1:2)) |>
  as_tibble() |> 
  bind_cols(dplyr::select(dfekm, rgb, name)) |>
  ggplot(aes(x = NmMDS1, y = NmMDS2)) +
    geom_point(col = dfekm$rgb, size = 4) +
    geom_text_repel(aes(label = name))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "**Metadata:** Many programs and workflows for biological sequence analysis or assays separate the environmental and contextual information, which they call **metadata**, from the assay data or sequence reads. We discourage such practice as the exact connections between the samples and covariates are important. A lost connection between the assays and covariates makes later analyses impossible. Covariates such as clinical history, time, batch or location are important and should be considered components of the data."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "**Bioconductor container:** These data are an example of an awkward way of combining batch information with the actual data. The `day` information has been combined with the array data and encoded as a number and could be confused with a continuous variable. We will see in the next section a better practice for storing and manipulating heterogeneous data using a Bioconductor container called *SummarizedExperiment*."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: vsn28Exprs
IBDchip = readRDS("../data/vsn28Exprd.rds")
library("ade4")
library("factoextra")
library("sva")

## -----------------------------------------------------------------------------
#| label: IBDchip
class(IBDchip)
dim(IBDchip)
tail(IBDchip[,1:3])
table(IBDchip[nrow(IBDchip), ])

## -----------------------------------------------------------------------------
#| label: assayIBD
assayIBD = IBDchip[-nrow(IBDchip), ]
day      = factor(IBDchip[nrow(IBDchip), ])

## -----------------------------------------------------------------------------
#| label: fig-screepc12
#| fig-width: 4
#| fig-height: 4
#| fig-show: asis
#| fig-cap: "The screeplot shows us that the samples can be usefully represented in a two dimensional embedding."
rankthreshPCA = function(x, threshold = 3000) {
  ranksM = apply(x, 2, rank)
  ranksM[ranksM < threshold] = threshold
  ranksM = threshold - ranksM
  dudi.pca(t(ranksM), scannf = FALSE, nf = 2)
}
pcaDay12 = rankthreshPCA(assayIBD[, day != 3])
fviz_eig(pcaDay12, bar_width = 0.6) + ggtitle("")

## -----------------------------------------------------------------------------
#| label: fig-rankthreshPCA
#| fig-width: 5.2
#| fig-height: 3.5
#| fig-show: asis
#| fig-cap: "We have used colors to identify the different days and have kept the sample labels as well. We have also added convex hulls for each day. The group mean is identified as the point with the larger symbol (circle, triangle or square)."
day12 = day[ day!=3 ]
rtPCA1 = fviz(pcaDay12, element = "ind", axes = c(1, 2), geom = c("point", "text"),
  habillage = day12, repel = TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + ggtitle("") +
  coord_fixed()
rtPCA1

## -----------------------------------------------------------------------------
#| label: Threesetspca123-code
#| warning.known: !expr c("unlabeled data points")
#| fig-show: hide
pcaDay123 = rankthreshPCA(assayIBD)
fviz(pcaDay123, element = "ind", axes = c(1, 2), geom = c("point", "text"),
  habillage = day, repel = TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + 
  ggtitle("") + coord_fixed()

## -----------------------------------------------------------------------------
#| label: fig-Threesetspca123
#| warning.known: !expr c("unlabeled data points")
#| echo: false
#| layout-nrow: 1
#| fig-width: 4.5
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "When comparing the three day analysis to that of the first two days, we notice the inversion of signs in the coordinates on the second axis: this has no biological relevance. The important finding is that group 3 overlaps heavily with group 1 indicating that it was the protocol change on Day 2 which created the variability."
#| fig-subcap:
#|   - ""
#|   - ""
rtPCA1
fviz(pcaDay123, element="ind", axes=c(1,2), geom=c("point","text"),
  habillage = day, repel=TRUE, palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "convex") + ggtitle("") +
  coord_fixed()

## -----------------------------------------------------------------------------
#| label: exCodeEllipse
#| fig-show: hide
fviz_pca_ind(pcaDay123, habillage = day, labelsize = 3,
  palette = "Dark2", addEllipses = TRUE, ellipse.level = 0.69)

## -----------------------------------------------------------------------------
#| label: fig-screepc123
#| echo: false
#| column: margin
#| fig-width: 4
#| fig-height: 4
#| fig-cap: "The eigenvalue screeplot the case of 3 groups is extremely similar to that with two groups shown in @fig-screepc12."
fviz_eig(pcaDay123, bar_width=0.6) + ggtitle("")

## -----------------------------------------------------------------------------
#| label: fig-CombatIBD
#| column: margin
#| fig-width: 5.2
#| fig-height: 3.5
#| fig-cap: "The modified data with the batch effects removed now show three batch-groups heavily overlapping and centered almost at the origin."
model0 = model.matrix(~1, day)
combatIBD = ComBat(dat = assayIBD, batch = day, mod = model0)
pcaDayBatRM = rankthreshPCA(combatIBD)
fviz(pcaDayBatRM, element = "ind", geom = c("point", "text"),
  habillage = day, repel=TRUE, palette = "Dark2", addEllipses = TRUE,
  ellipse.type = "convex", axes =c(1,2)) + coord_fixed() + ggtitle("")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "A confusing notational similarity occurs here, in the SummarizedExperiment framework a `DataFrame` is not the same as a *data.frame.*"

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: SummarizedExperiment
#| cache: false
library("SummarizedExperiment")
treatment  = factor(ifelse(grepl("Cntr|^C", colnames(IBDchip)), "CTL", "IBS"))
sampledata = DataFrame(day = day, treatment = treatment)
chipse = SummarizedExperiment(assays  = list(abundance = assayIBD),
                              colData = sampledata)

## -----------------------------------------------------------------------------
#| label: checktreatment
#| echo: false
# check whether the 'grep'ed treatment status agree with
# a manual 'hardcoded' version provided by Susan previously.
stopifnot(identical(chipse$treatment, factor(c("IBS", "CTL")[
  c(1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1)])))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "You can explore composite objects using the Environment pane in RStudio. You will see that in `chipse`, some of the slots are empty."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: SEfilter
chipse[, day == 2]

## -----------------------------------------------------------------------------
#| label: SingleCell
#| cache: false
corese = readRDS("../data/normse.rds")
norm = assays(corese)$normalizedValues

## -----------------------------------------------------------------------------
#| label: coldatacore
length(unique(colData(corese)$Batch))

## -----------------------------------------------------------------------------
#| label: fig-screeplotnorm
#| column: margin
#| fig-width: 4
#| fig-height: 4
#| fig-cap: "Screeplot of the PCA of the normalized data."
respca = dudi.pca(t(norm), nf = 3, scannf = FALSE)
plotscree(respca, 15)
PCS = respca$li[, 1:3]

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "We have set up colors for the clusters as in the workflow, (the code is not shown here)."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: setupcolors
#| echo: false
library("RColorBrewer")
publishedClusters = colData(corese)[, "publishedClusters"]
batch = colData(corese)$Batch
col_clus = c(
  "transparent", "#1B9E77", "antiquewhite2", "cyan", 
  "#E7298A", "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", 
  "darkorchid2", "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")
names(col_clus) = sort(unique(publishedClusters))

## -----------------------------------------------------------------------------
#| label: screeplotnorm-2
library("rgl")
batch = colData(corese)$Batch
plot3d(PCS,aspect=sqrt(c(84,24,20)),col=col_clus[batch])
plot3d(PCS,aspect=sqrt(c(84,24,20)),
col = col_clus[as.character(publishedClusters)])

## -----------------------------------------------------------------------------
#| label: fig-3dplotsnorm
#| echo: false
#| layout-nrow: 1
#| fig-width: 4
#| fig-height: 2
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Two-dimensional screenshots of three-dimensional **[rgl](https://cran.r-project.org/web/packages/rgl/)** plots. The points are colored according to batch numbers in (a), and according to the original clustering in (b). We can see that the batch effect has been effectively removed and that the cells show the original clustering."
#| fig-subcap:
#|   - ""
    'plotnormpcabatch1.png', 'plotnormpclust1.png')))

## -----------------------------------------------------------------------------
#| label: tbl-HIV
#| echo: false
#| column: margin
#| tbl-cap-location: top
#| tbl-cap: "Sample by mutation matrix."
HIV <- data.frame(
  Patient = c('AHX112', 'AHX717', 'AHX543'), 
  Mut1 = c(0, 1, 1),
  Mut2 = c(0, 0, 0),
  Mut3 = c(0, 1, 0),
  '...' = rep(' ', 3))
knitr::kable(HIV, format = 'html')

## -----------------------------------------------------------------------------
#| label: tbl-crossHIV
#| echo: false
#| column: margin
#| tbl-cap-location: top
#| tbl-cap: "Cross-tabulation of the HIV mutations showing two-way co-occurrences."
crossHIV <- data.frame(
  Patient = c('Mut1', 'Mut2', 'Mut3'), 
  Mut1 = c(853, 29, 10),
  Mut2 = c(29, 853, 52),
  Mut3 = c(10, 52, 853),
  '...' = rep(' ', 3))
knitr::kable(crossHIV, format = 'html')

## -----------------------------------------------------------------------------
#| label: fig-HIVnnrti
#| column: margin
#| fig-cap: "The dependencies between HIV mutations is clearly a three dimensional phenomenon, the three first eigenvalues show a clear signal in the data."
cooc = read.delim2("../data/coccurHIV.txt", header = TRUE, sep = ",")
cooc[1:4, 1:11]
HIVca = dudi.coa(cooc, nf = 4, scannf = FALSE)
fviz_eig(HIVca, geom = "bar", bar_width = 0.6) + ggtitle("")

## -----------------------------------------------------------------------------
#| label: fig-HIV3d
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: CA123-HIV
library("rgl")
CA1=HIVca$li[,1];CA2=HIVca$li[,2];CA3=HIVca$li[,3]
plot3d(CA1,CA2,CA3,aspect=FALSE,col="purple")

## -----------------------------------------------------------------------------
#| label: fig-HIVca
#| layout-nrow: 1
#| fig-height: 3
#| fig-width: 4.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Two planar maps of the mutations defined with the horizontal axis corresponding to the first eigenvector of the CA and the vertical axis being the second axis in (a), and the third in (b); notice the difference in heights."
#| fig-subcap:
#|   - ""
#|   - ""
fviz_ca_row(HIVca,axes = c(1, 2),geom="text", col.row="purple",
  labelsize=3)+ggtitle("") + xlim(-0.55, 1.7) + ylim(-0.53,1.1) +
  theme_bw() +  coord_fixed()
fviz_ca_row(HIVca,axes = c(1, 3), geom="text",col.row="purple",
    labelsize=3)+ggtitle("")+ xlim(-0.55, 1.7)+ylim(-0.5,0.6) +
    theme_bw() + coord_fixed()

## -----------------------------------------------------------------------------
#| label: HIVnnrtiMut13a
#| fig-show: hide
fviz_ca_row(HIVca, axes=c(1, 3), geom="text", col.row="purple", labelsize=3) +
  ggtitle("") + theme_minimal() + coord_fixed()

## -----------------------------------------------------------------------------
#| label: chisquaredtest
#| warning.known: !expr c("approximation may be incorrect")
HairColor = HairEyeColor[,,2]
chisq.test(HairColor)

## -----------------------------------------------------------------------------
#| label: tbl-HairEye
#| echo: false
#| column: margin
#| tbl-cap-location: top
#| tbl-cap: "Cross tabulation of students hair and eye color."
knitr::kable(HairColor, format = 'html')

## -----------------------------------------------------------------------------
#| label: ExpectedEyes
rowsums = as.matrix(apply(HairColor, 1, sum))
rowsums
colsums = as.matrix(apply(HairColor, 2, sum))
t(colsums)
HCexp = rowsums %*%t (colsums) / sum(colsums)

## -----------------------------------------------------------------------------
#| label: fig-HCexp
#| echo: false
#| column: margin
#| fig-show: hide
#| fig-width: 4
#| fig-height: 3.8
#| fig-cap: "Here is a schematic representation of the expected table `HCexp`. We see that it has the 'rectangular' property chracteristic of rank one matrices we saw in @sec-multivariate. The boxes are all white."
mosaicplot(HCexp, shade=TRUE, las=1, type="pearson", cex.axis=0.7, main="")

## -----------------------------------------------------------------------------
sum((HairColor  - HCexp)^2/HCexp)

## -----------------------------------------------------------------------------
#| label: fig-MosaicHair
#| out-width: 75%
#| fig-width: 5
#| fig-height: 5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Visualization of the departure from independence. Now, the boxes are proportional in size to the actual observed counts and we no longer have a 'rectangular' property. The departure from independence is measured in Chisquared distance for each of the boxes and colored according to whether the residuals are large and positive. Dark blue indicates a positive association, for instance between blue eyes and blonde hair, red indicates a negative association such as in the case of blond hair and brown eyes."
round(t(HairColor-HCexp))
library("vcd")
mosaicplot(HairColor, shade=TRUE, las=1, type="pearson", cex.axis=0.7, main="")

## -----------------------------------------------------------------------------
#| label: fig-HCscatter
#| out-width: 75%
#| fig-width: 4
#| fig-height: 3
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "The CA plot gives a representation of a large proportion of the chisquare distance between the data and the values expected under independence. The first axis shows a contrast between black haired and blonde haired students, mirrored by the brown eye, blue eye opposition. In CA the two categories play symmetric roles and we can interpret the proximity of Blue eyes and Blond hair has meaning that there is strong co-occurence of these categories."
HC = as.data.frame.matrix(HairColor)
coaHC = dudi.coa(HC,scannf=FALSE,nf=2)
round(coaHC$eig[1:3]/sum(coaHC$eig)*100)
fviz_ca_biplot(coaHC, repel=TRUE, col.col="brown", col.row="purple") +
  ggtitle("") + ylim(c(-0.5,0.5))

## -----------------------------------------------------------------------------
#| label: VeganCCA
#| fig-width: 4
#| fig-height: 3
#| fig-show: hide
library("vegan")
res.ca = vegan::cca(HairColor)
plot(res.ca, scaling=3)

## -----------------------------------------------------------------------------
#| label: ProustProxy
#| echo: false

## -----------------------------------------------------------------------------
#| label: LakeCAr
load("../data/lakes.RData")
lakelike[1:3,1:8]
reslake=dudi.coa(lakelike,scannf=FALSE,nf=2)
round(reslake$eig[1:8]/sum(reslake$eig),2)

## -----------------------------------------------------------------------------
#| label: fig-LakeCAr
#| layout-nrow: 1
#| fig-width: 4
#| fig-height: 2.5
#| fig-margin: false
#| fig-cap: "The locations near the lake are ordered along an arch as shown in (a). In the biplot (b), we can see which plants are most frequent at which locations by looking at the red triangles closest to the blue points."
#| fig-subcap:
#|   - ""
#|   - ""
fviz_ca_row(reslake,repel=TRUE)+ggtitle("")+ylim(c(-0.55,1.7))
fviz_ca_biplot(reslake,repel=TRUE)+ggtitle("")+ylim(c(-0.55,1.7))

## -----------------------------------------------------------------------------
#| label: fig-CellTree
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: Descriptors
#| eval: false
#| echo: false
#| results: hide
## # Provenance tracking, keep this for the record:
## Nor = read.csv("../data/nbt.3154-S3.csv",row.names=1)
## dim(Nor)
## blom = as.matrix(Nor)
## desc1=unlist(strsplit(rownames(blom),"_"))
## desc=desc1[seq(1,7867,2)]
## gr4sfg=which(substr(rownames(blom),1,5)=="4SFGA")
## gr4sf=which(substr(rownames(blom),1,4)=="4SGA")
## gr1=which(substr(rownames(blom),1,2)=="PS")
## gr2=which(substr(rownames(blom),1,2)=="NP")
## gr3=which(substr(rownames(blom),1,2)=="HF")
## colscells=c("blue","green","orange","red","purple")
## colnb=rep(0,3934)
## colnb[gr1]=1
## colnb[gr2]=2
## colnb[gr3]=3
## colnb[gr4sf]=4
## colnb[gr4sfg]=5
## typesort=rep(0,3934)
## typesort[ nchar(desc) < 5 & substr(rownames(blom), 3, 3) == "A"] = "sortA"
## typesort[ nchar(desc) < 5 & substr(rownames(blom), 3, 3) == "B"] = "sortB"
## typesort[ nchar(desc) >= 5 ] = "sortA"
## ftable(typesort)
## celltypes=as.factor(c("PS","NP","HF","4SG","4SGF-")[colnb])
## cellcol = colscells[colnb]
## colCells = DataFrame(celltypes=celltypes, cellcol=colscells[colnb])
## Moignard= SummarizedExperiment(assays=list(assayCells = blom), rowData=colCells)
## saveRDS(Moignard,file="../data/Moignard.rds")

## -----------------------------------------------------------------------------
#| label: Distances
#| cache: false
Moignard = readRDS("../data/Moignard.rds")
cellt = rowData(Moignard)$celltypes
colsn = c("red", "purple", "orange", "green", "blue")
blom = assay(Moignard)
dist2n.euclid = dist(blom)
dist1n.l1     = dist(blom, "manhattan")

## -----------------------------------------------------------------------------
#| label: CMDS
ce1Mds = cmdscale(dist1n.l1,     k = 20, eig = TRUE)
ce2Mds = cmdscale(dist2n.euclid, k = 20, eig = TRUE)
perc1  = round(100*sum(ce1Mds$eig[1:2])/sum(ce1Mds$eig))
perc2  = round(100*sum(ce2Mds$eig[1:2])/sum(ce2Mds$eig))

## -----------------------------------------------------------------------------
#| label: fig-CMDSplotscree
#| column: margin
#| layout-ncol: 1
#| fig-height: 3
#| fig-width: 2.5
#| fig-cap: "Screeplots from MDS on $\\ell_1$ (a) and $L_2$ (b) distances. We see that the eigenvalues are extremely similar and both point to a $2$ dimensional phenomenon."
#| fig-subcap:
#|   - ""
#|   - ""
plotscree(ce1Mds, m = 4)
plotscree(ce2Mds, m = 4)

## -----------------------------------------------------------------------------
#| label: fig-CMDSplotL2
#| layout-nrow: 1
#| fig-height: 5
#| fig-width: 5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Moignard cell data colored according to the cell types (blue: PS, green: NP, yellow: HF, red: 4SG, purple: 4SFG$^-$) in the two dimensional MDS plots created. In (a) using $\\ell_1$ distances and in (b) using the L2 distances."
#| fig-subcap:
#|   - ""
#|   - ""
c1mds = ce1Mds$points[, 1:2] |>
        `colnames<-`(paste0("L1_PCo", 1:2)) |>
        as_tibble()
ggplot(c1mds, aes(x = L1_PCo1, y = L1_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values = colsn) + guides(color = "none")
c2mds = ce2Mds$points[, 1:2] |>
        `colnames<-`(paste0("L2_PCo", 1:2)) |>
        as_tibble()
ggplot(c2mds, aes(x = L2_PCo1, y = L2_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
   scale_colour_manual(values = colsn) + guides(color = "none")

## -----------------------------------------------------------------------------
#| label: colorlegend
#| echo: false
#| column: margin
#| fig-width: 2
#| fig-height: 3
## Hack from https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
ggpcolor = ggplot(c1mds,aes(x=L1_PCo1,y=L1_PCo2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
  scale_colour_manual(values=colsn, name = "cell type")
g_legend = function(a) {
     gt = ggplot_gtable(ggplot_build(a))
     leg = which(sapply(gt$grobs, function(x) x$name) == "guide-box")
     gt$grobs[[leg]]
}
grid.draw(g_legend(ggpcolor))

## -----------------------------------------------------------------------------
#| label: fig-tsnecells
#| dev: jpeg
#| layout-nrow: 1
#| fig-width: 5
#| fig-height: 4
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "The four cell populations studied here are representative of three sequential states (PS,NP,HF) and two possible final branches (4SG and 4SFG$^{-}$). The plot on the left was obtained by choosing 2 dimensions for t-sne at a perplexity of 30. The lower plot has obtained by choosing 3 dimensions, we can see that this third t-SNE axis represented here as the horizontal axis."
#| fig-subcap:
#|   - ""
#|   - ""
library("Rtsne")
restsne = Rtsne(blom, dims = 2, perplexity = 30, verbose = FALSE,
                max_iter = 900)
dftsne = restsne$Y[, 1:2] |>
         `colnames<-`(paste0("axis", 1:2)) |>
         as_tibble()
ggplot(dftsne,aes(x = axis1, y = axis2, color = cellt)) +
  geom_point(aes(color = cellt), alpha = 0.6) +
   scale_color_manual(values = colsn) + guides(color = "none")
restsne3 = Rtsne(blom, dims = 3, perplexity = 30, verbose = FALSE,
                 max_iter = 900)
dftsne3 = restsne3$Y[, 1:3] |>
          `colnames<-`(paste0("axis", 1:3)) |> 
          as_tibble()
ggplot(dftsne3,aes(x = axis3, y = axis2, group = cellt)) +
      geom_point(aes(color = cellt), alpha = 0.6) +
      scale_colour_manual(values = colsn) + guides(color = "none")

## -----------------------------------------------------------------------------
#| label: fig-tsne3d
#| echo: false
#| layout-nrow: 1
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Moignard cell data colored according to the cell types (blue: PS, green: NP, yellow: HF, red: 4SG, purple: 4SFG$^-$) in the three-dimensional t-SNE layouts. We can see that the purple cells (4SFG$^-$) segregate at the outer shell on the top of the point cloud."
#| fig-subcap:
#|   - ""
    c('tsnemoignard3scrop.png', 'tsnemoignard3crop.png')))

## -----------------------------------------------------------------------------
#| label: fig-ukrainetsne
#| fig-width: 6
#| fig-height: 4
#| out-width: 70%
#| fig-cap: "t-SNE map based of Ukraine."
ukraine_tsne = Rtsne(ukraine_dists, is_distance = TRUE, perplexity = 8)
ukraine_tsne_df = tibble(
  PCo1 = ukraine_tsne$Y[, 1],
  PCo2 = ukraine_tsne$Y[, 2],
  labs = attr(ukraine_dists, "Labels")
)
ggplot(ukraine_tsne_df, aes(x = PCo1, y = PCo2, label = labs)) +
  geom_point() + geom_text_repel(col = "#0057b7") + coord_fixed() 

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "There are some precautions to be taken when using the Mantel coefficient, see a critical review in @Guillot:2013."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "We will see many examples of regularization and danger of overfitting in @sec-supervised."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: multitable-setup
library("genefilter")
load("../data/microbe.rda")
metab = read.csv("../data/metabolites.csv", row.names = 1) |> as.matrix()

## -----------------------------------------------------------------------------
#| label: multitable-filtering
library("phyloseq")
metab   = metab[rowSums(metab == 0) <= 3, ]
microbe = prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe = filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab   = log(1 + metab, base = 10)
X       = log(1 + as.matrix(otu_table(microbe)), base = 10)

## -----------------------------------------------------------------------------
#| label: RVtest
colnames(metab) = colnames(X)
pca1 = dudi.pca(t(metab), scal = TRUE, scann = FALSE)
pca2 = dudi.pca(t(X), scal = TRUE, scann = FALSE)
rv1 = RV.rtest(pca1$tab, pca2$tab, 999)
rv1

## -----------------------------------------------------------------------------
#| label: multitable-sparse-cca
library("PMA")
ccaRes = CCA(t(X), t(metab), penaltyx = 0.15, penaltyz = 0.15, 
             typex = "standard", typez = "standard")
ccaRes

## -----------------------------------------------------------------------------
#| label: multitablepluginpca
#| echo: false
combined = cbind(t(X[ccaRes$u != 0, ]),
                 t(metab[ccaRes$v != 0, ]))
pcaRes = dudi.pca(combined, scannf = FALSE, nf = 3)
# annotation
genotype    = substr(rownames(pcaRes$li), 1, 2)
sampleType  = substr(rownames(pcaRes$l1), 3, 4)
featureType = grepl("\\.", colnames(combined))
featureType = ifelse(featureType, "Metabolite", "OTU")
sampleInfo  = data.frame(pcaRes$li, genotype, diet=sampleType)
featureInfo = data.frame(pcaRes$c1, feature = substr(colnames(combined), 1, 6))

## -----------------------------------------------------------------------------
#| label: fig-multitableinterpretpca
#| warning.known: !expr c("unlabeled data points")
#| echo: false
#| out-width: 75%
#| fig-width: 7
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "A PCA triplot produced from the CCA selected features from muliple data types (metabolites and OTUs)."
ggplot() +  geom_point(data = sampleInfo,
  aes(x = Axis1, y = Axis2, col = diet, shape = genotype), size = 3) +
  geom_label_repel(data = featureInfo,
  aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = featureType),
      size = 2, segment.size = 0.3,
      label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = featureInfo,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = featureType),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed()+
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pcaRes$eig[1] / sum(pcaRes$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pcaRes$eig[2] / sum(pcaRes$eig), 2)),
       fill = "Feature Type", col = "Sample Type")

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "**Notational overload for CCA**: Originally invented by  @terBraak:1985 and called Canonical Correspondence analysis, we will call this method Constrained Correspondence Analysis and abbreviate it CCpnA to avoid confusion with Canonical Correlation Analysis (CCA). However several R packages, such as **[ade4](https://cran.r-project.org/web/packages/ade4/)** and **[vegan](https://cran.r-project.org/web/packages/vegan/)** use the name `cca` for their correspondence analyses function."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: ccpna-correspondence-analysis
ps1=readRDS("../data/ps1.rds")
ps1p=filter_taxa(ps1, function(x) sum(x) > 0, TRUE)
psCCpnA = ordinate(ps1p, "CCA",
                 formula = ps1p ~ ageBin + family_relationship)

## -----------------------------------------------------------------------------
#| label: ccpna-join-data
#| echo: false
library("dplyr")
tax = data.frame(tax_table(ps1p),stringsAsFactors = FALSE)
tax$seq = rownames(tax)
mainOrders = c("Clostridiales", "Bacteroidales",
               "Lactobacillales", "Coriobacteriales")
tax$Order[!(tax$Order %in% mainOrders)] = "Other"
tax$Order = factor(tax$Order, levels = c(mainOrders, "Other"))
tax$otu_id = seq_len(ncol(otu_table(ps1p)))
scoresCCpnA = vegan::scores(psCCpnA)
sites = data.frame(scoresCCpnA$sites)
sites$SampleID = rownames(sites)
sites = left_join(sites, as(sample_data(ps1p), "data.frame"))
species = data.frame(scoresCCpnA$species)
species$otu_id = seq_along(colnames(otu_table(ps1p)))
species = left_join(species, tax)

## -----------------------------------------------------------------------------
#| label: fig-ccpnaplotage
#| warning.known: !expr c("unlabeled data points")
#| column: margin
#| fig-width: 4
#| fig-height: 6
#| fig-margin: false
#| fig-cap: "The mouse and taxa scores generated by CCpnA. The sites (mice samples) are triangles; species are circles, respectively. The separate panels indicate different age groups."
evalProp = 100 * psCCpnA$CCA$eig[1:2] / sum(psCCpnA$CA$eig)
ggplot() +
 geom_point(data = sites,aes(x =CCA2, y =CCA1),shape =2,alpha=0.5) +
 geom_point(data = species,aes(x =CCA2,y =CCA1,col = Order),size=1)+
 geom_text_repel(data = dplyr::filter(species, CCA2 < (-2)),
                   aes(x = CCA2, y = CCA1, label = otu_id),
                   size = 2, segment.size = 0.1) +
 facet_grid(. ~ ageBin) +
 guides(col = guide_legend(override.aes = list(size = 2))) +
 labs(x = sprintf("Axis2 [%s%% variance]", round(evalProp[2])),
      y = sprintf("Axis1 [%s%% variance]", round(evalProp[1]))) +
 scale_color_brewer(palette = "Set1") + theme(legend.position="bottom")

## -----------------------------------------------------------------------------
#| label: fig-ccpnaplotlitter
#| warning.known: !expr c("unlabeled data points")
#| echo: false
#| out-width: 50%
#| fig-width: 4
#| fig-height: 6
#| fig-margin: false
#| fig-cap: "The analogue to @fig-ccpnaplotage, faceting by litter membership rather than age bin."
ggplot() +
  geom_point(data = sites,   aes(x = CCA2, y = CCA1), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA2, y = CCA1, col = Order), size = 1) +
  geom_text_repel(data =  dplyr::filter(species, CCA2 < (-2)),
                  aes(x = CCA2, y = CCA1, label = otu_id),
                  size = 2, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  labs(x = sprintf("Axis2 [%s%% variance]", round(evalProp[2])),
       y = sprintf("Axis1 [%s%% variance]", round(evalProp[1]))) +
  scale_color_brewer(palette = "Set1") + theme(legend.position="bottom")

## -----------------------------------------------------------------------------
#| label: Threesetscode1
ibd.pres = ifelse(assayIBD[, 1:28] > 8.633, 1, 0)

## -----------------------------------------------------------------------------
#| label: fig-Threesetscoa
#| layout-nrow: 1
#| fig-height: 4
#| fig-cap: "Correspondence analysis on binary data."
#| fig-subcap:
#|   - "$\\text{}$"
#|   - "$\\text{}$"
IBDca = dudi.coa(ibd.pres, scannf = FALSE, nf = 4)
fviz_eig(IBDca, geom = "bar", bar_width = 0.7) +
    ylab("Percentage of chisquare") + ggtitle("")
fviz(IBDca, element = "col", axes = c(1, 2), geom = "point",
     habillage = day, palette = "Dark2", addEllipses = TRUE, color = day,
     ellipse.type = "convex", alpha = 1, col.row.sup =  "blue",
     select = list(name = NULL, cos2 = NULL, contrib = NULL),
     repel = TRUE)

## -----------------------------------------------------------------------------
#| label: tbl-colors
#| echo: false
#| tbl-cap-location: top
#| tbl-cap: "Contingency table of co-occurring terms from search engine results."
d1 <- t(data.frame(
  quiet = c(2770, 2150, 2140, 875, 1220, 821, 2510),
  angry = c(2970, 1530, 1740, 752, 1040, 710, 1730),
  clever = c(1650, 1270, 1320, 495, 693, 416, 1420),
  depressed = c(1480, 957, 983, 147, 330, 102, 1270),
  happy = c(19300, 8310, 8730, 1920, 4220, 2610, 9150),
  lively = c(1840, 1250, 1350, 659, 621, 488, 1480),
  perplexed = c(110,  71,  80,  19,  23,  15, 109),
  virtuous = c(179,  80, 102,  20,  25,  17, 165)))
colnames(d1) <- c('black','blue','green','grey','orange','purple','white')
knitr::kable(d1, format='html')

## -----------------------------------------------------------------------------
#| label: fig-ColorBiplot-1
#| echo: false
#| out-width: 50%
#| fig-height: 5
#| fig-width: 4.5
#| fig-cap: "Correspondence Analysis allows for a symmetrical graphical representation of two categorical variables, in this case colors and emotions for a contingency table of co-occurrences such as @tbl-colors."
colorsentiment = read.csv("../data/colorsentiment.csv")
colsent = xtabs(colorsentiment[,3] ~ colorsentiment[,2] + colorsentiment[,1])
coldf = data.frame(unclass(colsent))
coldf = round(coldf / 1000)
# xtable::xtable(round(coldf),display=rep("d", 8))
colorfactor = names(coldf)
veganout = vegan::cca(coldf)
colorfactor[c(4,7)] = c("darkgrey", "grey")
ordiplot(veganout, scaling = 3, type = "none", xlim =c(-1.2, 0.75), ylim =c(-0.7, 1))
text(veganout, "sites", pch = 21, col = "red", bg = "yellow", scaling = 3)
text(veganout, "species", pch = 21, col = colorfactor, bg = "black", cex=1.2, scaling = 3)

## -----------------------------------------------------------------------------
#| label: PlatoTableImage
#| echo: false

## -----------------------------------------------------------------------------
#| label: platoca
platof = read.table("../data/platof.txt", header = TRUE)
platof[1:4, ]

## -----------------------------------------------------------------------------
#| label: fig-platoca
#| column: margin
#| layout-ncol: 1
#| fig-cap: "Biplot of Plato's sentence endings."
#| fig-subcap:
#|   - ""
#|   - ""
resPlato = dudi.coa(platof, scannf = FALSE, nf = 2)
fviz_ca_biplot(resPlato, axes=c(2, 1)) + ggtitle("")
fviz_eig(resPlato, geom = "bar", width = 0.6) + ggtitle("")

## -----------------------------------------------------------------------------
#| label: PercentageInertia
names(resPlato)
sum(resPlato$eig)
percentageInertia=round(100*cumsum(resPlato$eig)/sum(resPlato$eig))
percentageInertia
percentageInertia[2]

## -----------------------------------------------------------------------------
#| label: MakeMatrices
load("../data/lakes.RData")
lakelike[ 1:3, 1:8]
lakelikeh[1:3, 1:8]
e_coa  = dudi.coa(lakelike,  scannf = FALSE, nf = 2)
e_pca  = dudi.pca(lakelike,  scannf = FALSE, nf = 2)
eh_coa = dudi.coa(lakelikeh, scannf = FALSE, nf = 2)
eh_pca = dudi.pca(lakelikeh, scannf = FALSE, nf = 2)

## -----------------------------------------------------------------------------
#| label: adescatter-show
#| output: false
scatter(e_pca)
scatter(e_coa)
s.label(e_pca$li)
s.label(e_coa$li)

s.label(eh_pca$co)
s.label(eh_pca$li)
s.label(eh_coa$li)
s.label(eh_coa$co)

## -----------------------------------------------------------------------------
#| label: RawData
moignard_raw = as.matrix(read.csv("../data/nbt.3154-S3-raw.csv", row.names = 1))
dist2r.euclid = dist(moignard_raw)
dist1r.l1     = dist(moignard_raw, "manhattan")
cells1.cmds = cmdscale(dist1r.l1,     k = 20, eig = TRUE)
cells2.cmds = cmdscale(dist2r.euclid, k = 20, eig = TRUE)
sum(cells1.cmds$eig[1:2]) / sum(cells1.cmds$eig)
sum(cells2.cmds$eig[1:2]) / sum(cells2.cmds$eig)

## -----------------------------------------------------------------------------
#| label: KernelD1
#| results: hide
library("kernlab")
laplacedot1 = laplacedot(sigma = 1/3934)
rbfdot1     = rbfdot(sigma = (1/3934)^2 )
Klaplace_cellsn   = kernelMatrix(laplacedot1, blom)
KGauss_cellsn     = kernelMatrix(rbfdot1, blom)
Klaplace_rawcells = kernelMatrix(laplacedot1, moignard_raw)
KGauss_rawcells   = kernelMatrix(rbfdot1, moignard_raw)

## -----------------------------------------------------------------------------
#| label: KernelD2
#| results: hide
dist1kr = 1 - Klaplace_rawcells
dist2kr = 1 - KGauss_rawcells
dist1kn = 1 - Klaplace_cellsn
dist2kn = 1 - KGauss_cellsn

cells1.kcmds = cmdscale(dist1kr, k = 20, eig = TRUE) 
cells2.kcmds = cmdscale(dist2kr, k = 20, eig = TRUE) 

percentage = function(x, n = 4) round(100 * sum(x[seq_len(n)]) / sum(x[x>0]))
kperc1 = percentage(cells1.kcmds$eig)
kperc2 = percentage(cells2.kcmds$eig)

cellsn1.kcmds = cmdscale(dist1kn, k = 20, eig = TRUE) 
cellsn2.kcmds = cmdscale(dist2kn, k = 20, eig = TRUE)

## -----------------------------------------------------------------------------
#| label: fig-KernelMDSplots
#| results: hide
#| layout-nrow: 1
#| fig-cap: "Kernel multidimensional scaling."
#| fig-subcap:
#|   - "$\\text{}$"
#|   - "$\\text{}$"
colc = rowData(Moignard)$cellcol
library("scatterplot3d")
scatterplot3d(cellsn2.kcmds$points[, 1:3], color=colc, pch = 20,
   xlab = "Axis k1", ylab = "Axis k2", zlab = "Axis k3", angle=15)
scatterplot3d(cellsn2.kcmds$points[, 1:3], color=colc, pch = 20,
   xlab = "Axis k1", ylab = "Axis k2", zlab = "Axis k3", angle = -70)

## -----------------------------------------------------------------------------
#| label: CodeB
#| results: hide
library("rgl")
plot3d(cellsn2.kcmds$points[, 1:3], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis3")
plot3d(cellsn2.kcmds$points[, c(1,2,4)], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis4")
# Using an L1 distance instead.
plot3d(cellsn1.kcmds$points[, 1:3], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis3")
plot3d(cellsn1.kcmds$points[, c(1,2,4)], col = colc, size = 3,
       xlab = "Axis1", ylab = "Axis2", zlab = "Axis4")

## -----------------------------------------------------------------------------
#| label: lpc3d
library("LPCM")
library("diffusionMap")
dmap1 = diffuse(dist1n.l1, neigen = 10)
combs = combn(4, 3)
lpcplots = apply(combs, 2, function(j) lpc(dmap1$X[, j], scale = FALSE))

## -----------------------------------------------------------------------------
#| label: PLOT3DLPC
#| output: false
library("rgl")
for (i in seq_along(lpcplots))
  plot(lpcplots[[i]], type = "l", lwd = 3,
  xlab = paste("Axis", combs[1, i]),
  ylab = paste("Axis", combs[2, i]),
  zlab = paste("Axis", combs[3, i]))

## -----------------------------------------------------------------------------
#| label: fig-diffusionmap3
#| echo: false
#| column: margin
#| layout-nrow: 1
    c('TripleArm.png','SmoothLineP134h7.png')))

## -----------------------------------------------------------------------------
#| label: scatter3smoothedline
#| fig-show: hide
outlpce134 = lpc(dmap1$X[,c(1,3,4)], scale=FALSE, h=0.5)
plot3d(dmap1$X[,c(1,3,4)], col=colc, pch=20, 
       xlab="Axis1", ylab="Axis3", zlab="Axis4")
plot3d(outlpce134$LPC, type="l", lwd=7, add=TRUE)

outlpce134 = lpc(dmap1$X[,c(1,3,4)], scale=FALSE, h=0.7)
plot3d(outlpce134$LPC, type="l", lwd=7,
       xlab="Axis1", ylab="Axis3", zlab="Axis4")
plot3d(dmap1$X[,c(1,3,4)], col=colc, 
       xlab="", ylab="", zlab="", add=TRUE)

## -----------------------------------------------------------------------------
#| label: fig-dmap134
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: Diffuse
#| fig-show: hide
library("diffusionMap")
dmap2 = diffuse(dist2n.euclid, neigen = 11)
dmap1 = diffuse(dist1n.l1, neigen = 11)
plot(dmap2)

## -----------------------------------------------------------------------------
#| label: scp3d
#| fig-show: hide
library("scatterplot3d")
scp3d = function(axestop = 1:3, dmapRes = dmap1, color = colc,
           anglea = 20, pch = 20)
scatterplot3d(dmapRes$X[, axestop], color = colc,
    xlab = paste("Axis",axestop[1]), ylab = paste("Axis", axestop[2]),
    zlab = paste("Axis",axestop[3]), pch = pch, angle = anglea)

## -----------------------------------------------------------------------------
#| label: dmap3dplots-eval
#| output: false
#| echo: false
# TODO: duplicate of the below as a workaround;
# else white space is rendered between lines
scp3d()
scp3d(anglea=310)
scp3d(anglea=210)
scp3d(anglea=150)

## -----------------------------------------------------------------------------
#| label: dmap3dplots-show
#| eval: false
## scp3d()
## scp3d(anglea=310)
## scp3d(anglea=210)
## scp3d(anglea=150)

## -----------------------------------------------------------------------------
#| label: CodeD
#| fig-show: hide
# interactive plot
library("rgl")
plot3d(dmap1$X[,1:3], col=colc, size=3)
plot3d(dmap1$X[,2:4], col=colc, size=3)
