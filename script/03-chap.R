
## -----------------------------------------------------------------------------
#| label: introplot
#| eval: false
#| echo: false
## ## produces the xkcd-plotting image used below
## #library("xkcd")
## #library("showtext")
## #library("sysfonts")
## #library("tibble")
## 
## introplotdata = tibble(
##   y = c(seq(-8, 1, length=25)^2, rep(1, 5), seq(1, 5,length=25)^2)^2,
##   x = seq(1, 55, length.out = length(y)))
## 
## dataman = tibble(
##   x = 30,
##   y = 400,
##   scale = 100,
##   ratioxy = 0.1,
##   angleofspine =  -pi/2 ,
##   anglerighthumerus = -pi/6,
##   anglelefthumerus  = pi * 7/6,
##   anglerightradius = 0,
##   angleleftradius = 0,
##   angleleftleg  = 19*pi/12,
##   anglerightleg = 17*pi/12,
##   angleofneck   = 1.4*pi)
## 
## mapping = do.call(aes_string, colnames(dataman) %>% (function(x) setNames(as.list(x), x)))
## 
## ggplot(introplotdata) + geom_line(aes(x = x, y = y), size = 2) +
##    xkcdaxis(c(0, 50), c(0, 1000)) + xlab("Time to make plot in minutes") +
##    ylab("Time to understand plot in minutes") + xkcdman(mapping, dataman) +
##    theme(axis.title.x = element_text(margin = margin(15, 0, 0, 0)))

## -----------------------------------------------------------------------------
#| label: xkcd-plotting
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-graphics-plotter
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-graphics-basicplotting1
#| column: margin
#| fig-cap: "Plot of concentration vs. density for an ELISA assay of DNase."
head(DNase)
plot(DNase$conc, DNase$density)

## -----------------------------------------------------------------------------
#| label: fig-graphics-basicplotting2
#| column: margin
#| fig-cap: "Same data as in @fig-graphics-basicplotting1 but with better axis labels and a different plot symbol."
plot(DNase$conc, DNase$density,
  ylab = attr(DNase, "labels")$y,
  xlab = paste(attr(DNase, "labels")$x, attr(DNase, "units")$x),
  pch = 3,
  col = "blue")

## -----------------------------------------------------------------------------
#| label: fig-graphics-basicplotting3
#| column: margin
#| layout-ncol: 1
#| fig-width: 3
#| fig-height: 3
#| fig-margin: false
#| fig-cap: "\\(a) Histogram of the density from the ELISA assay, and (b) boxplots of these values stratified by the assay run. The boxes are ordered along the axis in lexicographical order because the runs were stored as text strings. We could use R's type conversion functions to achieve numerical ordering."
#| fig-subcap:
#|   - ""
#|   - ""
hist(DNase$density, breaks=25, main = "")
boxplot(density ~ Run, data = DNase)

## -----------------------------------------------------------------------------
#| label: fig-graphics-cells
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: loadHiiragi
#| cache: false
library("Hiiragi2013")
data("x")
dim(Biobase::exprs(x))

## -----------------------------------------------------------------------------
#| label: xpData
head(pData(x), n = 2)

## -----------------------------------------------------------------------------
#| label: groupSize
library("dplyr")
groups = group_by(pData(x), sampleGroup) |>
  summarise(n = n(), color = unique(sampleColour))
groups

## -----------------------------------------------------------------------------
#| label: explainpipe
#| eval: false
## f(x) |> g(y) |> h()
## h(g(f(x), y))

## -----------------------------------------------------------------------------
#| label: fig-graphics-figredobasicplottingwithggplot
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "Our first **[ggplot2](https://cran.r-project.org/web/packages/ggplot2/)** figure, similar to the base graphics @fig-graphics-basicplotting1."
library("ggplot2")
ggplot(DNase, aes(x = conc, y = density)) + geom_point()

## -----------------------------------------------------------------------------
#| label: fig-graphics-qplot1
#| column: margin
#| fig-width: 5
#| fig-height: 4
#| fig-cap: "A barplot, produced with the `ggplot` function from the table of group sizes in the mouse single cell data."
ggplot(groups, aes(x = sampleGroup, y = n)) +
  geom_bar(stat = "identity")

## -----------------------------------------------------------------------------
#| label: checkgeombar
#| echo: false
## check an assertion made in the text above
stopifnot(formals(ggplot2::geom_bar)$stat=="count")

## -----------------------------------------------------------------------------
#| label: groupColor
groupColor = setNames(groups$color, groups$sampleGroup)

## -----------------------------------------------------------------------------
#| label: fig-graphics-qplot2
#| column: margin
#| fig-width: 5
#| fig-height: 4
#| fig-cap: "Similar to @fig-graphics-qplot1, but with colored bars and better bar labels."
ggplot(groups, aes(x = sampleGroup, y = n, fill = sampleGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = groupColor, name = "Groups") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## -----------------------------------------------------------------------------
#| child: devil.qmd

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: ggplotobject
gg = ggplot(DNase, aes(x = conc, y = density)) + geom_point()

## -----------------------------------------------------------------------------
#| label: ggpprintobject
#| eval: false
## gg
## print(gg)

## -----------------------------------------------------------------------------
#| label: plotsave1
ggplot2::ggsave("DNAse-histogram-demo.pdf", plot = gg)

## -----------------------------------------------------------------------------
#| label: plotsave2
#| echo: false
#| results: hide
file.remove("DNAse-histogram-demo.pdf")

## -----------------------------------------------------------------------------
#| label: loadlib
#| echo: false
library("mouse4302.db")

## -----------------------------------------------------------------------------
#| label: findprobepairs
#| eval: false
#| echo: false
## # I used this code to find the below two probes
## idx = order(rowVars(Biobase::exprs(x)), decreasing=TRUE)[seq_len(2000)]
## cc  = cor(t(Biobase::exprs(x)[idx,]))
## cco = order(cc)[seq(1, 1001, by=2) ]
## jj2 = rownames(Biobase::exprs(x))[ idx[ (cco-1) %/% length(idx) + 1 ] ]
## jj1 = rownames(Biobase::exprs(x))[ idx[ (cco-1) %%  length(idx) + 1 ] ]
## dftx = as.data.frame(t(Biobase::exprs(x)))
## par(ask = TRUE)
## for(i in seq(along = cco)) {
##   df = AnnotationDbi::select(mouse4302.db,
##    keys = c(jj1[i], jj2[i]), keytype = "PROBEID",
##    columns = c("SYMBOL", "GENENAME"))
##   print(ggplot(dftx, aes( x = get(jj1[i]), y = get(jj2[i]))) +
##   geom_point(shape = 1) +
##   xlab(paste(jj1[i], df$SYMBOL[1])) +
##   ylab(paste(jj2[i], df$SYMBOL[2])) +
##   ggtitle(round(cc[jj1[i], jj2[i]], 3)) +
##   geom_smooth(method = "loess"))
## }

## -----------------------------------------------------------------------------
#| label: fig-graphics-scp2layers1
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "A scatterplot with three layers that show different statistics of the same data: points (`geom_point`), a smooth regression line and a confidence band (the latter two from `geom_smooth`)."
dftx = data.frame(t(Biobase::exprs(x)), pData(x))
ggplot( dftx, aes( x = X1426642_at, y = X1418765_at )) +
  geom_point( shape = 1 ) +
  geom_smooth( method = "loess" )

## -----------------------------------------------------------------------------
#| label: checkclassdftx
#| echo: false
stopifnot(is(dftx, "data.frame"))

## -----------------------------------------------------------------------------
#| echo: false
.one <- AnnotationDbi::select(mouse4302.db, 
                              keys = "1418765_at", 
                              keytype = "PROBEID", 
                              columns = "SYMBOL")$SYMBOL
.two <- AnnotationDbi::select(mouse4302.db, 
                              keys = "1426642_at", 
                              keytype = "PROBEID", 
                              columns = "SYMBOL")$SYMBOL

## -----------------------------------------------------------------------------
#| label: fig-graphics-scp2layers2
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: !expr paste("As @fig-graphics-scp2layers1, but in addition with points colored by the time point and cell lineage (as defined in @fig-graphics-qplot2). We can now see that the expression values of the gene", .one, "(targeted by the probe 1418765_at, along the y-axis) are consistently high in the early time points, whereas its expression goes down in the EPI samples at days 3.5 and 4.5. In the FGF4-KO, this decrease is delayed - at E3.5, its expression is still high. Conversely, the gene", .two, "(1426642_at, x-axis) is off in the early timepoints and then goes up at days 3.5 and 4.5. The PE samples (green) show a high degree of cell-to-cell variability.")
ggplot(dftx, aes(x = X1426642_at, y = X1418765_at))  +
  geom_point(aes(color = sampleGroup), shape = 19) +
  scale_color_manual(values = groupColor, guide = "none") +
  geom_smooth(method = "loess")

## -----------------------------------------------------------------------------
#| child: devil.qmd

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: mouse4302.db
#| results: hide
library("mouse4302.db")

## -----------------------------------------------------------------------------
#| label: select
AnnotationDbi::select(mouse4302.db,
   keys = c("1426642_at", "1418765_at"), keytype = "PROBEID",
   columns = c("SYMBOL", "GENENAME"))

## -----------------------------------------------------------------------------
#| label: fig-graphics-hists
#| column: margin
#| fig-width: 3.5
#| fig-height: 2.5
#| fig-cap: "Histogram of probe intensities for one particular sample, cell number 20, which was from day E3.25."
dfx = as.data.frame(Biobase::exprs(x))
ggplot(dfx, aes(x = `20 E3.25`)) + geom_histogram(binwidth = 0.2)

## -----------------------------------------------------------------------------
#| label: figbpgg1
pb = ggplot(groups, aes(x = sampleGroup, y = n))

## -----------------------------------------------------------------------------
#| label: fig-graphics-figbpempty
#| column: margin
#| fig-width: 3.2
#| fig-height: 2.5
#| fig-cap: "`pb`: without a geometric object, the plot remains empty."
class(pb)
pb

## -----------------------------------------------------------------------------
#| label: fig-graphics-bpgg3
#| column: margin
#| fig-width: 5
#| fig-height: 4
#| fig-cap: "The graphics object `bp` in its full glory."
pb = pb + geom_bar(stat = "identity")
pb = pb + aes(fill = sampleGroup)
pb = pb + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pb = pb + scale_fill_manual(values = groupColor, name = "Groups")
pb

## -----------------------------------------------------------------------------
#| label: fig-graphics-bpgg7
#| column: margin
#| fig-width: 5
#| fig-height: 4
#| fig-cap: "A barplot in a polar coordinate system."
pb.polar = pb + coord_polar() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") + ylab("")
pb.polar

## -----------------------------------------------------------------------------
#| label: genes2ps1
selectedProbes = c( Fgf4 = "1420085_at", Gata4 = "1418863_at",
                   Gata6 = "1425463_at",  Sox2 = "1416967_at")

## -----------------------------------------------------------------------------
#| label: genes2ps2
#| eval: false
#| echo: false
## # How I found the selectedProbes:
## AnnotationDbi::select(mouse4302.db,
##    keys = c("Fgf4", "Sox2", "Gata6", "Gata4"), keytype = "SYMBOL",
##    columns = c("PROBEID"))

## -----------------------------------------------------------------------------
#| label: genes2ps3
#| echo: false
selectedProbes2 = AnnotationDbi::select(mouse4302.db,
   keys = selectedProbes, keytype = "PROBEID", columns = c("SYMBOL"))
stopifnot(identical(sort(selectedProbes2$SYMBOL), sort(names(selectedProbes))),
          all(selectedProbes[selectedProbes2$SYMBOL] == selectedProbes2$PROBEID))

## -----------------------------------------------------------------------------
#| label: melt
library("reshape2")
genes = melt(Biobase::exprs(x)[selectedProbes, ],
             varnames = c("probe", "sample"))

## -----------------------------------------------------------------------------
#| label: symbol
genes$gene =
  names(selectedProbes)[match(genes$probe, selectedProbes)]
head(genes)

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedbp1
#| column: margin
#| fig-width: 3
#| fig-height: 3.75
#| fig-cap: "Barplots showing the means of the distributions of expression measurements from four probes."
ggplot(genes, aes(x = gene, y = value)) +
  stat_summary(fun = mean, geom = "bar")

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedbp2
#| column: margin
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: "Barplots with error bars indicating standard error of the mean."
library("Hmisc")
ggplot(genes, aes( x = gene, y = value, fill = gene)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               width = 0.25)

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedboxpl
#| column: margin
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: "Boxplots."
p = ggplot(genes, aes( x = gene, y = value, fill = gene))
p + geom_boxplot()

## -----------------------------------------------------------------------------
#| label: fig-graphics-oneddot
#| layout-nrow: 1
#| fig-width: 5
#| fig-height: 5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "\\(a) Dot plots, made using `geom_dotplot` from **[ggplot2](https://cran.r-project.org/web/packages/ggplot2/)**. (b) Beeswarm plots, made using `geom_beeswarm` from **[ggbeeswarm](https://cran.r-project.org/web/packages/ggbeeswarm/)**."
#| fig-subcap:
#|   - ""
#|   - ""
p + geom_dotplot(binaxis = "y", binwidth = 1/6,
       stackdir = "center", stackratio = 0.75,
       aes(color = gene))
library("ggbeeswarm")
p + geom_beeswarm(aes(color = gene))

## -----------------------------------------------------------------------------
#| label: fig-graphics-oneddens
#| column: margin
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: "Density plots."
ggplot(genes, aes( x = value, color = gene)) + geom_density()

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedviolin
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "Violin plots."
p + geom_violin()

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedridge4
#| column: margin
#| fig-width: 3.5
#| fig-height: 3
#| fig-cap: "Ridgeline plots."
library("ggridges")
ggplot(genes, aes(x = value, y = gene, fill = gene)) + 
  geom_density_ridges()

## -----------------------------------------------------------------------------
#| label: onedridge42-show
#| warning.known: !expr c("rows containing non-finite values")
#| fig-show: hide
top42 = order(rowMeans(Biobase::exprs(x)), decreasing = TRUE)[1:42]
g42 = melt(Biobase::exprs(x)[rev(top42), ], varnames = c("probe", "sample"))
ggplot(g42, aes(x = value, y = probe)) 

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedridge42
#| echo: false
#| warning.known: !expr c("rows containing non-finite values")
#| column: margin
#| fig-height: 4.5
#| fig-width: 3.5
#| fig-cap: "Like @fig-graphics-onedridge4, with more genes."
ggplot(g42, aes(x = value, y = probe)) + 
  geom_density_ridges() + theme(legend.position = "none",
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) + xlim(13, 15) 

## -----------------------------------------------------------------------------
#| label: fig-graphics-ecdfexample
#| column: margin
#| out-width: 75%
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "Sorted values of `simdata` versus their index. This is the empirical cumulative distribution function of `simdata`."
simdata = rnorm(70)
tibble(index = seq(along = simdata),
          sx = sort(simdata)) %>%
ggplot(aes(x = sx, y = index)) + geom_step()

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedecdf
#| column: margin
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-cap: "Empirical cumulative distribution functions (ECDF)."
ggplot(genes, aes( x = value, color = gene)) + stat_ecdf()

## -----------------------------------------------------------------------------
#| label: fig-graphics-Lawrence-TCGA-Nature-2013-Fig1
#| echo: false
#| out-width: 75%
#| fig-margin: false
#| fig-cap-location: margin

## -----------------------------------------------------------------------------
#| label: modes
#| eval: false
#| echo: false
## # I used the functon bimodality_coefficient from the modes package to identify the most
## # bimodal looking array, number 64
## j0 = which.max(vapply(seq_len(ncol(x)), function(j){
##        modes  ::    bimodality_coefficient(Biobase::exprs(x)[, j])
##     }, numeric(1)))

## -----------------------------------------------------------------------------
#| label: fig-graphics-onedtrsf
#| warning.known: !expr c("rows containing non-finite values", "rows containing missing values")
#| layout-nrow: 1
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Histograms of the same data, with and without logarithm transform. (a) The data are shown on the scale on which they are stored in the data object `x`, which resulted from logarithm (base 2) transformation of the microarray fluorescence intensities [@Irizarry:Biostat:2003]; (b) after re-exponentiating them back to the fluorescence scale. For better use of space, we capped the $x$-axis range at 1500."
#| fig-subcap:
#|   - ""
#|   - ""
ggplot(dfx, aes(x = `64 E4.5 (EPI)`)) + geom_histogram(bins = 100)
ggplot(dfx, aes(x = 2 ^ `64 E4.5 (EPI)`)) + 
  geom_histogram(binwidth = 20) + xlim(0, 1500)

## -----------------------------------------------------------------------------
#| label: fig-graphics-twodsp1
#| column: margin
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: !expr paste("Scatterplot of", nrow(dfx), "expression measurements for two of the samples.")
scp = ggplot(dfx, aes(x = `59 E4.5 (PE)` ,
                      y = `92 E4.5 (FGF4-KO)`))
scp + geom_point()

## -----------------------------------------------------------------------------
#| label: fig-graphics-twodsp2
#| column: margin
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: "As @fig-graphics-twodsp1, but with semi-transparent points to resolve some of the overplotting."
scp  + geom_point(alpha = 0.1)

## -----------------------------------------------------------------------------
#| label: fig-graphics-twodsp3
#| column: margin
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: "As @fig-graphics-twodsp1, but rendered as a contour plot of the 2D density estimate."
scp + geom_density2d()

## -----------------------------------------------------------------------------
#| label: fig-graphics-twodsp4
#| layout-nrow: 1
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap: "2D density plots. (a) As @fig-graphics-twodsp3, but with smaller smoothing bandwidth and tighter binning for the contour lines. (b) With color filling."
#| fig-subcap:
#|   - ""
#|   - ""
scp + geom_density2d(h = 0.5, bins = 60)
library("RColorBrewer")
colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100))))

scp + stat_density2d(h = 0.5, bins = 60,
          aes( fill = after_stat(level)), geom = "polygon") +
  colorscale + coord_fixed()

## -----------------------------------------------------------------------------
#| label: fig-graphics-twodsp6
#| layout-nrow: 1
#| fig-width: 3.5
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Hexagonal binning. (a) Default parameters. (b) Finer bin sizes and customized color scale."
#| fig-subcap:
#|   - ""
#|   - ""
scp + geom_hex() + coord_fixed()
scp + geom_hex(binwidth = c(0.2, 0.2)) + colorscale +
  coord_fixed()

## -----------------------------------------------------------------------------
#| label: fig-graphics-banking
#| column: margin
#| layout-ncol: 1
#| fig-width: 3.75
#| fig-height: 3.75
#| fig-cap: "The sunspot data. In (a), the plot shape is roughly quadratic, a frequent default choice. In (b), a technique called **banking** was used to choose the plot shape. (Note: the placement of the tick labels is not great in this plot and would benefit from customization.)"
#| fig-subcap:
#|   - ""
#|   - ""
library("ggthemes")
sunsp = tibble(year   = time(sunspot.year),
               number = as.numeric(sunspot.year))
sp = ggplot(sunsp, aes(x = year, y = number)) + geom_line()
sp
ratio = with(sunsp, bank_slopes(year, number))
sp + coord_fixed(ratio = ratio)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Sometimes this is also called *trellis* or *lattice* graphics, in an allusion to how these arrays of plots look like. The first major R package to implement faceting was **[lattice](https://cran.r-project.org/web/packages/lattice/)**. In this book, we'll use the faceting functionalities provided through **[ggplot2](https://cran.r-project.org/web/packages/ggplot2/)**."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "The first three lines of this code chunk are not strictly necessary – they're just reformatting the `lineage` column of the `dftx` dataframe to make it more concise."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-graphics-facet1
#| fig-width: 8
#| fig-height: 2
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "An example of **faceting**: the same data as in @fig-graphics-scp2layers1, but now split by the categorical variable `lineage`."
library("magrittr")
dftx$lineage %<>% sub("^$", "no", .)
dftx$lineage %<>% factor(levels = c("no", "EPI", "PE", "FGF4-KO"))

ggplot(dftx, aes(x = X1426642_at, y = X1418765_at)) +
  geom_point() + facet_grid( . ~ lineage )

## -----------------------------------------------------------------------------
#| label: fig-graphics-facet2
#| out-width: 75%
#| fig-width: 8
#| fig-height: 6
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "**Faceting**: the same data as in @fig-graphics-scp2layers1, split by the categorical variables `Embryonic.day` (rows) and `lineage` (columns)."
ggplot(dftx,
  aes(x = X1426642_at, y = X1418765_at)) + geom_point() +
   facet_grid( Embryonic.day ~ lineage )

## -----------------------------------------------------------------------------
#| label: fig-graphics-facet3
#| column: margin
#| fig-width: 4
#| fig-height: 4
#| fig-cap: "**Faceting**: the same data as in @fig-graphics-scp2layers1, split by the continuous variable `X1450989_at` and arranged by `facet_wrap`."
ggplot(mutate(dftx, Tdgf1 = cut(X1450989_at, breaks = 4)),
   aes(x = X1426642_at, y = X1418765_at)) + geom_point() +
   facet_wrap( ~ Tdgf1, ncol = 2 )

## -----------------------------------------------------------------------------
#| label: plotly
#| eval: false
## library("plotly")
## plot_ly(economics, x = ~ date, y = ~ unemploy / pop)

## -----------------------------------------------------------------------------
#| label: fig-graphics-rglvolcano
#| column: margin
#| echo: false

## -----------------------------------------------------------------------------
#| label: volcano
data("volcano")
volcanoData = list(
  x = 10 * seq_len(nrow(volcano)),
  y = 10 * seq_len(ncol(volcano)),
  z = volcano,
  col = terrain.colors(500)[cut(volcano, breaks = 500)]
)
library("rgl")
with(volcanoData, persp3d(x, y, z, color = col))

## -----------------------------------------------------------------------------
#| label: volcanocheck
#| echo: false
.volcanocut = cut(volcano, breaks = 500)
stopifnot(!any(is.na(.volcanocut)), all(as.integer(.volcanocut) %in% 1:500))

## -----------------------------------------------------------------------------
#| label: fig-graphics-simplecolorpie
#| column: margin
#| echo: !expr -c(1)
#| fig-width: 3
#| fig-height: 3
#| fig-cap: "The first eight colors in the base R color palette."
par(mai = rep(0, 4))
pie(rep(1, 8), col=1:8)

## -----------------------------------------------------------------------------
#| label: fig-graphics-ggplot2colorpie
#| column: margin
#| fig-width: 4
#| fig-height: 3
#| fig-cap: "The first eight colors in the **[ggplot2](https://cran.r-project.org/web/packages/ggplot2/)** color palette."
ggplot(tibble(u = factor(1:8), v = 1), 
       aes(x = "",  y = v, fill = u)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) + theme_void()

## -----------------------------------------------------------------------------
#| label: fig-graphics-RColorBrewer
#| echo: !expr -c(1)
#| column: margin
#| fig-width: 6
#| fig-height: 10
#| fig-cap: "RColorBrewer palettes."
par(mai = c(0, 0.8, 0, 0))
display.brewer.all()

## -----------------------------------------------------------------------------
#| label: color3
head(brewer.pal.info)
table(brewer.pal.info$category)

## -----------------------------------------------------------------------------
#| label: color4
brewer.pal(4, "RdYlGn")

## -----------------------------------------------------------------------------
#| label: fig-graphics-colorRampPalette
#| echo: !expr -c(1)
#| column: margin
#| fig-width: 3
#| fig-height: 0.7
#| fig-cap: "A quasi-continuous color palette derived by interpolating between the colors `darkorange3`, `white` and `darkblue`."
par(mai = rep(0.1, 4))
mypalette  = colorRampPalette(
    c("darkorange3", "white","darkblue")
  )(100)
head(mypalette)
image(matrix(1:100, nrow = 100, ncol = 10), col = mypalette,
        xaxt = "n", yaxt = "n", useRaster = TRUE)

## -----------------------------------------------------------------------------
#| label: fig-graphics-heatmap
#| fig-width: 7
#| fig-height: 7
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "A heatmap of relative expression values, i.e., logarithmic fold change compared to the average expression of that gene (row) across all samples (columns). The color scale uses a diverging palette whose midpoint is at 0."
library("pheatmap")
topGenes = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:500]
rowCenter = function(x) { x - rowMeans(x) }
pheatmap(rowCenter(dfx[topGenes, ]), 
  show_rownames = FALSE, 
  show_colnames = FALSE, 
  breaks = seq(-5, +5, length = 101),
  annotation_col = pData(x)[, c("sampleGroup", "Embryonic.day", "ScanDate", "genotype") ],
  annotation_colors = list(
    sampleGroup = groupColor,
    genotype = c(`FGF4-KO` = "chocolate1", `WT` = "azure2"),
    Embryonic.day = setNames(brewer.pal(9, "Blues")[c(3, 6, 9)], c("E3.25", "E3.5", "E4.5")),
    ScanDate = setNames(brewer.pal(nlevels(x$ScanDate), "YlGn"), levels(x$ScanDate))
  )
)

## -----------------------------------------------------------------------------
#| label: groupColor2
groupColor[1]

## -----------------------------------------------------------------------------
#| label: hexvals
#| echo: false
hexvals = sapply(1:3, function(i) substr(groupColor[1], i*2, i*2+1))
decvals = strtoi(paste0("0x", hexvals))

## -----------------------------------------------------------------------------
#| label: somecolors
#| echo: false
#| results: hide
library("colorspace")
library("grid")

plothcl = function(h, c, l, what, x0 = 0.5, y0 = 0.5, default.units = "npc", ...) {
  switch(what,
         "c" = {
           stopifnot(length(l)==1)
           n = length(c)
         },
         "l" = {
           stopifnot(length(c)==1)
           n = length(l)
         },
         stop("Sapperlot"))

  cr = seq(0.1, 0.5, length = n+1)
  dr = 0.05 / n

  for (j in seq_len(n)) {
    r = c(cr[j]+dr, cr[j+1]-dr)
    for(i in 1:(length(h)-1)){
      phi = seq(h[i], h[i+1], by=1)/180*pi
      px = x0 + c(r[1]*cos(phi), r[2]*rev(cos(phi)))
      py = y0 + c(r[1]*sin(phi), r[2]*rev(sin(phi)))
      mycol = switch(what,
        "c" = hcl(h=mean(h[i+(0:1)]), c=c[j], l=l),
        "l" = hcl(h=mean(h[i+(0:1)]), c=c, l=l[j]))
      grid::grid.polygon(px, py, 
                         gp=gpar(col=mycol, fill=mycol),
                         default.units=default.units,...)
    }
  }
}

## -----------------------------------------------------------------------------
#| label: fig-graphics-hcl
#| echo: false
#| layout-nrow: 1
#| fig-asp: 1
#| fig-height: 3.5
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Circles in HCL colorspace. (a) Luminance $L$ fixed at $75$, while the angular coordinate $H$ (hue) varies from 0 to 360 and the radial coordinate $C$ takes the values $0, 10, ..., 60$. (b) Constant chroma $C$ at $50$, $H$ as above, and the radial coordinate is luminance $L$ at values $10, 20, ..., 90$."
#| fig-subcap:
#|   - ""
#|   - ""
plothcl( h = seq(0, 360, by=3), c = seq(5, 75, by=10), l = 75,   what="c")
grid.newpage()
plothcl( h = seq(0, 360, by=3), c = 55, l = seq(20, 100, by=10), what="l")

## -----------------------------------------------------------------------------
#| label: fig-graphics-MA
#| layout-nrow: 1
#| fig-width: 3
#| fig-height: 3
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "The effect of rank transformation on the visual perception of dependency."
#| fig-subcap:
#|   - ""
#|   - ""
gg = ggplot(tibble(A = Biobase::exprs(x)[, 1], M = rnorm(length(A))),
            aes(y = M))
gg + geom_point(aes(x = A), size = 0.2)
gg + geom_point(aes(x = rank(A)), size = 0.2)

## -----------------------------------------------------------------------------
#| label: fig-graphics-otherfont
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-graphics-mathnot
#| column: margin
#| fig-height: 2.5
#| fig-cap: "Volume $\\Omega$ of the $\\nu$-dimensional sphere with radius $\\rho=1$, for $\\nu=1,...,15$."
volume = function(rho, nu)
            pi^(nu/2) * rho^nu / gamma(nu/2+1)

ggplot(tibble(nu    = 1:15,
  Omega = volume(1, nu)), aes(x = nu, y = Omega)) +
geom_line() +
xlab(expression(nu)) + ylab(expression(Omega)) +
geom_text(label =
"Omega(rho,nu)==frac(pi^frac(nu,2)~rho^nu, Gamma(frac(nu,2)+1))",
  parse = TRUE, x = 6, y = 1.5)

## -----------------------------------------------------------------------------
#| label: fig-graphics-timesfont
#| column: margin
#| fig-height: 3
#| fig-cap: "As @fig-graphics-onedecdf, with a different font."
ggplot(genes, aes( x = value, color = gene)) + stat_ecdf() +
  theme(text = element_text(family = "Times"))

## -----------------------------------------------------------------------------
#| label: otherfont
#| eval: false
#| echo: false
## ggplot(genes, aes( x = value, color = gene)) + stat_ecdf() + theme(text = element_text(family = "Bauhaus 93"))

## -----------------------------------------------------------------------------
#| label: fig-graphics-EBI-genomebrowser-rnaseq
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-graphics-ideogram-1
#| column: margin
#| fig-height: 1.5
#| fig-cap: "Chromosome 1 of the human genome: ideogram plot."
library("ggbio")
data("hg19IdeogramCyto", package = "biovizBase")
plotIdeogram(hg19IdeogramCyto, subchr = "chr1")

## -----------------------------------------------------------------------------
#| label: fig-graphics-darned1
#| warning.known: !expr c("need valid seqlengths information for accurate mapping")
#| column: margin
#| fig-cap: "Karyogram with RNA editing sites. `exReg` indicates whether a site is in the coding region (C), 3'- or 5'-UTR."
library("GenomicRanges")
data("darned_hg19_subset500", package = "biovizBase")
autoplot(darned_hg19_subset500, layout = "karyogram",
         aes(color = exReg, fill = exReg))

## -----------------------------------------------------------------------------
#| label: fig-graphics-darned2
#| out-width: 50%
#| fig-width: 5
#| fig-height: 6
#| fig-cap: "Improved version of @fig-graphics-darned1."
data("ideoCyto", package = "biovizBase")
dn = darned_hg19_subset500
seqlengths(dn) = seqlengths(ideoCyto$hg19)[names(seqlengths(dn))]
dn = keepSeqlevels(dn, paste0("chr", c(1:22, "X")))
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg))

## -----------------------------------------------------------------------------
#| label: whatisdarned1
darned_hg19_subset500[1:2,]

## -----------------------------------------------------------------------------
#| label: whatisdarned2
#| echo: false
stopifnot(is(darned_hg19_subset500, "GRanges"), identical(start(darned_hg19_subset500),end(darned_hg19_subset500)))

## -----------------------------------------------------------------------------
#| label: theme_bw
#| eval: false
#| fig-width: 3
#| fig-height: 3
## ggcars = ggplot(mtcars, aes(x = hp, y = mpg)) + geom_point()
## ggcars
## ggcars + theme_bw()
## ggcars + theme_minimal()

## -----------------------------------------------------------------------------
#| label: xkcdgraph
#| echo: false
