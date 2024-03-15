
## -----------------------------------------------------------------------------
#| label: Darwin-Tree-1837-web
#| echo: false

## ----warning=FALSE, message=FALSE---------------------------------------------
#| label: fig-chemo
#| echo: false
#| column: margin
#| fig-cap: "A small protein-protein network that represents pairwise relationships between proteins."
# Changed everywhere ../data/ to file.path(DATA, )
#dats = read.table("../data/small_chemokine.txt", header = TRUE)
dats = read.table(file.path(DATA, "small_chemokine.txt"), header = TRUE)
library("ggtree")
library("igraph")
library("tibble")
library("ggplot2")
library("reshape")
library("ggraph")
library("ggrepel")
library("tidygraph")
library("dplyr")
gr = graph_from_data_frame(dats[,c("node1", "node2")], directed = FALSE)
E(gr)$weight = 1
V(gr)$size = centr_degree(gr)$res
ggd<-ggraph(gr, layout = 'fr') + 
  geom_edge_link(color="black",  alpha=1/2, linewidth = 1) + 
  geom_node_point(size=3,alpha=1/2, color="orange") +
  geom_node_text(aes(label = vertex_attr(gr)$name), color="#8856a7",size=6, repel = TRUE) +
  theme_bw() + theme(legend.position="none")

ggd

## -----------------------------------------------------------------------------
#| label: fromtotableexample
#| echo: FALSE
data.frame(from=c("A","B","A","C","E"), to=c("B","C","E","D","F"))

## -----------------------------------------------------------------------------
#| label: fig-igraphplot
#| fig-margin: false
#| fig-width: 6
#| fig-height: 6
#| out-width: 50%
#| fig-cap: "A small undirected graph with numbered nodes."
library("igraph")
edges = matrix(c(1,3, 2,3, 3,4, 4,5, 4,6), byrow = TRUE, ncol = 2)
g1 = graph_from_edgelist(edges, directed = FALSE)
vertex_attr(g1, name = "name") = 1:6
plot(g1, vertex.size = 25, edge.width = 5, vertex.color = "coral")

## -----------------------------------------------------------------------------
#| label: adjacencyplot1
#| echo: false
ggplotadjacency = function(a) {
  n = nrow(a)
  p = ncol(a)
  fromto  = reshape2::melt(a)
  stopifnot(identical(nrow(fromto), n*p))
  fromto$value = as.factor(fromto$value)
  cols = c("white", "darkblue")
  ggplot(data = fromto, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(colour = "black") +
    coord_fixed(ratio = 1, ylim = c(0.5, n + 0.5), xlim = c(0.5, p + 0.5)) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(name = "" , breaks = 1:p, labels = paste(1:p)) +
    scale_y_reverse(  name = "" , breaks = n:1, labels = paste(n:1)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 14),
      legend.key = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "white"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank() 
    )
}

## -----------------------------------------------------------------------------
#| label: fig-adjmatrix
#| echo: false
#| out-width: 50%
#| fig-width: 4.5
#| fig-height: 4
#| fig-margin: false
#| fig-cap: "The adjacency matrix of the graph shown in @fig-igraphplot is is a symmetric $n \\times n$ matrix of $0$s and $1$s, where is $n$ is the number of nodes."
ggplotadjacency(as_adj(g1, sparse = FALSE))

## -----------------------------------------------------------------------------
#| label: graphex1wh
edges = "1,3\n2,3\n3,4\n4,6\n4,5"
df = read.csv(textConnection(edges), header = FALSE)
sg = graph_from_data_frame(df, directed = FALSE)
sg

## -----------------------------------------------------------------------------
#| label: fig-bipartite
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-finches-1
#| out-width: 85%
#| fig-width: 9
#| fig-height: 6
#| fig-margin: false
#| fig-cap: "The `finches` graph. There are many ways to improve the layout, including better taking into account the bipartite nature of the graph."
#
# Replace the data directory with the correct path
#finch = readr::read_csv("../data/finch.csv", comment = "#", col_types = "cc")
finch = readr::read_csv(file.path(DATA, "finch.csv"), comment = "#", col_types = "cc")
finch
library("network")
finch.nw  = as.network(finch, bipartite = TRUE, directed = FALSE)
is.island = nchar(network.vertex.names(finch.nw)) == 1
plot(finch.nw, vertex.cex = 2.5, displaylabels = TRUE, 
     vertex.col = ifelse(is.island, "forestgreen", "gold3"),
     label= sub(" finch", "", network.vertex.names(finch.nw)))
finch.nw |> as.matrix() |> t() |> (\(x) x[, order(colnames(x))])()

## -----------------------------------------------------------------------------
#| label: fig-ggnetworkong1-1
#| fig-width: 3
#| fig-height: 3
#| out-width: 33%
#| fig-cap: "A **[ggraph](https://cran.r-project.org/web/packages/ggraph/)** example."
library("ggraph")
ggraph(g1, layout = "nicely") + 
  geom_edge_link() + 
  geom_node_point(size=6,color="#8856a7") + 
  geom_node_text(label=vertex_attr(g1)$name,  color="white")

## -----------------------------------------------------------------------------
#| label: T1net
#| fig-show: hide
library("markovchain")
statesNames = c("A", "C", "G","T")
T1MC = new("markovchain", states = statesNames, transitionMatrix =
  matrix(c(0.2,0.1,0.4,0.3,0,1,0,0,0.1,0.2,0.2,0.5,0.1,0.1,0.8,0.0),
         nrow = 4,byrow = TRUE, dimnames = list(statesNames, statesNames)))
plot(T1MC, edge.arrow.size = 0.4, vertex.color = "purple",
     edge.arrow.width = 2.2, edge.width = 5, edge.color = "blue",
     edge.curved = TRUE, edge.label.cex = 2.5, vertex.size= 32,
     vertex.label.cex = 3.5, edge.loop.angle = 3,
     vertex.label.family = "sans", vertex.label.color = "white")

## -----------------------------------------------------------------------------
#| label: fig-fourstateMC
#| echo: false
#| column: margin

## ----warning=FALSE, message=FALSE---------------------------------------------
#| label: fig-completechemokine
#| warning.known: !expr c("unlabeled data points")
#| echo: !expr -c(2,9)
#| fig-width: 4
#| fig-height: 3
#| out-width: 75%
#| fig-margin: false
#| fig-cap: "Perturbed chemokine subnetwork uncovered in  @YuGXNA using differential gene expression patterns in sorted T-cells. Notice the clique-like structure of the genes CXCR3, CXCL13, CCL19, CSCR5 and CCR7 in the right hand corner."
oldpar = par(mar = c(5.1,4.7,4.1,2.6))
datf= read.table(file.path(DATA, "string_graph.txt"), header = TRUE)
grs = graph_from_data_frame(datf[, c("node1", "node2")], directed = FALSE)
E(grs)$weight = 1
V(grs)$size = centralization.degree(grs)$res
ggraph(grs) +
  geom_edge_arc(color = "black",  strength = 0.05, alpha = 0.8)+
  geom_node_point(size = 2.5, alpha = 0.5, color = "orange") +
  geom_node_label(aes(label=vertex_attr(grs)$name), size = 3, alpha = 0.9, color = "#8856a7",repel=TRUE) 

par(oldpar)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "A long unstructured laundry list of possibly differentially expressed genes can be daunting."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: testGSEABase
#| eval: false
## library("GSEABase")
## ## This requires a login to the website.
## fl   =  "/path/to/msigdb_v5.1.xml"
## gss  =  getBroadSets(fl)
## organism(gss[[1]])
## table(sapply(gss, organism))

## -----------------------------------------------------------------------------
#| label: tbl-hypergeom
#| echo: false
#| column: margin
#| tbl-cap-location: top
#| tbl-cap: "Although there are the same number of each category of gene found in the *significant* set, both the simulation below and the theory of testing in two-way tables shows us that the blue category is enriched."
d1 <- data.frame(Yellow = c(25, 500), Blue = c(25, 100), Red = c(25, 400))
row.names(d1) <- c('Significant', 'Universe')
knitr::kable(d1)

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "So-called 'exact' tests because they are nonparametric and based on exhaustive enumerations: **not** because we are sure of the answer – this is statistics after all."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: simulationFishertest
universe = c(rep("Yellow", 500), rep("Blue", 100), rep("Red", 400))
countblue = replicate(20000, {
  pick75 = sample(universe, 75, replace = FALSE)
  sum(pick75 == "Blue")
})
summary(countblue)

## -----------------------------------------------------------------------------
#| label: fig-histblue
#| echo: false
#| out-width: 50%
#| fig-width: 6
#| fig-height: 4
#| fig-cap: !expr paste("We can see that even in", length(countblue), "simulations, no blue count comes close to being 25. We can reject such an event as having happened by chance and conclude that the blue are **enriched**.")
ggplot(data.frame(countblue), aes(x = countblue)) +
  geom_histogram(binwidth = 1, colour = "white", fill = "purple", center = 0.5, linewidth = 1)

## -----------------------------------------------------------------------------
#| label: checkcountblue
#| echo: false
stopifnot(all(countblue<22))

## -----------------------------------------------------------------------------
#| label: fig-GOplotEC
#| warning.known: !expr c("using size for a discrete variable is not advised", "rows containing missing values")
#| fig-width: 11
#| fig-height: 12
#| dpi: 80
#| dev: jpeg
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "This graph shows the correspondence between GO terms and significantly changed genes in a study on differential expression in endothelial cells from two steady state tissues (brain and heart, see  @Nolan:2013). After normalization a differential expression analysis was performed giving a list of genes. A gene-annotation enrichment analysis of the set of differentially expressed genes (adjusted p-value < 0.05) was then performed with the **[GOplot](https://cran.r-project.org/web/packages/GOplot/)** package."
library("GOplot")
data("EC")
circ  =  circle_dat(EC$david, EC$genelist)
chord =  chord_dat(circ, EC$genes, EC$process)
GOChord(chord, limit = c(0, 5))

## -----------------------------------------------------------------------------
#| label: BioNetsimple
library("BioNet")
library("DLBCL")
data("dataLym")
data("interactome")
interactome
pval = dataLym$t.pval
names(pval)  =  dataLym$label
subnet = subNetwork(dataLym$label, interactome)
subnet = rmSelfLoops(subnet)
subnet

## -----------------------------------------------------------------------------
#| label: auxiliaryggplotfit
#| include: false
## Function to qqplot the output from fitBumModel
ggplotqqbu = function(x) {
  n = length(x$pvalues)
  probs = (rank(sort(x$pvalues)) - 0.5)/n
  quantiles = unlist(sapply(probs, uniroot, f = BioNet:::.pbum.solve,
        interval = c(0, 1), lambda = x$lambda, a = x$a)[1, ])
  df = data.frame(fitted = quantiles, observed = sort(x$pvalues))
  ggplot(df, aes(x = fitted, y = observed)) +
      xlab("Theoretical quantiles") + ylab("Observed quantiles") +
      geom_point(size = 0.3, alpha = 0.3, color = "red")+
      geom_segment(aes(x = 0, y = 0, xend = 1, yend= 1 ), color = "blue") +
      coord_fixed(ratio = 1)
}

hist1.bum = function(x, breaks = 50){
  n = length(x$pvalues)
  bumdata = seq(from = 0.006, to = 1, length=n)
  ys = x$lambda + (1 - x$lambda) * x$a * bumdata^(x$a -1)
  dat = data.frame(pvalues = x$pvalues, xxs = bumdata, y = ys)
  y0 = BioNet:::piUpper(x)
  ggplot(dat, aes(x = pvalues)) +
    geom_histogram(aes( y= after_stat(density)), binwidth = .02, fill = "orange", alpha = 0.75) +
    geom_hline(yintercept = y0, color = "blue3", alpha = 0.5) +
    geom_line(aes(x = xxs, y = y), col = "red3", alpha = 0.5, linewidth = exp(0.5)) +
    xlab("p-values") +
    annotate("text", x = -0.03, y = y0 + 0.5, label = "pi[0]", parse = TRUE, size = 8)
}

## -----------------------------------------------------------------------------
#| label: originalBioNet
#| include: false
#| eval: false
## ## Original analysis as done by the authors using both p-values.
## library("BioNet")
## library("DLBCL")
## data("dataLym")
## data("interactome")
## interactome
## pvals  =  cbind(t = dataLym$t.pval, s = dataLym$s.pval)
## rownames(pvals)  =  dataLym$label
## pval  =  aggrPvals(pvals, order = 2, plot = FALSE)
## subnet  =  subNetwork(dataLym$label, interactome)
## subnet  =  rmSelfLoops(subnet)
## subnet

## -----------------------------------------------------------------------------
#| label: fitBUM
fb = fitBumModel(pval, plot = FALSE)
fb
scores=scoreNodes(subnet, fb, fdr = 0.001)

## -----------------------------------------------------------------------------
#| label: fig-plotFITBum
#| echo: false
#| column: margin
#| fig-asp: 1
#| fig-height: 4
#| fig-cap: "The qqplot shows the quality of the fit of beta-uniform mixture model to the data. The red points have the theoretical quantiles from the beta distribution as the x coordinates the observed quantiles and the y coordinates. The blue line shows that this model fits nicely."
ggp=ggplotqqbu(fb)
print(ggp)

## -----------------------------------------------------------------------------
#| label: fig-histFITBum
#| echo: false
#| column: margin
#| fig-asp: 1
#| fig-height: 4
#| fig-cap: "A histogram of the mixture components for the p-values, the beta in red and the uniform in blue, $\\pi_0$ is the mixing proportion assigned to the null component whose distribution should be uniform."
ggh = hist1.bum(fb)
print(ggh)

## -----------------------------------------------------------------------------
#| label: subgraphHeinz
hotSub  =  runFastHeinz(subnet, scores)
hotSub
logFC=dataLym$diff
names(logFC)=dataLym$label

## -----------------------------------------------------------------------------
#| label: fig-plotBioNet
#| out-width: 75%
#| fig-width: 8
#| fig-height: 8
#| dpi: 80
#| fig-margin: false
#| fig-cap: "The subgraph found as maximally enriched for differential expression between ABC and GCB B-cell lymphoma. The nodes are colored in red and green: green shows an upregulation in ACB and red an upregulation in GBC. The shape of the nodes depicts the score: rectangles indicate a negative score, circles a positive score."
plotModule(hotSub, layout = layout.davidson.harel, scores = scores,
                  diff.expr = logFC)

## -----------------------------------------------------------------------------
#| label: fig-IntroTree
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: fig-HIVtree
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "This is called the **molecular clock** hypothesis, if we do not make this assumption we run into what is known as non-identifiability (ie we can't tell the difference between the many possible mutational histories given the observed data)."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "We use an error term written here as $o(h)$, we read this little $o$ of $h$, which just means that this error terms grows much slower (i.e., sublinear) than $h$."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: tbl-JCKimura
#| echo: false
#| tbl-cap-location: margin
#| tbl-cap: "Two examples of rate matrices, on the left:the Jukes-Cantor (`JC69`) model, on the right is shown the Kimura (`K80`) two parameter model."
knitr::kable(col.names = c("",""), align = "lr", data.frame(
  "$Q = \\begin{array}{lcccc} & A & T & C & G \\\\ 
  A & -3\\alpha & \\alpha & \\alpha & \\alpha \\\\ 
  T & \\alpha & -3\\alpha & \\alpha & \\alpha \\\\ 
  C & \\alpha & \\alpha & -3\\alpha & \\alpha \\\\ 
  G & \\alpha & \\alpha & \\alpha & -3\\alpha \\\\ \\end{array}$", 
  "$Q = \\begin{array}{lcccc} & A & T & C & G \\\\
  A & -\\alpha-2 \\beta & \\beta & \\beta & \\alpha \\\\ 
  T & \\beta & -\\alpha-2 \\beta & \\alpha & \\beta \\\\
  C & \\beta & \\alpha & -\\alpha-2 \\beta & \\beta \\\\
  G & \\alpha & \\beta & \\beta & -\\alpha-2 \\beta \\\\ \\end{array}$"))

## -----------------------------------------------------------------------------
#| child: devil.qmd
#| fig-cap: "Vocabulary overload here! : Transitions in this context mean mutational changes within the purines (A<->G]) or within the pyrimidines (C <-> T); whereas when we talked about Markov chains earlier our **transition** matrix contains all probabilities of any state changes."

## -----------------------------------------------------------------------------
#| echo: false

## -----------------------------------------------------------------------------
#| label: Atree
library("phangorn")
library("ggtree")
load(file.path(DATA,"tree1.RData"))

## -----------------------------------------------------------------------------
#| label: fig-Atree1b
#| column: margin
#| fig-asp: 1
#| fig-height: 5
#| fig-cap: "This is the tree we use as our *true* parameter. We generate nucleotides one at a time from the root and `dropping' them down the tree. With some probability proportional to the edge lengths, mutations occur down the branches."
ggtree(tree1, lwd = 2, color = "darkgreen", alpha = 0.8, right = TRUE) +
  geom_tiplab(size = 7, angle = 90, offset = 0.05) +
  geom_point(aes(shape = isTip, color = isTip), size = 5, alpha = 0.6)

## -----------------------------------------------------------------------------
#| label: fig-ggtreeAlignment
#| column: margin
#| fig-height: 4.5
#| fig-margin: false
#| fig-cap: "The tree on the left was used to generate the sequences on the right according to a Jukes Cantor model. The nucleotide frequencies generated at the root were quite unequal, with `A` and `C` being generated more rarely. As the sequences percolate down the tree, mutations occur, they are more likely to occur on the longer branches. "
seqs6 = simSeq(tree1, l = 60, type = "DNA", bf = c(1, 1, 3, 3)/8, rate = 0.1)
seqs6
mat6df = data.frame(as.character(seqs6))
p = ggtree(tree1, lwd = 1.2) + geom_tiplab(aes(x = branch), size = 5, vjust = 2)
gheatmap(p, mat6df[, 1:60], offset = 0.01, colnames = FALSE)

## -----------------------------------------------------------------------------
#| label: Magnify
#| echo: false

## -----------------------------------------------------------------------------
#| label: fig-igraphsteiner
#| echo: false
#| column: margin
#| fig-cap: "A Steiner tree, the inner points are represented as squares. The method for creating the shortest tree that passes through all outer 1,2,5,6 is to create two inside (\"ancester\") points 3 and 4."
plot(g1, vertex.size=20, edge.width=5, vertex.color="coral",
     vertex.shape=c("circle","square")[c(1,1,2,2,1,1)])

## -----------------------------------------------------------------------------
#| label: fig-njtree1
#| warning.known: !expr c("tree contained negative edge length")
#| column: margin
#| fig-cap: "Trees built with a neighbor joining algorithm are very fast to compute and are often used as initial values for more expensive estimation procedures such as the maximum likelihood or parsimony."
tree.nj = nj(dist.ml(seqs6, "JC69"))
ggtree(tree.nj) + geom_tiplab(size = 7) 

## -----------------------------------------------------------------------------
#| label: pmltree1
fit = pml(tree1, seqs6, k = 4)

## -----------------------------------------------------------------------------
#| label: readseqs
library("dada2")
seqtab = readRDS(file.path(DATA,"seqtab.rds"))
seqs = getSequences(seqtab)
names(seqs) = seqs

## -----------------------------------------------------------------------------
#| label: tax
#| eval: false
## fastaRef = "../tmp/rdp_train_set_16.fa.gz"
## taxtab = assignTaxonomy(seqtab, refFasta = fastaRef)

## -----------------------------------------------------------------------------
#| label: taxtabload
taxtab = readRDS(file.path(DATA,"taxtab16.rds"))
dim(taxtab)

## -----------------------------------------------------------------------------
#| label: taxtabanswer
head(taxtab) |> `rownames<-`(NULL)

## -----------------------------------------------------------------------------
#| label: readDNAalign
readLines(file.path(DATA,"mal2.dna.txt")) |> head(12) |> cat(sep="\n")

## -----------------------------------------------------------------------------
#| label: alignwithdecipher
library("DECIPHER")
alignment = AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

## -----------------------------------------------------------------------------
#| label: treeFIT
phangAlign = phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm = phangorn::dist.ml(phangAlign)
treeNJ = phangorn::NJ(dm)   # Note: tip order != sequence order
fit = phangorn::pml(treeNJ, data = phangAlign)
fitGTR = update(fit, k = 4, inv = 0.2)
fitGTR = phangorn::optim.pml(fitGTR, model = "GTR", optInv = TRUE,
         optGamma = TRUE,  rearrangement = "stochastic",
         control = phangorn::pml.control(trace = 0))

## -----------------------------------------------------------------------------
#| label: samdatprepare
samples = read.csv(file.path(DATA, "MIMARKS_Data_combined.csv"), header = TRUE)
samples$SampleID = paste0(gsub("00", "", samples$host_subject_id), 
                          "D", samples$age-21) 
samples = samples[!duplicated(samples$SampleID), ] 
stopifnot(all(rownames(seqtab) %in% samples$SampleID))
rownames(samples) = samples$SampleID 
keepCols = c("collection_date", "biome", "target_gene", "target_subfragment", 
  "host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
  "diet", "family_relationship", "genotype", "SampleID") 
samples = samples[rownames(seqtab), keepCols] 

## -----------------------------------------------------------------------------
#| label: samdatcheck
#| echo: false
# For 'sex', samples_2 has oddly the values 'FALSE' 
samples_2 = read.csv(file.path(DATA,"MIMARKS_Data_clean.csv"), header = TRUE)
for (k in setdiff(keepCols, "sex")) stopifnot(identical(samples[[k]], samples_2[[k]]))

## -----------------------------------------------------------------------------
#| label: phyloseqOBJ
library("phyloseq")
pso = phyloseq(tax_table(taxtab), 
               sample_data(samples),
               otu_table(seqtab, taxa_are_rows = FALSE), 
               phy_tree(fitGTR$tree))

## -----------------------------------------------------------------------------
#| label: onelineprune
prune_samples(rowSums(otu_table(pso)) > 5000, pso)

## -----------------------------------------------------------------------------
#| label: what-phyla
prevalence = apply(X = otu_table(pso),
                   MARGIN = ifelse(taxa_are_rows(pso), yes = 1, no = 2),
                   FUN = function(x) {sum(x > 0)})
prevdf = data.frame(Prevalence = prevalence,
                    TotalAbundance = taxa_sums(pso),
                    tax_table(pso))
tab = table(prevdf$Phylum)
keepPhyla = names(tab)[tab>5]
prevdf1   = subset(prevdf,   Phylum %in% keepPhyla)
ps2v      = subset_taxa(pso, Phylum %in% keepPhyla)

## -----------------------------------------------------------------------------
#| label: deseq-transform
#| warning: false
# warning: !expr c("DESeqDataSet.se, design = design, ignoreRank.: some variables in design formula are characters, converting to factors")
library("DESeq2")
ps1 = readRDS(file.path(DATA,"ps1.rds"))
ps_dds = phyloseq_to_deseq2(ps1, design = ~ ageBin + family_relationship)
geometricmean = function(x)
   if (all(x == 0)) { 0 } else { exp(mean(log(x[x != 0]))) }
geoMeans = apply(counts(ps_dds), 1, geometricmean)
ps_dds = estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds = estimateDispersions(ps_dds)
abund = getVarianceStabilizedData(ps_dds)

## -----------------------------------------------------------------------------
#| label: shorten
rownames(abund) = substr(rownames(abund), 1, 5) |> make.names(unique = TRUE)

## -----------------------------------------------------------------------------
#| label: structssi-unadjp
library("structSSI")
el = phy_tree(ps1)$edge
el0 = el
el0 = el0[rev(seq_len(nrow(el))), ]
el_names = c(rownames(abund), seq_len(phy_tree(ps1)$Nnode))
el[, 1] = el_names[el0[, 1]]
el[, 2] = el_names[el0[, 2]]
unadj_p = treePValues(el, abund, sample_data(ps1)$ageBin)

## -----------------------------------------------------------------------------
#| label: structssi-test
#| results: false
hfdr_res = hFDR.adjust(unadj_p, el, 0.75)
summary(hfdr_res)
#plot(hfdr_res, height = 5000) # not run: opens in a browser

## -----------------------------------------------------------------------------
#| label: fig-structssi-hfdr
#| echo: false
#| column: margin

## -----------------------------------------------------------------------------
#| label: structssi-tax
library("dplyr")
options(digits = 3)
tax = tax_table(ps1)[, c("Family", "Genus")] |> data.frame()
tax$seq = rownames(abund)
hfdr_res@p.vals$seq = rownames(hfdr_res@p.vals)
left_join(tax, hfdr_res@p.vals[,-3]) |>
  arrange(adjp) |> head(9) |> dplyr::select(1,2,4,5)

## -----------------------------------------------------------------------------
#| label: ST-setup
#| echo: false
library("igraph")
library("ggraph")
pts = structure(c(0, 0, 1, 1, 1.5, 2, 0, 1, 1, 0, 0.5, 0.5),
                .Dim = c(6L, 2L))
matxy = pts
distxy = stats::dist(matxy)
g = graph_from_adjacency_matrix(as.matrix(distxy), weighted = TRUE)
mst1 = igraph::mst(g)

## ----warning=TRUE-------------------------------------------------------------
#| label: fig-graphs-MST
#| echo: false
#| layout-ncol: 1
#| column: margin
#| fig-height: 2.5
#| fig-cap: "Two spanning trees for the same set of six vertices. The blue graph is the minimum spanning tree, if the Euclidean distances between the points in the 2D plane are used."
#| fig-subcap:
#|   - ""
#|   - ""
gred = igraph::make_ring(6) - igraph::edges(6) 

ggraph(gred,layout = pts) +
  geom_edge_link(color = "red", alpha = 0.8 , 
                 width= 2) +
  geom_node_point(color="black", size=4) +
  theme_graph(background = "white")

ggraph(mst1, layout = pts) +
  geom_edge_link(color = "steelblue", alpha = 0.8, 
                 width= 2) +
  geom_node_point(color="black", size = 4) +
  theme_graph(background = "white")

## -----------------------------------------------------------------------------
#| label: fig-HIVMSTi
#| fig-width: 7.5
#| fig-height: 6
#| dpi: 80
#| fig-margin: false
#| fig-cap: "The minimum spanning tree computed from DNA distances between HIV sequences from samples taken in 2009 and whose country of origin was known, data as published in the `HIVdb` database [@HIVdb]."
load(file.path(DATA, "dist2009c.RData"))
country09 = attr(dist2009c, "Label")
mstree2009 = ape::mst(dist2009c)
gr09 = graph_from_adjacency_matrix(mstree2009, mode = "undirected")
ggraph(gr09, layout="fr") +
  geom_edge_link(color = "black",alpha=0.5) +
  geom_node_point(aes(color = vertex_attr(gr09)$name), size = 2) +
  geom_node_text(aes(label = vertex_attr(gr09)$name), color="black",size=2) +
  theme_void() +
  guides(color=guide_legend(keyheight=0.1,keywidth=0.1,
      title="Countries"))

## -----------------------------------------------------------------------------
#| label: fig-graphs-networklabelrepel
#| warning.known: !expr c("too many overlaps")
#| fig-width: 6
#| fig-height: 6
#| dpi: 80
#| fig-margin: false
#| fig-cap: "Solution to @prp-graphs-networklabelrepel."
library("ggraph")
ggraph(gr09, layout="fr") +
  geom_edge_link(color = "black",alpha=0.5) +
  geom_node_point(aes(color = vertex_attr(gr09)$name), size = 2) +
  geom_node_label(aes(label = vertex_attr(gr09)$name), color="black",size=2,repel=TRUE) +
  theme_void() +
  guides(color=guide_legend(keyheight=0.1,keywidth=0.1,
      title="Countries"))

## -----------------------------------------------------------------------------
#| label: fig-HIVmap
#| fig-width: 8
#| fig-height: 6
#| dpi: 80
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "A minimum spanning tree between HIV cases. The geographic locations of the cases were jittered to reduce overlapping. The DNA sequence distances between the HIV strains were used as input to an undirected minimum spanning tree algorithm."
library("rworldmap")
mat = match(country09, countriesLow$NAME)
coords2009 = data.frame(
  lat = countriesLow$LAT[mat],
  lon = countriesLow$LON[mat],
  country = country09)
layoutCoordinates = cbind(
  x = jitter(coords2009$lon, amount = 15),
  y = jitter(coords2009$lat, amount = 8))
labc = names(table(country09)[which(table(country09) > 1)])
matc = match(labc, countriesLow$NAME)
dfc = data.frame(
  latc = countriesLow$LAT[matc],
  lonc = countriesLow$LON[matc],
  labc)
dfctrans = dfc
dfctrans[, 1] = (dfc[,1] + 31) / 93
dfctrans[, 2] = (dfc[,2] + 105) / 238
Countries = vertex_attr(gr09)$name 

ggraph(gr09, layout=layoutCoordinates) +
  geom_node_point(aes(color=Countries),size = 3, alpha=0.75) +
  geom_edge_arc(color = "black", alpha = 0.5, strength=0.15) +
  geom_label(data=dfc,aes(x=lonc,y=latc,label=labc,fill=labc),colour="white",alpha=0.8,size=3,show.legend=F) +
  theme_void()  

## -----------------------------------------------------------------------------
#| label: fig-WWtest
#| echo: false
#| fig-width: 6
#| fig-height: 1
#| fig-margin: false
#| fig-cap-location: margin
#| fig-cap: "Seeing the number of runs in a one-dimensional, two-sample, nonparametric Wald-Wolfowitz test can indicate whether the two groups have the same distributions."
dfbr=data.frame(measure=c(rnorm(15,0.9),rnorm(15,1.8)),
  group=as.factor(c(rep("men",15),rep("women",15))))
ggplot(dfbr,aes(x=measure,group=group,y=0)) + ylim(-0.25,+0.25) +
  geom_point(aes(col=group,x=measure,y=0,shape=group),size=5,alpha=0.6)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw() + geom_hline(yintercept = 0) +
  theme(panel.border = element_blank(),
  axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank() )+ coord_fixed()

## -----------------------------------------------------------------------------
#| label: preparegraph
ps1  = readRDS(file.path(DATA,"ps1.rds"))
sampledata = data.frame( sample_data(ps1))
d1 = as.matrix(phyloseq::distance(ps1, method="jaccard"))
gr = graph_from_adjacency_matrix(d1,  mode = "undirected", weighted = TRUE)
net = igraph::mst(gr)
V(net)$id = sampledata[names(V(net)), "host_subject_id"]
V(net)$litter = sampledata[names(V(net)), "family_relationship"]

## -----------------------------------------------------------------------------
#| label: fig-mstplot
#| column: margin
#| fig-width: 4
#| fig-height: 4.5
#| fig-cap: "The minimum spanning tree based on Jaccard dissimilarity and annotated with the mice ID and litter factors"

ggraph(net, layout="fr")+
  geom_edge_arc(color = "darkgray") +
  geom_node_point(aes(color = id, shape = litter)) + 
  theme(legend.position="bottom")

## -----------------------------------------------------------------------------
#| label: MSTJaccardplain
library("phyloseqGraphTest")
gt = graph_perm_test(ps1, "host_subject_id", distance="jaccard",
                     type="mst",  nperm=1000)
gt$pval

## -----------------------------------------------------------------------------
#| label: fig-mstJaccard
#| column: margin
#| fig-width: 4
#| fig-height: 2.1
#| fig-cap: "The permutation histogram of the number of pure edges in the network obtained from the minimal spanning tree with Jaccard similarity."
plot_permutations(gt)

## -----------------------------------------------------------------------------
#| label: ggnetworkphyl
net = make_network(ps1, max.dist = 0.35)
sampledata = data.frame(sample_data(ps1))
V(net)$id = sampledata[names(V(net)), "host_subject_id"]
V(net)$litter = sampledata[names(V(net)), "family_relationship"]

## ----warning=FALSE------------------------------------------------------------
#| label: fig-ggnetworkplotJ
#| fig-width: 5
#| fig-height: 4
#| out-width: 75%
#| fig-cap-location: margin
#| fig-cap: "A co-occurrence network created by using a threshold on the Jaccard dissimilarity matrix. The colors represent which mouse the sample came from; the shape represents which litter the mouse was in."
ggraph(net, layout="fr") +
  geom_edge_link(color = "darkgray") +
  geom_node_point(aes(color = id, shape = litter)) + 
    theme(plot.margin = unit(c(0, 5, 2, 0), "cm"))+
    theme(legend.position = c(1.4, 0.3),legend.background = element_blank(),
          legend.margin=margin(0, 3, 0, 0, "cm"))+
         guides(color=guide_legend(ncol=2))+
  theme_graph(background = "white")

## -----------------------------------------------------------------------------
#| label: mst
gt = graph_perm_test(ps1, "family_relationship",
        grouping = "host_subject_id",
        distance = "jaccard", type = "mst", nperm= 1000)
gt$pval

## -----------------------------------------------------------------------------
#| label: fig-mstpermplotNest
#| column: margin
#| fig-width: 4
#| fig-height: 2.5
#| fig-cap: "The permutation histogram obtained from the minimal spanning tree with Jaccard similarity."
plot_permutations(gt)

## -----------------------------------------------------------------------------
#| label: knn1test
gtnn1 = graph_perm_test(ps1, "family_relationship",
                      grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
gtnn1$pval

## -----------------------------------------------------------------------------
#| label: fig-knn-1-plot
#| fig-width: 4
#| fig-height: 3
#| out-width: 75%
#| fig-cap: "The graph obtained from a nearest-neighbor graph with Jaccard similarity."
plot_test_network(gtnn1)

## -----------------------------------------------------------------------------
#| label: adjacencyplot2
#| eval: false
#| ref-label: 'adjacencyplot1'
## NA

## -----------------------------------------------------------------------------
#| label: adjmatrix2
#| eval: false
#| ref-label: 'adjmatrix1'
## NA

## -----------------------------------------------------------------------------
#| label: fig-PSB-MC-s
#| echo: false
#| out-width: 75%

## -----------------------------------------------------------------------------
#| label: plot-MC-plus
#| eval: false
## library("markovchain")
## # Make Markov chain object
## mcPreg  =  new("markovchain", states = CSTs,
##               transitionMatrix = trans, name="PregCST")
## mcPreg
## # Set up igraph of the markov chain
## netMC  =  markovchain:::.getNet(mcPreg, round = TRUE)

## -----------------------------------------------------------------------------
#| label: CSTMarkov
#| eval: false
## wts  =  E(netMC)$weight/100
## edgel  =  get.edgelist(netMC)
## elcat  =  paste(edgel[,1], edgel[,2])
## elrev  =  paste(edgel[,2], edgel[,1])
## edge.curved  =  sapply(elcat, function(x) x %in% elrev)
## samples_def  =  data.frame(sample_data(ps))
## samples_def  =  samples_def[samples$Preterm | samples$Term,] # Only those definitely assigned, i.e. not marginal
## premat  =  table(samples_def$CST, samples_def$Preterm)
## rownames(premat)  =  markovchain::states(mcPreg)
## colnames(premat)  =  c("Term", "Preterm")
## premat
## premat  =  premat/rowSums(premat)
## vert.CSTclrs  =  CSTColors

## -----------------------------------------------------------------------------
#| label: CSTMC
#| eval: false
## default.par  =  par(no.readonly = TRUE)
## # Define color scale
## # Plotting function for markov chain
## plotMC  =  function(object, ...) {
##     netMC  =  markovchain:::.getNet(object, round = TRUE)
##     plot.igraph(x = netMC, ...)
## }
## # Color bar for the markov chain visualization, gradient in strength of preterm association
## color.bar  =  function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title=NULL) {
##     scale = (length(lut)-1)/(max-min)
##     cur.par = par(no.readonly = TRUE)
##     par(mar = c(0, 4, 1, 4) + 0.1, oma = c(0, 0, 0, 0) + 0.1)
##     par(ps = 10, cex = 0.8)
##     par(tcl=-0.2, cex.axis=0.8, cex.lab = 0.8)
##     plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab=", yaxt='n', ylab=", main=title)
##     axis(1, c(0, 0.5, 1))
##     for (i in 1:(length(lut)-1)) {
##       x = (i-1)/scale + min
##       rect(x,0,x+1/scale,10, col=lut[i], border=NA)
##     }
## }
## 
## pal  =  colorRampPalette(c("grey50", "maroon", "magenta2"))(101)
## vert.clrs  =  sapply(states(mcPreg), function(x) pal[1+round(100*premat[x,"Preterm"])])
## vert.sz  =  4 + 2*sapply(states(mcPreg),
##               function(x) nrow(unique(sample_data(ps)[sample_data(ps)$CST==x,"SubjectID"])))
## vert.sz  =  vert.sz * 0.85
## vert.font.clrs  =  c("white", "white", "white", "white", "white")
## 
## # E(netMC) to see edge list, have to define loop angles individually by the # in edge list, not vertex
## edge.loop.angle = c(0, 0, 0, 0, 3.14, 3.14, 0, 0, 0, 0, 3.14, 0, 0, 0, 0, 0)-0.45
## layout  =  matrix(c(0.6,0.95, 0.43,1, 0.3,0.66, 0.55,0.3, 0.75,0.65), nrow = 5, ncol = 2, byrow = TRUE)
## 
## # Color by association with preterm birth
## layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,10))
## color.bar(pal, min=0, max=1, nticks=6, title="Fraction preterm")
## par(mar=c(0,1,1,1)+0.1)
## edge.arrow.size=0.8
## edge.arrow.width=1.4
## edge.width = (15*wts + 0.1)*0.6
## edge.labels  =  as.character(E(netMC)$weight/100)
## edge.labels[edge.labels<0.4]  =  NA  # labels only for self-loops
## plotMC(mcPreg, edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
##        edge.width=edge.width, edge.curved=edge.curved,
##        vertex.color=vert.clrs, vertex.size=(vert.sz),
##        vertex.label.font = 2, vertex.label.cex = 1,
##        vertex.label.color = vert.font.clrs, vertex.frame.color = NA,
##        layout=layout, edge.loop.angle = edge.loop.angle)
## par(default.par)

## -----------------------------------------------------------------------------
#| label: fig-ccnb1img
#| echo: false
#| out-width: 75%

## -----------------------------------------------------------------------------
#| label: ccnb1data
dat = read.table(file.path(DATA,"ccnb1datsmall.txt"), header = TRUE, comment.char = "", stringsAsFactors = TRUE)
v = levels(unlist(dat[,1:2]))        # vertex names
n = length(v)                        # number of vertices
e = matrix(match(as.character(unlist(dat[,1:2])), v),ncol=2) # edge list
w = dat$coexpression                 # edge weights

## -----------------------------------------------------------------------------
#| label: ccnbet
M = matrix(0, n, n)
M[e] = w
M = M + t(M)
dimnames(M) = list(v, v)
A = 1*(M > 0)

## -----------------------------------------------------------------------------
#| label: plotgraph
#| eval: false
## library(igraph)
## net = network(e, directed=FALSE)
## par(mar=rep(0,4))
## plot(net, label=v)

## -----------------------------------------------------------------------------
#| label: fig-heatmapCCNB1
#| fig-width: 6
#| fig-cap: "This represents the adjacency of the CCNB1 network -- 2 step neighborhood with co-expression levels $\\geq$ 0.900, generated from R (darker is closer to 1, we ignore values < 0.9)."
breaks  =  c(0, seq(0.9, 1, length=11))
cols  =  grey(1-c(0,seq(0.5,1,length=10)))
ccnb1ind  =  which(v == "CCNB1")
vcols  =  rep("white",n)
vcols[ccnb1ind]  =  "blue"
vcols[which(M[,ccnb1ind]>0 | M[ccnb1ind,])]  =  "red"
par(mar = rep(0, 4))
heatmap(M, symm = TRUE, ColSideColors = vcols, RowSideColors = vcols,
        col = cols, breaks = breaks,  frame = TRUE)
legend("topleft", c("Neighbors(CCNB1)", "CCNB1"),
       fill = c("red","blue"),
       bty = "n", inset = 0, xpd = TRUE,  border = FALSE)

## -----------------------------------------------------------------------------
#| label: readDNA
#| eval: false
## library("ape")
## library("phangorn")
## GAG=read.dna(file.path(DATA,"DNA_GAG_20.txt"))

## -----------------------------------------------------------------------------
#| label: knn-2-plot
gt = graph_perm_test(ps1, "family_relationship", distance = "bray", 
                     grouping = "host_subject_id", type = "knn", knn = 2)
gt$pval

## -----------------------------------------------------------------------------
#| label: fig-knn-2-plot
#| fig-width: 4.5
#| fig-height: 3
#| layout-nrow: 1
#| fig-cap: "The graph (a) and permutation histogram (b) obtained from a two nearest-neighbor graph with Jaccard similarity."
#| fig-subcap:
#|   - "$\\text{}$"
#|   - "$\\text{}$"
plot_test_network(gt)
permdf = data.frame(perm=gt$perm)
obs = gt$observed
ymax = max(gt$perm)
ggplot(permdf, aes(x = perm)) + geom_histogram(bins = 20) +
  geom_segment(aes(x = obs, y = 0, xend = obs, yend = ymax/10), color = "red") +
  geom_point(aes(x = obs, y = ymax/10), color = "red") + xlab("Number of pure edges")
