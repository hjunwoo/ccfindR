## ------------------------------------------------------------------------
library(ccfindR)

## ---- eval=T,echo=T------------------------------------------------------
# A toy matrix for count data
set.seed(1)
mat <- matrix(rpois(n = 80, lambda = 2), nrow = 4, ncol = 20)
ABC <- LETTERS[1:4]
abc <- letters[1:20]

## ---- eval=F,echo=T------------------------------------------------------
#  # create scNMFSet object
#  sc <- scNMFSet(count = mat)

## ---- eval=T-------------------------------------------------------------
# read count and annotations
genes <- data.frame(ABC)
rownames(genes) <- ABC
cells <- data.frame(abc)
rownames(cells) <- abc
sc <- scNMFSet(count=mat,genes=genes,cells=cells)
sc

## ---- eval=T-------------------------------------------------------------
# read sparse matrix
mat <- Matrix::readMM('pbmc/matrix.mtx')
rownames(mat) <- 1:nrow(mat)
colnames(mat) <- 1:ncol(mat)
sc <- scNMFSet(count=mat)
sc

## ---- eval=T-------------------------------------------------------------
# read 10x files
sc <- read_10x(dir = 'pbmc', count = 'matrix.mtx', 
               genes = 'genes.tsv', barcodes = 'barcodes.tsv')
sc

## ------------------------------------------------------------------------
# slots and subsetting
count(sc)[1:7,1:3]
count(sc)[7,2] <- 0
count(sc)[1:7,1:3]
count(sc)[7,2] <- 2
head(genes(sc))
head(cells(sc))

sc2 <- sc[1:20,1:70]  # subsetting of object
genes(sc2)

## ---- cache=T, fig.width=4, fig.height=4,fig.cap='Quality control filtering of cells. Histogram of UMI counts is shown. Cells can be selected (red) by setting lower and upper thresholds of the UMI count.'----
sc <- filter_cells(sc, umi.min = 10^2.8, umi.max = 10^3.6)

## ---- cache=T, fig.width=4, fig.height=4, fig.cap='Selection of genes for clustering. The scatter plot shows distributions of expression variance to mean ratio (VMR) and the number of cells expressed. Minimum VMR and a range of cell number can be set to select genes (red). Symbols in orange are marker genes provided as input, selected irrespective of expression variance.'----
markers <- c('CD4','CD8A','CD8B','CD19','CD3G','CD3D',
             'CD3Z','CD14')
sc2 <- filter_genes(sc, markers = markers, vmr.min = 2, 
            min.cells.expressed = 100, rescue.genes = FALSE)

## ---- cache=T, eval=T,echo=T, fig.width=4, fig.height=4, fig.cap='Additional selection of genes with modes at nonzero counts. Symbols in blue represent genes rescued.'----
sc_rescue <- filter_genes(sc, markers = markers, vmr.min = 2, min.cells.expressed = 100,
                          rescue.genes = TRUE, progress.bar = FALSE)

## ---- cache=T------------------------------------------------------------
rownames(genes(sc_rescue)) <- rownames(count(sc_rescue)) <- genes(sc_rescue)[,2]
sc <- sc_rescue

## ---- cache=T------------------------------------------------------------
set.seed(2)
sc <- factorize(sc, ranks = 3, nrun = 1, ncnn.step = 1, 
                criterion='connectivity', verbose = 3)

## ---- cache=T------------------------------------------------------------
sc <- factorize(sc, ranks = 3, nrun = 10, verbose = 2)

## ---- cache=T------------------------------------------------------------
sc <- factorize(sc, ranks = 3:8, nrun = 5, verbose = 1, progress.bar = FALSE)

## ------------------------------------------------------------------------
measure(sc)

## ---- fig.width=6.5, fig.height=3, fig.cap='Factorization quality measures as functions of the rank. Residual is (negative) the likelihood function being optimized. Dispersion measures the degree of bimodality in consistency matrix. Cophenetic correlation measures the degree of agreement between consistency matrix and hierarchical clustering.'----
plot(sc)

## ---- cache=T------------------------------------------------------------
set.seed(1)
sb <- sc_rescue
sb <- vb_factorize(sb, verbose = 3, Tol = 2e-4, hyper.update.n0 = 5)

## ---- cache=T------------------------------------------------------------
sb <- vb_factorize(sb, ranks = 2:8, verbose = 1, nrun = 5, progress.bar = FALSE)

## ------------------------------------------------------------------------
head(measure(sb))

## ---- fig.width=4, fig.height=4, fig.cap='Dependence of log evidence with rank.'----
plot(sb)

## ------------------------------------------------------------------------
ranks(sb)

head(basis(sb)[[which(ranks(sb)==5)]]) # basis matrix W for rank 5

## ---- cache=T, fig.width = 4, fig.height=6, fig.cap='Heatmap of basis matrix elements. Marker genes selected in rows, other than those provided as input, are based on the degree to which each features strongly in a particular cluster only and not in the rest. Columns represent the clusters.'----
gene_map(sb, markers = markers, rank = 5, max.per.cluster = 4, gene.name=genes(sb)[,2],
         cexRow = 0.7)

## ------------------------------------------------------------------------
cell_type <- c('CD8+_T','B_cells','CD4+_T','NK','Monocytes')
colnames(basis(sb)[[which(ranks(sb) == 5)]]) <- cell_type
rownames(coeff(sb)[[which(ranks(sb) == 5)]]) <- cell_type

## ---- fig.width = 4, fig.height=4, fig.cap = 'Heatmap of cluster coefficient matrix elements. Rows indicate clusters and columns the cells.'----
cell_map(sb, rank = 5)

## ---- fig.width = 6.5, fig.height = 4, fig.cap='tSNE-based visualization of coefficient matrix elements of cells with colors indicating predicted cluster assignment. The bar plot shows the cell counts of each cluster.'----
visualize_clusters(sb, rank = 5, cex = 0.7)

## ---- fig.width = 4, fig.height = 4, fig.cap = 'Hierarchical tree of clusters derived from varying ranks. The rank increases from 2 to 5 horizontally and nodes are labeled by cluster IDs which bifurcated in each rank.'----
tree <- build_tree(sb, rmax = 5)
tree <- rename_tips(tree, rank = 5, tip.labels = cell_type)
plot_tree(tree, cex = 0.8, show.node.label = TRUE)

