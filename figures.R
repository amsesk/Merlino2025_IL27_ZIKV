#Interleukin-27 is antiviral at the maternal-fetal interface 
#MM (2025)

# %% load packages----
library(rhdf5) #handling hdf5 file formats
library(tidyverse)
library(tximport) #load Kallisto results
library(EnsDb.Hsapiens.v86)
library(edgeR) #for DGEList object and normalization 
library(limma)
library(gt)
library(DT)
library(RColorBrewer)
library(gplots)
#

# %% Change dir to repo and create output directories, if needed
REPO="./"
setwd(REPO)
output_folders = c("figures", "tables", "rds", "data")
purrr::map(output_folders, \(o) {
             if (!dir.exists(o)) {
               dir.create(o)
             }
})

# %% Import + data cleaning----
#read in your study design
targets <- read_tsv("metadata/studydesign.txt")
sampleLabels <- targets$sample

# %% set file paths to your mapped data
abunds_paths <- glue::glue("data/{targets$mapped}_abundance.tsv")

# %% get annotations 
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) 
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# %% import Kallisto transcript counts
Txi_gene <- tximport(abunds_paths,
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #transcript level for TRUE and gene level for FALSE
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
# %% Filter and Normalize data
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts

# %% Make a DGElist from counts and get CPM
myDGEList <- DGEList(myCounts)
# colnames(myDGEList) = sampleLabels
cpm <- cpm(myDGEList) 

# %%filter
keepers <- rowSums(cpm>1)>=3 #Keeps genes with greater than 1 CPM in at least 3 samples
myDGEList.filtered <- myDGEList[keepers,]

# %% normalize
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM") 

log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
sampleLabels <- targets$sample
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
write_tsv(log2.cpm.filtered.norm.df,"tables/Log2CPM_filter1CPM.txt")

saveRDS(myDGEList.filtered.norm, "rds/Merlino2025_DGEList.filtered.norm.rds")

# %% cd to repo and make output folders
myDGEList.filtered.norm = readRDS("rds/Merlino2025_DGEList.filtered.norm.rds")
targets <- read_tsv("metadata/studydesign.txt")
sampleLabels <- targets$sample

# %% contrast matrix----
condition <- factor(targets$condition, levels = c("Vehicle", "Isotype", "AntiIFNL", "AntiIL27", "IFNLam", "IL27"))
experiment <- factor(targets$experiment, levels = c("REP3", "REP1", "REP2"))

condition
experiment

design <- model.matrix(~0 + condition + experiment) #correct for batch effect
rownames(design) <- sampleLabels

design
rownames(design) == colnames(myDGEList.filtered.norm)

design
myDGEList.filtered.norm
sampleLabels
colnames(myDGEList.filtered.norm) = sampleLabels

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE) #model mean-variance relationship

# %% fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design) #linear model fit


contrast.matrix2 <- makeContrasts(
  pIL27aIL27 = conditionIL27 - conditionAntiIL27, #coef=1
  pIFNLaIFNL = conditionIFNLam - conditionAntiIFNL, #coef=2
  levels=design)

fits2 <- contrasts.fit(fit, contrast.matrix2) #linear model fit
ebFit2 <- eBayes(fits2) #bayesian stats



# %% Heatmap with Anti v Plus comparisons ----

myheatcolors <- colorRampPalette(colors=c("blue","white","red"))(100)

results <- decideTests(ebFit2, method="separate", adjust.method="BH", p.value=0.05, lfc=0.5)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[, 2] != 0 ,]
dim(diffGenes)

# %%
pdf(file = "figures/HeatmapComparingStimsToAntis_test.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 10)
heatmap.2(diffGenes, 
          col=myheatcolors, 
          scale='row', 
          #labRow=rownames(diffGenes),
          labRow = FALSE,
          density.info="none", 
          trace="none",  
          cexRow=1, 
          cexCol=1, 
          margins=c(6,8)) 
dev.off()

# %%

diffGenes_matrix.df <- as_tibble(diffGenes, rownames = "geneID")
write_tsv(diffGenes_matrix.df, "tables/DEGs_StimsToAntis.txt")

# %% Volcano Plots: Plus-Anti----
pIL27aIL27.df <- topTable(ebFit2, adjust ="BH", coef=1, number=40000, sort.by="logFC") %>%
  as_tibble(rownames = "geneID")
pIFNLaIFNL.df <- topTable(ebFit2, adjust ="BH", coef=2, number=40000, sort.by="logFC") %>%
  as_tibble(rownames = "geneID")


shared_genes <- c("IFIT1", "PARP9", "DTX3L", "PARP14", "PSME1")

dots_pIL27aIL27 <- pIL27aIL27.df %>%
  dplyr::filter(geneID %in% shared_genes) 
highlighted_labels_pIL27aIL27 <- dots_pIL27aIL27$geneID

dots_pIFNLaIFNL <- pIFNLaIFNL.df %>%
  dplyr::filter(geneID %in% shared_genes) 
highlighted_labels_pIFNLaIFNL <- dots_pIFNLaIFNL$geneID

# %%
vplot.pIL27aIL27 <- ggplot(pIL27aIL27.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2,
             alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = .5, linetype="longdash", colour="#bf0000", linewidth=1) +
  geom_vline(xintercept = -.5, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 0.5, xmax = 7.5, ymin = -log10(0.05), ymax = 4.5, alpha=.2, fill="#bf0000") +
  annotate("rect", xmin = -0.5, xmax = -7.5, ymin = -log10(0.05), ymax = 4.5, alpha=.2, fill="#2C467A") +
  xlim(-7.5, 7.5) +
  ylim(0,4.5) +
  theme_bw() +
  geom_point(data = dots_pIL27aIL27,
             color = "#bf0000",
             size=3) +
  geom_text(data = dots_pIL27aIL27,
            size = 5,
            fontface = "bold",
            aes(label = highlighted_labels_pIL27aIL27),
            nudge_y = 0.12,
            color = "#bf0000") +
  theme(legend.position = "none")

# %%
vplot.pIFNLaIFNL <- ggplot(pIFNLaIFNL.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2,
             alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = .5, linetype="longdash", colour="#bf0000", linewidth=1) +
  geom_vline(xintercept = -.5, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 0.5, xmax = 7.5, ymin = -log10(0.05), ymax = 4.5, alpha=.2, fill="#bf0000") +
  annotate("rect", xmin = -0.5, xmax = -7.5, ymin = -log10(0.05), ymax = 4.5, alpha=.2, fill="#2C467A") +
  xlim(-7.5, 7.5) +
  ylim(0,4.5) +
  #labs(title="plusIFNL - antiIFNL") +
  theme_bw() +
  geom_point(data = dots_pIFNLaIFNL,
             color = "#bf0000",
             size=3) +
  geom_text(data = dots_pIFNLaIFNL,
            size = 5,
            fontface = "bold",
            aes(label = highlighted_labels_pIFNLaIFNL),
            nudge_y = 0.12,
            color = "#bf0000") +
  theme(legend.position = "none")

# %%
pdf(file = "figures/IL27_volcano.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
vplot.pIL27aIL27
dev.off()

pdf(file = "figures/IFNL_volcano.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
vplot.pIFNLaIFNL
dev.off()

# %% heatmap with more genes----
library(ComplexHeatmap)
library(RColorBrewer) 
library(circlize)
to_specific <- c("CGA", 
                 "ERVW-1", 
                 "TP63", 
                 "HSD3B1", 
                 "XAGE2", 
                 "GATA3")
do_specific <- c("MUC1",
                 "MUC5A",
                 "MUC5B",
                 "SOX17")

il27_specific <- c("IL27",
                   "IL27RA",
                   "EBI3", 
                   "IL6ST")
ifnl_specific <- c("IFNL1",
                   "IFNL2",
                   "IFNL3",
                   "IL10RB", 
                   "IFNLR1") 

comb_list <- c(to_specific, do_specific, ifnl_specific, il27_specific)
comb <- factor(comb_list,
               levels = comb_list)

gene_groups <- factor(rep(c("TO",  "DO","IFNL","IL27"), 
                          times = c(length(to_specific), length(do_specific), length(ifnl_specific), length(il27_specific))),
                      levels = c("TO",  "DO", "IFNL", "IL27"))

colnames(log2.cpm.filtered.norm) <- targets$sample

mat <- log2.cpm.filtered.norm[match(comb, rownames(log2.cpm.filtered.norm)), ]
rownames(mat) <- comb


# %% Define colors and breaks for the heatmap
myCol <- colorRampPalette(c('white', '#bf0000'))(50)
myBreaks <- seq(-5, 15, length.out = 50)

hm <- Heatmap(mat, 
              name = "Log2CPM",
              col = colorRamp2(myBreaks, myCol),
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              cluster_column_slices = FALSE,
              show_column_names = TRUE,               
              column_names_gp = gpar(fontsize = 10),   
              column_gap = unit(1, "mm"),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 10),
              column_title_rot = 0,
              use_raster = FALSE,
              row_split = gene_groups
) 

png("figures/hm_ifnl.png", width=6.5, height=7, unit="in", res = 800)
draw(hm, padding = unit(c(1, 5, 5, 1), "mm")) 
dev.off()

# %% GSEA----
library(tidyverse)
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots
library(msigdbr) # access to msigdb collections directly within R

# use the msigdb package to access up-to-date collections
#msigdbr_species()
hs_gsea <- msigdbr(species = "Homo sapiens") #gets all collections/signatures with human gene IDs

hs_gsea_c2 <- msigdbr(species = "Homo sapiens", 
                      collection = "C2",
                      subcollection = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols 

hs_gsea_hall <- msigdbr(species = "Homo sapiens", 
                        collection = "H") %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols 

hs_gsea_gobp <- msigdbr(species = "Homo sapiens", 
                        collection = "C5",
                        subcollection = "GO:BP") %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols 


# %% prepare data and run gsea for il27----
gseadata_pIL27aIL27.df <- mutate(log2.cpm.filtered.norm.df,
                                 pIL27.AVG = (IL27_1 + IL27_2 + IL27_3)/3, 
                                 aIL27.AVG = (AntiIL27_1 + AntiIL27_2 + AntiIL27_3)/3,
                                 LogFC = (pIL27.AVG - aIL27.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

gseadata_pIL27aIL27_sub.df <- mutate(log2.cpm.filtered.norm.df,
                                     pIL27.AVG = (IL27_1 + IL27_2 + IL27_3)/3, 
                                     aIL27.AVG = (AntiIL27_1 + AntiIL27_2 + AntiIL27_3)/3,
                                     LogFC = (pIL27.AVG - aIL27.AVG)) %>% 
  filter(LogFC >= 0.5) %>%
  mutate_if(is.numeric, round, 2)

# Pull gene symbols and LogFC for the enrichment analysis
gseadata_pIL27aIL27.df_sub <- dplyr::select(gseadata_pIL27aIL27_sub.df, geneID, LogFC)
gseadata_pIL27aIL27.gsea <- gseadata_pIL27aIL27.df_sub$LogFC
names(gseadata_pIL27aIL27.gsea) <- as.character(gseadata_pIL27aIL27.df_sub$geneID)
gseadata_pIL27aIL27.gsea <- sort(gseadata_pIL27aIL27.gsea, decreasing = TRUE)

# %% run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) 
GSEA_pIL27aIL27.res <- GSEA(gseadata_pIL27aIL27.gsea, TERM2GENE=hs_gsea_gobp, verbose=FALSE)
GSEA_pIL27aIL27.df <- as_tibble(GSEA_pIL27aIL27.res@result)

#prepare data and run gsea for ifnL----
gseadata_pIFNLaIFNL.df <- mutate(log2.cpm.filtered.norm.df,
                                 pIFNL.AVG = (IFNLam1 + IFNLam2 + IFNLam3)/3, 
                                 aIFNL.AVG = (AntiIFNL_1 + AntiIFNL_2 + AntiIFNL_3)/3,
                                 LogFC = (pIFNL.AVG - aIFNL.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

gseadata_pIFNLaIFNL_sub.df <- mutate(log2.cpm.filtered.norm.df,
                                     pIFNL.AVG = (IFNLam1 + IFNLam2 + IFNLam3)/3, 
                                     aIFNL.AVG = (AntiIFNL_1 + AntiIFNL_2 + AntiIFNL_3)/3,
                                     LogFC = (pIFNL.AVG - aIFNL.AVG)) %>% 
  filter(LogFC >= 0.5) %>%
  mutate_if(is.numeric, round, 2)

# %% Pull gene symbols and LogFC for the enrichment analysis
gseadata_pIFNLaIFNL.df_sub <- dplyr::select(gseadata_pIFNLaIFNL_sub.df, geneID, LogFC)
gseadata_pIFNLaIFNL.gsea <- gseadata_pIFNLaIFNL.df_sub$LogFC
names(gseadata_pIFNLaIFNL.gsea) <- as.character(gseadata_pIFNLaIFNL.df_sub$geneID)
gseadata_pIFNLaIFNL.gsea <- sort(gseadata_pIFNLaIFNL.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) 
GSEA_pIFNLaIFNL.res <- GSEA(gseadata_pIFNLaIFNL.gsea, TERM2GENE=hs_gsea_gobp, verbose=FALSE)
GSEA_pIFNLaIFNL.df <- as_tibble(GSEA_pIFNLaIFNL.res@result)

write.csv(GSEA_pIFNLaIFNL.df, file = "tables/GSEA_pIFNLaIFNL_sub.csv")
# %%

