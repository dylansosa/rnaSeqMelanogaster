library(tidyverse)
library(gridExtra)
library(tximport)
library(tximportData)
library(GenomicFeatures)
library(DESeq2)
library(vsn)
library(DEGreport)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(ashr)
library(apeglm)
library(IHW)
library(vsn)
library(LaCroixColoR)

setwd("//Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Yakuba  Long Lab/thesis/1_compare_rnaExpression/rnaseq_brianOliver_forComparison/")
setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/YakubaUChicago-1/thesis/1_compare_rnaExpression/rnaseq_brianOliver_forComparison')

# something with mountain duck changing the wd somtimes (adding '-1')

# interesting: https://www.hadriengourle.com/tutorials/rna/

# import transcript abundances from Salmon
txdb <- makeTxDbFromGFF('annotation/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff')
keytypes(txdb)
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID") 
head(txdf)

# sampleInfo <- read.table('2_8spp_GSE99574/samples_melanogaster.txt')
sampleInfo <- read.table('2_8spp_GSE99574/samples_melanogaster2.txt')
# excluding whole_body, so my numbers are slightly lower than what I commented out.
head(sampleInfo)

colnames(sampleInfo) <- c('runs','sex','tissue')
nrow(sampleInfo)
head(sampleInfo)

sampleInfo$sex <- factor(sampleInfo$sex)
sampleInfo$tissue <- factor(sampleInfo$tissue)
# to compare sex and tissue!
# sampleInfo$sex2 <- factor(sampleInfo$sex, levels = c("male","female"))
# identical so i don't need to use it?

dataDir <- '2_8spp_GSE99574/salmon'
salmonOutput <- file.path(dataDir,sampleInfo$runs, 'quant.sf')
salmonOutput
names(salmonOutput) <- sampleInfo$runs
salmonOutput

tx2gene <- txdf[,2:1]
tx2gene

dmel.txi <- tximport(salmonOutput, type="salmon", tx2gene=tx2gene,
                       countsFromAbundance = "lengthScaledTPM")
# lengthScaledTPM first multiplies TPM by feature length and then scales up to library size. 
# These are then quantities that are on the same scale as original counts, except no longer correlated with feature length across samples.
# get the gene-level counts 
# or additionally scaled using the average transcript length over samples and the library size (lengthScaledTPM). if using scaledTPM or lengthScaledTPM, 
# then the counts are no longer correlated with average transcript length, and so the length offset matrix should not be used.
# dmel.txi$length
dmel.txi$abundance

dmel.dds <- DESeqDataSetFromTximport(dmel.txi, sampleInfo, ~ sex + tissue)
# compare tissue expression for each sex
keep <- rowSums(counts(dmel.dds)) >= 10
dmel.dds <- dmel.dds[keep,]
# pre-filter
# assay(dmel.dds)
dmel.de <-  DESeq(dmel.dds)
# run the three steps of deseq2:

resultsNames(dmel.de)
# [2] "sex_male_vs_female" 

# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta
plotDispEsts(dmel.de)

a <- 0.01
contrast <- c('sex','male','female')
# contrast2 <- c('sex','female','male')
# just causes the values to be reverse, i.e. the LFC >0 for male becomes LFC < 0 when using female as comparison
dmel.res <- results(dmel.de, alpha= a, contrast = contrast)
# get LFC for male vs female sampleInfo in all tissue
dmel.res                                                                                       

dmel.res.lfc.ashr <- lfcShrink(dmel.de, coef = 2, type='ashr')
dmel.res.lfc.apeglm <- lfcShrink(dmel.de, coef = 2, type='apeglm')
dmel.res.lfc.norm <- lfcShrink(dmel.de, coef = 2, type='normal')

dmel.resOrdered <- dmel.res[order(dmel.res$pvalue),]
summary(dmel.res)
# LFC > 0 (up)       : 5616, 34%
# LFC < 0 (down)     : 1133, 7%
# up in male, down in female
# these are not transformed ?
# are result so must be transformed ?

  # post removing whole_body:
# out of 16249 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 4699, 29%
# LFC < 0 (down)     : 688, 4.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 1261, 7.8%

  # if I flip contrast to be female vs male
# out of 16249 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 688, 4.2%
# LFC < 0 (down)     : 4699, 29%
# outliers [1]       : 0, 0%
# low counts [2]     : 1261, 7.8%

sum(dmel.res$padj < 0.01, na.rm=TRUE)
# [1] 6749
# [1] 5387 afrer removing whole_body

###########
resIHW <- results(dmel.de, filterFun=ihw, alpha = a)
summary(resIHW)
# out of 16294 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 5857, 36%
# LFC < 0 (down)     : 1189, 7.3%
sum(resIHW$padj < 0.01, na.rm=TRUE)
# [1] 7046

metadata(resIHW)$ihwResult
###########
par(mfrow=c(1,4), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)

plotMA(dmel.res)
plotMA(dmel.res.lfc.norm,ylim=c(-2,2))
plotMA(dmel.res.lfc.apeglm,ylim=c(-2,2))
plotMA(dmel.res.lfc.ashr,ylim=c(-2,2))

par(mfrow=c(1,1))
# clear canvas

dmel.vsd <- vst(dmel.de,blind = F)
head(assay(dmel.vsd),3)
# transforming raw counts for visualization of PCA
# raw counts are for testing DE

dmel.ntd <- normTransform(dmel.de)
meanSdPlot(assay(dmel.ntd))
meanSdPlot(assay(dmel.vsd))
# may not be possible to smooth variance
# if the samples are truly super different

# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
##################
df <- bind_rows(
  as_data_frame(log2(counts(dmel.de, normalized = TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  # as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(dmel.vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)  

select <- order(rowMeans(counts(dmel.de,normalized=TRUE)),
                decreasing=TRUE)[1:20] 

df <- as.data.frame(colData(dmel.de)[,c("sex","tissue")])

heatmap <- pheatmap(dmel.vsd.ordered[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T, annotation_col=df)
heatmap

pcaData <- plotPCA(dmel.vsd, intgroup=c("tissue", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

Dmel_sex_tissue_expressionPCA <- ggplot(pcaData, aes(PC1, PC2, color=tissue, shape=sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_shape_manual(values = c(0, 19)) +
  ggtitle("D. mel expression") +
  coord_fixed()

ggsave(Dmel_sex_tissue_expressionPCA,f='figures/Dmel_sex_tissue_expressionPCA.png')
##################
######### end deseq vignette and workflow

######## startt 
# https://github.com/hbctraining/DGE_workshop/blob/master/lessons/05_DGE_DESeq2_analysis2.md

### Set thresholds
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

summary(dmel.res)
summary(dmel.res.lfc.ashr)
# result:
# out of 16294 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 7482, 46%
# LFC < 0 (down)     : 1725, 11%
# outliers [1]       : 0, 0%
# low counts [2]     : 632, 3.9%
# these are transformed values?

# summary(dmel.res)
# result:
# out of 16294 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 5616, 34%
# LFC < 0 (down)     : 1133, 7%
# outliers [1]       : 0, 0%
# low counts [2]     : 632, 3.9%
# these are untransformed values


#########################
##### important 
dmel.res.tb <- as.data.frame(dmel.res) %>% 
  tibble::rownames_to_column("gene") %>% 
  arrange(padj)
# making same thing as above
nrow(dmel.res.tb)

dmel.res.sig.tb <- dmel.res.tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
###################################
nrow(dmel.res.sig.tb)

AplGeneData <- plotCounts(dmel.de, gene="Apl", intgroup=c('tissue','sex'), returnData=T,normalized = T)
AplExprPlot <- ggplot(AplGeneData, aes(x = tissue, y = count, color = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  ggtitle("Apl") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
AplExprPlot

ArtsGeneData <- plotCounts(dmel.de, gene="Arts", intgroup=c('tissue','sex'), returnData=T,normalized = T)
ArtsGeneExprPlot <- ggplot(ArtsGeneData, aes(x = tissue, y = count, color = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  ggtitle("Arts") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
ArtsGeneExprPlot

normalized_counts <-  counts(dmel.de, normalized=TRUE) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

dmel.res.sig.tb.names <- dmel.res.sig.tb %>% 
  arrange(padj) %>%	 #Arrange rows by padj values
  pull(gene) 

dmel.res.sig.tb.names

dmel.res.sig.tb.names.norm <- normalized_counts %>%
  filter(gene %in% dmel.res.sig.tb.names)

# Gathering the columns to have normalized counts to a single column
# 
# dmel.res.sig.top50.gathered <- dmel.res.sig.top50.norm %>%
#   gather(colnames(dmel.sig.norm)[2:ncol(dmel.sig.norm)], key = "runs", value = "normalized_counts")
# nrow(dmel.res.sig.top50.gathered)
# # 20 to pgenes * 480 runs = 9600 rows
# # grabbing all 480 runs (sampeles) normalized counts
# # this is why I am splicing 2:ncol()

dmel.res.sig.tb.gathered <- dmel.res.sig.tb.names.norm %>%
  gather(colnames(dmel.res.sig.tb.names.norm)[2:ncol(dmel.res.sig.tb.names.norm)], 
         key = "runs", 
         value = "normalized_counts")

dmel.res.sig.gather.annotated <- inner_join(sampleInfo,dmel.res.sig.tb.gathered)
# 5541 top genes * 416 runs (excluding whole_body per Deanna's advice) = 2305056
# 4396*416 after reloading everything, may have been using old variable for previous number 
nrow(dmel.res.sig.gather.annotated)
# 1828736
# 4396*416 after reloading everything, may have been using old variable for previous number 

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
################################################
################################################
################################################
################################################
######################core for making the data table

dmel.res.sig.tb <- dmel.res.sig.tb %>%
  mutate(threshold = padj < 0.01 & abs(log2FoldChange) >= 0.58)
nrow(dmel.res.sig.tb)
# [1] 4396

dmel.res.sig.gather.annotated <- rename(dmel.res.sig.gather.annotated,'gene' = 'GENEID')
# for joining

#########################################################################################################
################################### Age Information
# initYoung_essential <- read.table('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/proposal/data/essentialityScore/VDRC/newGeneListShengqian/ageData/essentialYoungGeneIDs.txt')
# initAncient_essential <- read.table('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/proposal/data/essentialityScore/VDRC/newGeneListShengqian/ageData/essentialAncientGeneIDs.txt')
initYoung_essential <- read.table('essentialYoungGeneIDs.txt')
initAncient_essential <- read.table('essentialAncientGeneIDs.txt')
names(initYoung_essential)[names(initYoung_essential) == 'V1'] <- 'LOCUSTAG'
names(initAncient_essential)[names(initAncient_essential) == 'V1'] <- 'LOCUSTAG'

initYoung_essential <- initYoung_essential %>% 
  mutate(age = 'young') %>% mutate(essentiality = 'essential') %>% dplyr::select(LOCUSTAG,age,essentiality)

initAncient_essential <- initAncient_essential %>% 
  mutate(age = 'ancient') %>% mutate(essentiality = 'essential') %>% dplyr::select(LOCUSTAG,age,essentiality)
# 
# initYoung_non <- read.table('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/proposal/data/essentialityScore/VDRC/newGeneListShengqian/ageData/nonEssentialYoungGeneIDs.txt')
# initAncient_non <- read.table('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/proposal/data/essentialityScore/VDRC/newGeneListShengqian/ageData/nonEssentialAncientGeneIDs.txt')
initYoung_non <- read.table('nonEssentialYoungGeneIDs.txt')
initAncient_non <- read.table('nonEssentialAncientGeneIDs.txt')
names(initYoung_non)[names(initYoung_non) == 'V1'] <- 'LOCUSTAG'
names(initAncient_non)[names(initAncient_non) == 'V1'] <- 'LOCUSTAG'

initYoung_non <- initYoung_non %>% 
  mutate(age = 'young') %>% mutate(essentiality = 'non-essential') %>% dplyr::select(LOCUSTAG,age,essentiality)

initAncient_non <- initAncient_non %>% 
  mutate(age = 'ancient') %>% mutate(essentiality = 'non-essential') %>% dplyr::select(LOCUSTAG,age,essentiality)

Dmel_IDs <- read.table('GENEID_validated_IDs_final_Dmel.txt')
colnames(Dmel_IDs) <- c('GENEID','FBgn','LOCUSTAG')
# error in naming cols here, i did the wrong order actually in the above files not Dmel_IDs
# now i can join

ancientJoined <- full_join(initAncient_essential,initAncient_non)
youngJoined <- full_join(initYoung_essential, initYoung_non)
# these are total ancient and young
geneAgesAndIDs <- full_join(ancientJoined,youngJoined)

Dmel_geneAgesIDsEssentiality <- left_join(Dmel_IDs,geneAgesAndIDs,'LOCUSTAG')
# 11354 have age and essentially annotation
# the remaining 6487 do not have this information, but do have deseq2 results
# confirmed this is correct now 
# 'correct' given that the essenitality annotation is just based on shengqian table and a crude cut off.
# need to improve this going forward. 
write_tsv(Dmel_geneAgesIDsEssentiality,'Dmel_geneAgesIDsEssentiality.tsv')
################################### Made table of gene IDs and annotations


dmel.res.sig.annotated <- inner_join(dmel.res.sig.gather.annotated,Dmel_geneAgesIDsEssentiality)
#nrow(unique(dmel.res.sig.annotated$GENEID))
# View(dmel.res.sig.annotated)

# write.csv(dmel.res.sig.annotated,
          # file = 'significant_DE_genes')

dmel.res.sig.annotated

sigGenes3 <- ggplot(dmel.res.sig.annotated, aes(x = tissue, y = normalized_counts, color = age, shape = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0),alpha = 0.6) +
  ggtitle("Significantly DE Genes in All Tissues") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
# sigGenes3
ggsave(sigGenes3,filename = 'figures/sigDEgenesAllTissues_MFalpha.png')

dmel.res.sig.annotated.young <- dmel.res.sig.annotated %>% 
  filter(age == 'young') 
# unique(dmel.res.sig.annotated.young$GENEID)
# 334 young genes are in this significant dataset
# nrow(dmel.res.sig.annotated.young)
# [1] 138944
# 138944 / 416 runs = 334 genes
sigYoungGenes <- ggplot(dmel.res.sig.annotated.young, aes(x = tissue, y = normalized_counts, color = essentiality, shape = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0),alpha = 0.6) +
  ggtitle("Significantly DE Young Genes in All Tissues") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
sigYoungGenes
ggsave(sigYoungGenes,filename = 'figures/sigYoungGenes.png')


dmel.res.sig.annotated.ancient <- dmel.res.sig.annotated %>% 
  filter(age == 'ancient') 
sigAncientGenes <- ggplot(dmel.res.sig.annotated.ancient, aes(x = tissue, y = normalized_counts, color = essentiality, shape = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0),alpha = 0.6) +
  ggtitle("Significantly DE Ancient Genes in All Tissues") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
sigAncientGenes
ggsave(sigAncientGenes,filename = 'figures/sigAncientGenes.png')


dmel.res.sig.annotated.essential <- dmel.res.sig.annotated %>% 
  filter(essentiality == 'essential') 
sigEssentialGenes <- ggplot(dmel.res.sig.annotated.essential, aes(x = tissue, y = normalized_counts, color = age, shape = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0),alpha = 0.6) +
  ggtitle("Significantly DE Essential Genes in All Tissues") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
sigEssentialGenes
ggsave(sigEssentialGenes,filename = 'figures/sigEssentialGenes.png')

################################################
######################core for making the data table for volcano
################################################################################################

dmel.res.tb <- dmel.res.tb %>%
  mutate(threshold = padj < 0.01 & abs(log2FoldChange) >= 0.58)


## Create a column to indicate which genes to label
dmel.res.tb <- dmel.res.tb %>% arrange(padj) %>% mutate(genelabels = "")
dmel.res.tb$genelabels <- dmel.res.tb$gene
dmel.res.tb.truncated <- dmel.res.tb %>% filter(padj < 0.1)
# just to plot fewer ?
nrow(dmel.res.tb.truncated)
dmel.res.tb.truncated <- rename(dmel.res.tb.truncated,'gene' = 'GENEID')
dmel.res.tb.truncated.annotated <- inner_join(dmel.res.tb.truncated,Dmel_geneAgesIDsEssentiality)
nrow(dmel.res.tb.truncated.annotated)

volcanoLog10 <- ggplot(dmel.res.tb.truncated.annotated, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = age, shape = essentiality)) +
  geom_point(data = dmel.res.tb.truncated.annotated[dmel.res.tb.truncated.annotated$threshold == "FALSE",], color = "orange", alpha=0.6) +
  ggtitle("Dmel expression (orange is non-significant Padj)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
volcanoLog10
ggsave(volcanoLog10,filename = 'figures/volcanoLog10.png')
# all genes that are significant
################################################################################################

# save.image(f = 'dmelRNAseq.Rdata')
# load('dmelRNAseq.Rdata')
################################################
# # Density plotting 

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)  

# ggplot(dmel.res.sig.annotated %>% filter(tissue =='gonad')) +
dmel.res.sig.annotated.group <- na.omit(dmel.res.sig.annotated) %>% mutate(group = paste0(age,'_',essentiality))
# nrow(na.omit(dmel.res.sig.annotated))
nrow(dmel.res.sig.annotated)
##### 1828736/416 = 4396


# dmel.res.sig.tb 
dmel.res.sig.annotated.notgathered <- inner_join(rename(dmel.res.sig.tb,gene='GENEID'),Dmel_geneAgesIDsEssentiality)
nrow(na.omit(dmel.res.sig.annotated.notgathered))
# [1] 4396
# 3333 when removing NA. This means those genes don't have age or essentiality annotation I believe
# this tells me the number of genes I have that are signficiantly express 
# padj < 0.01 & abs(log2FoldChange) >= 0.58
# ^ this number is based on the above line criteria

dmel.res.sig.annotated.notgathered %>% filter(age == 'young') %>%  nrow()
# [1] 334

dmel.res.sig.annotated.notgathered %>% filter(age == 'ancient') %>% nrow()
# [1] 2999

dmelGenesCountDensity <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                  mutate(group = paste0(age,'_',essentiality))) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel expression in various tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
  ylim(0,0.4) +
  facet_grid(. ~ tissue)
  # # this one 
dmelGenesCountDensity

youngGenesCountDensity <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                   mutate(group = paste0(age,'_',essentiality)) %>% 
                                   filter(age=='young')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel young gene (n=334) expression in various tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ tissue)
youngGenesCountDensity

ancientGenesCountDensity <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                     mutate(group = paste0(age,'_',essentiality)) %>% 
                                     filter(age=='ancient')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel ancient gene (n=2999) expression in various tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ tissue)
ancientGenesCountDensity

grid.arrange(dmelGenesCountDensity,youngGenesCountDensity,ancientGenesCountDensity,ncol=1)
dmelDensityGrid <- grid.arrange(dmelGenesCountDensity,youngGenesCountDensity,ancientGenesCountDensity,ncol=1)
ggsave(dmelDensityGrid, filename = 'figures/dmelDensityGrid.png')
#####################################################################################################################################
####### MAKING SEX SPECIFIC DENSITY PLOTS
youngGenesCountDensityMale <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                       mutate(group = paste0(age,'_',essentiality)) %>% 
                                       filter(age=='young') %>% 
                                       filter(sex=='male')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel expression in male sex tissues') +
  ylim(0,0.4) +
  facet_grid(.  ~ c('gonad','genitalia'))
youngGenesCountDensityFemale <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                         mutate(group = paste0(age,'_',essentiality)) %>% 
                                         filter(age=='young')%>% 
                                         filter(sex=='female')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel expression in female sex tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ c('gonad','genitalia'))

ancientGenesCountDensityMale <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                       mutate(group = paste0(age,'_',essentiality)) %>% 
                                       filter(age=='ancient') %>% 
                                       filter(sex=='male')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel expression in male sex tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ c('gonad','genitalia'))
ancientGenesCountDensityFemale <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                         mutate(group = paste0(age,'_',essentiality)) %>% 
                                         filter(age=='ancient')%>% 
                                         filter(sex=='female')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel expression in female sex tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ c('gonad','genitalia'))

dmelGenesCountDensitySexTissue <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                  mutate(group = paste0(age,'_',essentiality))) +
  aes(x =log2(normalized_counts+0.05), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel expression in sex tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
  ylim(0,0.4) +
  facet_grid(. ~ c('gonad','genitalia'))

youngGenesCountDensitySexTissue <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                   mutate(group = paste0(age,'_',essentiality)) %>% 
                                   filter(age=='young')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel young gene (n=334) expression in sex tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ c('gonad','genitalia'))
# youngGenesCountDensitySexTissue

ancientGenesCountDensitySexTissue <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
                                     mutate(group = paste0(age,'_',essentiality)) %>% 
                                     filter(age=='ancient')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel ancient gene (n=2999) expression in sex tissues') +
  ylim(0,0.4) +
  facet_grid(. ~ c('gonad','genitalia'))
# ancientGenesCountDensitySexTissue

# resultsNames(dmel.de)
             
dmelYoungGeneSexTissueGridMF <- grid.arrange(youngGenesCountDensitySexTissue,youngGenesCountDensityMale,youngGenesCountDensityFemale)
dmelAncientGeneSexTissueGridMF <- grid.arrange(ancientGenesCountDensitySexTissue,ancientGenesCountDensityMale,ancientGenesCountDensityFemale)

ggsave(dmelYoungGeneSexTissueGridMF, filename = 'figures/dmelYoungGeneSexTissueGridMF.png')
ggsave(dmelAncientGeneSexTissueGridMF, filename = 'figures/dmelAncientGeneSexTissueGridMF.png')


dmelGenesCountDensityBySex <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
         mutate(group = paste0(age,'_',essentiality))) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel sex-based expression various tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
  ylim(0,0.4) +
  facet_grid(sex ~ tissue)

dmelCountDensityAllAndBySex <- grid.arrange(dmelGenesCountDensity, dmelGenesCountDensityBySex)
ggsave(dmelCountDensityAllAndBySex, filename = 'figures/dmelCountDensityAllAndBySex.png')

dmelGenesCountDensityBySex
ggsave(dmelGenesCountDensityBySex,f='figures/dmelGenesCountDensityBySex.png')




############################
# for lab meeting 2/4
atlasGeneData <- plotCounts(dmel.de, gene="CG13541", intgroup=c('tissue','sex'), returnData=T,normalized = T)
atlasExprPlot <- ggplot(atlasGeneData, aes(x = tissue, y = count, color = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  ggtitle("atlas gene (CG13541) expression") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
atlasExprPlot
ggsave(atlasExprPlot,f='figures/atlasExprPlot.png')

atlas.volcanoLog10 <- ggplot(dmel.res.tb.truncated.annotated, aes(x = log2FoldChange, y = -log10(padj),label=LOCUSTAG)) +
  geom_point(aes(colour = age, shape = essentiality)) +
  geom_point(data = dmel.res.tb.truncated.annotated[dmel.res.tb.truncated.annotated$threshold == "FALSE",], color = "orange", alpha=0.6) +
  ggtitle("atlas gene (CG13541) expression, Dmel expression (orange is non-significant Padj)") +
  geom_label_repel(data = subset(dmel.res.tb.truncated.annotated, LOCUSTAG == "CG13541"), 
                  aes(label = LOCUSTAG), color = "black",
                  box.padding = 0.5, 
                  max.overlaps = Inf,
                  point.padding = 0.2,
                  nudge_x = .15,
                  nudge_y = .5,
                  segment.angle = 20) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
atlas.volcanoLog10
ggsave(atlas.volcanoLog10,f='figures/atlas.volcanoLog10.png')
#############################
#### for manyuan, individual plots showing sex expression for each of the four gene categories 
# geneGroups <- na.omit(dmel.res.sig.annotated) %>%
#   mutate(group = paste0(age,'_',essentiality)) %>%
#   distinct(.,group)
# geneGroups

# ggplot(na.omit(dmel.res.sig.annotated) %>%
#            mutate(group = paste0(age,'_',essentiality))) +
#     aes(x =log2(normalized_counts+1), fill=group=='young_essential') +
#     geom_density(alpha=0.35) +
#     scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
#     labs(title = 'Dmel sex-based expression various tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
#     ylim(0,0.4) +
#     facet_grid(sex ~ tissue)
# 
# ggplot(na.omit(dmel.res.sig.annotated) %>%
#          mutate(group = paste0(age,'_',essentiality))) +
#   aes(x =log2(normalized_counts+1), fill='ancient_essential') +
#   geom_density(alpha=0.35) +
#   scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
#   labs(title = 'Dmel sex-based expression invarious tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
#   ylim(0,0.4) +
#   facet_grid(sex ~ tissue)
# 
# ggplot(na.omit(dmel.res.sig.annotated) %>%
#          mutate(group = paste0(age,'_',essentiality))) +
#   aes(x =log2(normalized_counts+1), fill='young_non-essential') +
#   geom_density(alpha=0.35) +
#   scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
#   labs(title = 'Dmel sex-based expression various tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
#   ylim(0,0.4) +
#   facet_grid(sex ~ tissue)
# 
# ggplot(na.omit(dmel.res.sig.annotated) %>%
#          mutate(group = paste0(age,'_',essentiality))) +
#   aes(x =log2(normalized_counts+1), fill='young_essential') +
#   geom_density(alpha=0.35) +
#   scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
#   labs(title = 'Dmel sex-based expression various tissues of genes (n=3333) with padj < 0.01 & abs(log2FoldChange) >= 0.58') +
#   ylim(0,0.4) +
#   facet_grid(sex ~ tissue)


#########
# the two plots I sent to Manyuan:
youngGeneTissueExpressionMaleFemale <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
         mutate(group = paste0(age,'_',essentiality)) %>% 
         filter(age=='young')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel young gene (n=334) expression in 8 tissues') +
  ylim(0,0.4) +
  facet_grid(sex ~ tissue)
ggsave(youngGeneTissueExpressionMaleFemale,f='figures/forManyuan/youngGeneTissueExpressionMaleFemale.png')

ancientGeneTissueExpressionMaleFemale <- ggplot(na.omit(dmel.res.sig.annotated) %>% 
         mutate(group = paste0(age,'_',essentiality)) %>% 
         filter(age=='ancient')) +
  aes(x =log2(normalized_counts+1), fill=group) +
  geom_density(alpha=0.35) +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  labs(title = 'Dmel ancient gene (n=2999) expression in 8 tissues') +
  ylim(0,0.4) +
  facet_grid(sex ~ tissue)
ggsave(ancientGeneTissueExpressionMaleFemale,f='figures/forManyuan/ancientGeneTissueExpressionMaleFemale.png')


save.image(f = 'dmelRNAseq_relaxed.Rdata')
load('dmelRNAseq_relaxed.Rdata')


dmel.results <- dmel.res.sig.annotated.notgathered
write_tsv(dmel.results,'dmel_significant_results_4396.tsv')
