## Here is the R code for processing the temporal transcriptomics of Dehnel's phenomenon cortex
#Reads were preprocessed on Noctillio server, where they were trimmed using fastp and aligned/quantified using Kalisto
#First I want to get the quantification files onto R and convert transcript abundance to gene abundance

#First set up all your libraries
library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(readr)
library(tibble)
library( "genefilter" )
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(TCseq)
library(cluster)
library(EnhancedVolcano)

## STEP1: Import quant files
#Set up your working directories
indir <- "./script"
outputdir <- "../data/Liver"

#Create a transcript to gene key using the gff file for the shrew
sortxdb <- makeTxDbFromGFF("..data//refs/GCF_000181275.1_SorAra2.0_genomic.gtf.gz")
keytypes(sortxdb)k <- keys(sortxdb, keytype="TXNAME")       
sortx2gene <- select(sortxdb, k, "GENEID", "TXNAME")
sortx2gene <- sortx2gene[!duplicated(sortx2gene[,1]),]
sortx2gene <- na.omit(sortx2gene)
#gets rid of the XRs, which are some misc_rnas, and do not apppear to be associated with genes

###READ IN ABUNDANCE FILES AND CHANGE THEM TO GENE ABUNDANCES
cs_cor_samples <- read.table("../data/Cortex/SampleList", header = T)
cs_cor_files <-file.path("../data/Cortex/TranscriptAbundances/", cs_cor_samples$Sample_name, "abundance.tsv")
names(cs_cor_files) <- paste0("sample_", cs_cor_samples$Sample_name)
all(file.exists(cs_cor_files))
#write to tables and can save with hashed out code below, but note will overwrite
cs_cor.count.tsv <- tximport(cs_cor_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_cor.tpm.tsv <- tximport(cs_cor_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
#write.table(cs_cor.tpm.tsv$abundance, "../data/Cortex/GeneCounts/GeneCounts.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

## STEP2: Graph brain mass for panel A
#Individual region analysis of volume is being done with MRI for sepaate analysis, brain resolution here
SampleData_Rreadable <- read_csv("~/Dehnels_Seasonal_RNAseq/data/Liver/SampleData_Rreadable.csv")
brainmass <- SampleData_Rreadable$Brain_Mass
stagemass2 <- c(rep(0,5),rep(5,4),rep(8,5),rep(10,5),rep(12,5))
stage4mass <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
bmass_frame <- as.data.frame(cbind(brainmass,stagemass2,stage4mass))
bmass_frame$brainmass <-  as.numeric(bmass_frame$brainmass)
gg4<-ggplot(bmass_frame, aes(x = stage4mass, y = brainmass)) +
  geom_boxplot(size=.75)+
  xlab("Month")+
  ylab("BrainMass")+
  theme_bw()

## STEP3: Normalization and quality control
colnames(cs_cor.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
cor_stages <- factor(c(cs_cor_samples$Run))
cor_organs <- factor(c(cs_cor_samples$Condition))
cor_full1 <- factor(c(cs_cor_samples$Run))
cor_stages_organ_frame <-cbind(as.data.frame(cor_stages),as.data.frame(cor_organs),as.data.frame(cor_full1))
vis_cor_samples <- SampleData_Rreadable
ggplot(vis_cor_samples,aes(x=Cortex_RIN,y=Cortex_Reads_prefilt, color= Stage))+
  geom_point()+
  theme_bw()
ggplot(vis_cor_samples,aes(x=Cortex_RIN,y=Cortex_Reads_postfilt, color= Stage))+
  geom_point()+
  theme_bw()
ggplot(vis_cor_samples,aes(x=Cortex_RIN,y=Cortex_Pseudo_aligned_Percent, color= Stage))+
  geom_point()+
  theme_bw()
#everything looks pretty good
dds_cor_all <- DESeqDataSetFromMatrix(round(cs_cor.count.tsv$counts), DataFrame(cor_stages_organ_frame), ~ cor_full1)
mcols(dds_cor_all) <- cbind(mcols(dds_cor_all), row.names(cs_cor.count.tsv$counts))
rownames(dds_cor_all) <- row.names(cs_cor.count.tsv$counts)
dds_cor_all <- DESeq(dds_cor_all)
#pca and heatmap cortex
vst_dds_cor_all <- vst(dds_cor_all)
pcaData_cor_all<- plotPCA(vst_dds_cor_all,intgroup=c("cor_stages","cor_organs"), ntop=100, returnData=TRUE)
ggplot(pcaData_cor_all, aes(x = PC1, y = PC2, color = factor(cor_stages))) +
  geom_point(size=2)+
  theme_bw()
corsampleDists <- dist(t(assay(vst_dds_cor_all)))
corsampleDistMatrix <- as.matrix(corsampleDists)
colnames(corsampleDistMatrix) <- NULL
##make the heatmap cortex
pheatmap(corsampleDistMatrix, clustering_distance_rows=corsampleDists,
         clustering_distance_cols = corsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
cortopVarGenes <- head( order( rowVars( assay(vst_dds_cor_all) ), decreasing=TRUE ), 50 )
vst_dds_cor_all[ cortopVarGenes, ]
heatmap.2( assay(vst_dds_cor_all)[ cortopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

#STEP 4 : Run differential expression on stage 1 vs 3
lfc<-1.58
cor13res <- results(dds_cor_all, contrast = c("cor_full1","Stage3","Stage1"))
cor13resSig <- subset(cor13res,cor13res$padj<.05)
cor13resSigLog <- subset(cor13resSig,abs(cor13resSig$log2FoldChange)>=lfc)
cor13up <- subset(cor13resSigLog,(cor13resSigLog$log2FoldChange)>=0)
cor13down <- subset(cor13resSigLog,(cor13resSigLog$log2FoldChange)<=0)
####CREATE VOLCANO PLOT FOR FIGURE
cor13resX <-  cor13res
for (i in 1:length(cor13resX$padj)) {
  if  (cor13resX$padj[i]<1e-5 & !is.na (cor13resX$padj[i])) {
    cor13resX$padj[i] <- 1e-5
  }
  if (cor13resX$log2FoldChange[i]>4 & !is.na (cor13resX$log2FoldChange[i])) {
    cor13resX$log2FoldChange[i] <- 4
  }
  if (cor13resX$log2FoldChange[i]< -4 & !is.na (cor13resX$log2FoldChange[i])) {
    cor13resX$log2FoldChange[i] <- -4
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(cor13resX$log2FoldChange) == 4, 17,
  ifelse(cor13resX$padj==1e-5, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  cor13resX$padj > 0.05, 'grey',
  ifelse(cor13resX$log2FoldChange <= -1.58, 'red',
         ifelse(cor13resX$log2FoldChange >= 1.58, 'blue',
                ifelse(cor13resX$log2FoldChange >= 0, 'lightblue',
                       'pink'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'pink'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'lightblue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
####
cor13resX[rownames(hip34resX) %in% "APOA1",]
EnhancedVolcano(cor13resX,
                lab = rownames(cor13resX),
                xlim=c(-4 ,4),
                ylim=c(0,5),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = c('APOH','APOA1',"FAM163B","HPD","HSD11B1","LOC101542466","LOC101545806", "LOC101553347","MADCAM1","PCK1","PIGR"),
                #selectLab = c("SIRT1",'CREBBP','FOXO','G6PC','PPARA','RXRG','FOXO1','RXRA','PER1','PPARD','MBP','PLP','FABP5','SLC27A2','FGF21'),
                #selectLab = c('SIRT1','CREBBP','FOXO3','G6PC','PPARA','RXRG','FOXO1','RXRA','PER1','PPARD','MBP','PLP','FABP5','SLC27A2','FGF21','SIRT1'),
                pCutoff = .05,
                FCcutoff = 1.58,
                colCustom = keyvals,
                legendPosition = 'right',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.75)

#STEP5: Run TCSeq
#Timeseq cortex with LFC .5849 (*1.5 change)
con_cor3 <- c(rep("Stg1",5),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
cor3timeseq <- data.frame(sampleID = 1:24, group = c(1, 1, 1,1,1,2, 2, 2,2, 3, 3, 3,3,3,4, 4, 4, 4,4,5, 5, 5, 5,5),
                          timepoint = con_cor3)
gf <- data.frame(chr = c(rep('chr1', 15296), rep('chr2', 2000), rep('chr4', 2000)),
                 start = rep(100, 19296),
                 end = rep(500, 19296),
                 id = row.names(cs_cor.tpm.tsv$abundance))
cor3tca <- TCA(design = cor3timeseq, counts = round(DESeq2::counts(dds_cor_all, normalized=TRUE)), gf)
cor3tca <- DBanalysis(cor3tca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)
cor3tca <- timecourseTable(cor3tca, value = "expression",  lib.norm = FALSE, filter = TRUE,abs.fold = 0.5)
cor3_t <- tcTable(cor3tca)
cor3tca
#19k to 14195 to 940
#figure out how many clusters there likely are
# using timeclust2 with clusGap
clusGap(cor3_t,
        FUNcluster = timeclust2,
        algo = 'cm',
        K.max = 20,
        B = 20)
##12 for cortex .5LFC
ccor3tca <- timeclust(cor3tca, algo = "cm", k = 12, standardize = TRUE)
XXXccor3_px <-timeclustplot(ccor3tca, value = "z-score(PRKM)", cols = 2)
cor3cxxx<-clustResults(ccor3tca)

#### where do the significant genes reside
memb_matrix <- cor3cxxx@membership
mem_col <- as.data.frame(colnames(memb_matrix)[max.col(memb_matrix,ties.method="first")])
rownames(mem_col) <- rownames(memb_matrix)
mem_col_subset <-  mem_col[rownames(cor13resSigLog) %in% rownames(mem_col), ]
mem_col[rownames(mem_col) %in% rownames(cor13resSigLog), ]
cluster_cor13resSiglog <- subset(mem_col, rownames(mem_col) %in% rownames(cor13resSigLog))
colnames(cluster_cor13resSiglog) <- "membership"
cluster_cor13resSiglog$membership <- as.factor(cluster_cor13resSiglog$membership)
cluster_cor13resSiglog$membership <- factor(cluster_cor13resSiglog$membership, levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))
ggplot(cluster_cor13resSiglog, aes(x=membership))+
  geom_bar()+
  theme_bw()
#trying slightly different way to plot better
dribble <- cluster_cor13resSiglog %>%
  add_column(UpDown = NA)
dribble <- as.data.frame(dribble)
dribble
for(i in 1:length(cluster_cor13resSiglog$membership)) {
  xcc <- rownames(dribble)
  xcc <- xcc[i]
  dcc5 <-cor13resSig[rownames(cor13resSig) %in% xcc, ]
  LFCC <- dcc5$log2FoldChange
  print(LFCC)
  if (LFCC > 0) {
    dribble$UpDown[i] <- "Up"
  } else {
    dribble$UpDown[i] <- "Down"
  }
}
dribble
#
dribble$membership <- factor(dribble$membership, levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))
ggplot(dribble, aes(x=UpDown))+
  geom_bar(aes(fill = membership), width = 0.5)+ 
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values =c("#23000000", "#23E69F00", "#2356B4E9", "#23009E73", "#23F0E442", "#230072B2", "#23D55E00", "#23CC79A7", "#2331cefc", "#237ed239", "#23b46357", "#235c5f18", "#23c56b49"))+
  theme_bw()

#STEP 7: Figure
###READ IN ABUNDANCE FILES AND CHANGE THEM TO GENE ABUNDANCES for hippocampus
cs_hip_samples <- read.table("../data/Hippocampus/SampleList", header = T)
cs_hip_files <-file.path("../data/Hippocampus/TranscriptAbundances/", cs_hip_samples$Sample_name, "abundance.tsv")
names(cs_hip_files) <- paste0("sample_", cs_hip_samples$Sample_name)
all(file.exists(cs_hip_files))
#write to tables and can save with hashed out code below, but note will overwrite
cs_hip.count.tsv <- tximport(cs_hip_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_hip.tpm.tsv <- tximport(cs_hip_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
colnames(cs_hip.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
hip_stages <- factor(c(cs_hip_samples$Run))
hip_organs <- factor(c(cs_hip_samples$Condition))
hip_full1 <- factor(c(cs_hip_samples$Run))
hip_stages_organ_frame <-cbind(as.data.frame(hip_stages),as.data.frame(hip_organs),as.data.frame(hip_full1))
dds_hip <- DESeqDataSetFromMatrix(round(cs_hip.count.tsv$counts), DataFrame(hip_stages_organ_frame), ~ hip_full1_no1)
mcols(dds_hip) <- cbind(mcols(dds_hip), row.names(cs_hip.count.tsv$counts))
rownames(dds_hip) <- row.names(cs_hip.count.tsv$counts)
dds_hip <- DESeq(dds_hip)
###APOH###
a_cor <- DESeq2::counts(dds_cor_all, normalized=TRUE)
a_hip <- DESeq2::counts(dds_hip_all, normalized=TRUE)
apoh_cortex <-a_cor[rownames(a_cor) %in% "APOH", ]
apoh_hippoc <-a_hip[rownames(a_hip) %in% "APOH", ]
apoh <- c(apoh_cortex,apoh_hippoc)
corhip <- c(rep("Cortex",24),rep("Hippoc",23))
corhip_stage <- c(rep("Stg1",5),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5),rep("Stg1",4),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
corhip_frame_apoh <-cbind(as.data.frame(apoh),as.data.frame(corhip),as.data.frame(corhip_stage))
ggplot(corhip_frame_apoh, aes(x = corhip_stage, y = apoh, color = corhip)) +
  geom_boxplot(size=1)+
  scale_color_manual(values=c("#0b5da2ff","#E69F00"))+
  theme_bw()
###APOA1
apoa1_cortex <-a_cor[rownames(a_cor) %in% "APOA1", ]
apoa1_hippoc <-a_hip[rownames(a_hip) %in% "APOA1", ]
apoa1 <- c(apoa1_cortex,apoa1_hippoc)
corhip <- c(rep("Cortex",24),rep("Hippoc",23))
corhip_stage <- c(rep("Stg1",5),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5),rep("Stg1",4),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
corhip_frame_apoh <-cbind(as.data.frame(apoa1),as.data.frame(corhip),as.data.frame(corhip_stage))
ggplot(corhip_frame_apoa1, aes(x = corhip_stage, y = apoa1, color = corhip)) +
  geom_boxplot(size=1)+
  scale_color_manual(values=c("#0b5da2ff","#E69F00"))+
  theme_bw()