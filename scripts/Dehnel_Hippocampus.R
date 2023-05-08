## Here is the R code for processing the temporal transcriptomics of Dehnel's phenomenon hippocampus
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
library(WGCNA)
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
cs_hip_samples <- read.table("../data/Hippocampus/SampleList", header = T)
cs_hip_files <-file.path("../data/Hippocampus/TranscriptAbundances/", cs_hip_samples$Sample_name, "abundance.tsv")
names(cs_hip_files) <- paste0("sample_", cs_hip_samples$Sample_name)
all(file.exists(cs_hip_files))
#write to tables and can save with hashed out code below, but note will overwrite
cs_hip.count.tsv <- tximport(cs_hip_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_hip.tpm.tsv <- tximport(cs_hip_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
#write.table(cs_hip.tpm.tsv$abundance, "../data/Hippocampus/GeneCounts/GeneCounts.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

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

## STEP3: Normalization and QC checks
#Begin by setting up design matrix
colnames(cs_hip.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
hip_stages <- factor(c(cs_hip_samples$Run))
hip_organs <- factor(c(cs_hip_samples$Condition))
hip_full1 <- factor(c(cs_hip_samples$Run))
hip_stages_organ_frame <-cbind(as.data.frame(hip_stages),as.data.frame(hip_organs),as.data.frame(hip_full1))
ggplot(hip_stages_organ_frame, aes(x=hip_full1))+
  geom_bar(stat = 'count')+
  theme_bw()
#Hippocampus RNA Integrity
vis_hip_samples <-SampleData_Rreadable
ggplot(vis_hip_samples,aes(x=Hippocampus_RIN,y=Hippocampus_Reads_prefilt, color= Stage))+
  geom_point()+
  theme_bw()
#Remove outlier reads
cs_hip.count_NoOUTLIER_STG1.tsv <- cs_hip.count.tsv
cs_hip.count_NoOUTLIER_STG1.tsv
cs_hip.count_NoOUTLIER_STG1.tsv$counts <- subset(cs_hip.count_NoOUTLIER_STG1.tsv$counts, select=-sample_HFS156)
colnames(cs_hip.count_NoOUTLIER_STG1.tsv$counts) <- c("Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
hip_stages_noout <- hip_stages[-1]
hip_organs_noout <- hip_organs[-1]
hip_full1_noout <- hip_full1[-1]
#Revisualize RINs without sample 156
vis_hip_samples <- vis_hip_samples[-1,]
ggplot(vis_hip_samples,aes(x=Hippocampus_RIN,y=Hippocampus_Reads_prefilt, color= Stage))+
  geom_point()+
  theme_bw()
ggplot(vis_hip_samples,aes(x=Hippocampus_RIN,y=Hippocampus_Reads_postfilt, color= Stage))+
  geom_point()+
  theme_bw()
ggplot(vis_hip_samples,aes(x=Hippocampus_RIN,y=Hippocampus_Pseudo_aligned_Percent, color= Stage))+
  geom_point()+
  theme_bw()

###APPARENT MAPPING BIAS IN STAGE1 OR AGGRESSIVE STAGE 1 EFFECT, LETS LOOK AT A PCA 
hip_stages_organ_frame_noout <-cbind(as.data.frame(hip_stages_noout),as.data.frame(hip_organs_noout),as.data.frame(hip_full1_noout))
dds_hip_all_noout <- DESeqDataSetFromMatrix(round(cs_hip.count_NoOUTLIER_STG1.tsv$counts), DataFrame(hip_stages_organ_frame_noout), ~ hip_full1_noout)
mcols(dds_hip_all_noout) <- cbind(mcols(dds_hip_all_noout), row.names(cs_hip.count_NoOUTLIER_STG1.tsv$counts))
rownames(dds_hip_all_noout) <- row.names(cs_hip.count_NoOUTLIER_STG1.tsv$counts)
dds_hip_all_noout <- DESeq(dds_hip_all_noout)
vst_dds_hip_all <- DESeq2::vst(dds_hip_all_noout)
pcaData_hip_all_noout<- plotPCA(vst_dds_hip_all,intgroup=c("hip_stages_noout","hip_organs_noout"), ntop=300, returnData=TRUE)
ggplot(pcaData_hip_all_noout, aes(x = PC1, y = PC2, color = factor(hip_stages_noout))) +
  geom_point(size=3)+
  theme_bw()
###CONCLUSION: EITHER MAPPING EFFECT OR DEVELOPMENTAL EFFECT AS STAGE 1 VERY DISTANT TIME WISE FROM REST
###ACTIONABLE; REMOVE FROM DESEQ AS WILL INTRODUCE TOO MUCH VARIATION FROM MAPPING

#REMOVE STAGE 1
cs_hip.count_NoOUTLIER_NoSTG1.tsv <- cs_hip.count.tsv
colnames(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts
cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts <- subset(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts, select=-Stg1_1)
cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts <- subset(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts, select=-Stg1_2)
cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts <- subset(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts, select=-Stg1_3)
cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts <- subset(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts, select=-Stg1_4)
cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts <- subset(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts, select=-Stg1_5)
hip_stages_no1 <- factor(c(cs_hip_samples$Run))
#Gross code below, update at some point but need to remove all stage 1s
hip_stages_no1 <- hip_stages_no1[-1]
hip_organs_no1 <- factor(c(cs_hip_samples$Condition))
hip_stages_no1 <- hip_stages_no1[-1]
hip_organs_no1 <- factor(c(cs_hip_samples$Condition))
hip_stages_no1 <- hip_stages_no1[-1]
hip_organs_no1 <- factor(c(cs_hip_samples$Condition))
hip_stages_no1 <- hip_stages_no1[-1]
hip_organs_no1 <- factor(c(cs_hip_samples$Condition))
hip_stages_no1 <- hip_stages_no1[-1]
hip_organs_no1 <- factor(c(cs_hip_samples$Condition))
#Do x 5
hip_organs_no1 <- hip_organs_no1[-1]
hip_full1_no1 <- factor(c(cs_hip_samples$Run))
hip_organs_no1 <- hip_organs_no1[-1]
hip_full1_no1 <- factor(c(cs_hip_samples$Run))
hip_organs_no1 <- hip_organs_no1[-1]
hip_full1_no1 <- factor(c(cs_hip_samples$Run))
hip_organs_no1 <- hip_organs_no1[-1]
hip_full1_no1 <- factor(c(cs_hip_samples$Run))
hip_organs_no1 <- hip_organs_no1[-1]
hip_full1_no1 <- factor(c(cs_hip_samples$Run))
#Do x 5
hip_full1_no1 <- hip_full1_no1[-1]
hip_stages_organ_frame_no1 <-cbind(as.data.frame(hip_stages_no1),as.data.frame(hip_organs_no1),as.data.frame(hip_full1_no1))
hip_full1_no1 <- hip_full1_no1[-1]
hip_stages_organ_frame_no1 <-cbind(as.data.frame(hip_stages_no1),as.data.frame(hip_organs_no1),as.data.frame(hip_full1_no1))
hip_full1_no1 <- hip_full1_no1[-1]
hip_stages_organ_frame_no1 <-cbind(as.data.frame(hip_stages_no1),as.data.frame(hip_organs_no1),as.data.frame(hip_full1_no1))
hip_full1_no1 <- hip_full1_no1[-1]
hip_stages_organ_frame_no1 <-cbind(as.data.frame(hip_stages_no1),as.data.frame(hip_organs_no1),as.data.frame(hip_full1_no1))
hip_full1_no1 <- hip_full1_no1[-1]
hip_stages_organ_frame_no1 <-cbind(as.data.frame(hip_stages_no1),as.data.frame(hip_organs_no1),as.data.frame(hip_full1_no1))
#everything looks pretty good
dds_hip_all_no1 <- DESeqDataSetFromMatrix(round(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts), DataFrame(hip_stages_organ_frame_no1), ~ hip_full1_no1)
mcols(dds_hip_all_no1) <- cbind(mcols(dds_hip_all_no1), row.names(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts))
rownames(dds_hip_all_no1) <- row.names(cs_hip.count_NoOUTLIER_NoSTG1.tsv$counts)
dds_hip_all_no1 <- DESeq(dds_hip_all_no1)
#hippocampus pca
vst_dds_hip_all <- DESeq2::vst(dds_hip_all_no1)
pcaData_hip_all_no1<- plotPCA(vst_dds_hip_all,intgroup=c("hip_stages_no1","hip_organs_no1"), ntop=300, returnData=TRUE)
ggplot(pcaData_hip_all_no1, aes(x = PC1, y = PC2, color = factor(hip_stages_no1))) +
  geom_point(size=3)+
  theme_bw()
###HippocampusHeatmap###
hipsampleDists <- dist(t(assay(vst_dds_hip_all)))
hipsampleDistMatrix <- as.matrix(hipsampleDists)
colnames(hipsampleDistMatrix) <- NULL
##make the heatmap
pheatmap(hipsampleDistMatrix, clustering_distance_rows=hipsampleDists,
         clustering_distance_cols = hipsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
##trying to make a heatmap of genes
hiptopVarGenes <- head( order( rowVars( assay(vst_dds_hip_all) ), decreasing=TRUE ), 300 )
vst_dds_hip_all[ hiptopVarGenes, ]
heatmap.2( assay(vst_dds_hip_all)[ hiptopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

###STEP 4: RUN DESEQ
#no lfc
lfc <- 0
#Hippocampus 3-4
hip34res <- results(dds_hip_all_no1, contrast = c("hip_full1_no1","Stage4","Stage3"))
hip34resSig <- subset(hip34res,hip34res$padj<.05)
hip34resSigLog <- subset(hip34resSig,abs(hip34resSig$log2FoldChange)>=lfc)
hip34up <- subset(hip34resSigLog,(hip34resSigLog$log2FoldChange)>=0)
hip34down <- subset(hip34resSigLog,(hip34resSigLog$log2FoldChange)<=0)
#plot with volcano plot
hip34resX <-  hip34res
library(EnhancedVolcano)
for (i in 1:length(hip34resX$padj)) {
  if  (hip34resX$padj[i]<1e-5 & !is.na (hip34resX$padj[i])) {
    hip34resX$padj[i] <- 1e-5
  }
  if (hip34resX$log2FoldChange[i]>2 & !is.na (hip34resX$log2FoldChange[i])) {
    hip34resX$log2FoldChange[i] <- 2
  }
  if (hip34resX$log2FoldChange[i]< -2 & !is.na (hip34resX$log2FoldChange[i])) {
    hip34resX$log2FoldChange[i] <- -2
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(hip34resX$log2FoldChange) == 2, 17,
  ifelse(hip34resX$padj==1e-5, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  hip34resX$padj > 0.05, 'grey',
  ifelse(hip34resX$log2FoldChange <= -1.58, 'red',
         ifelse(hip34resX$log2FoldChange >= 1.58, 'blue',
                ifelse(hip34resX$log2FoldChange >= 0, 'lightblue',
                       'pink'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'pink'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'lightblue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
####
hip34resX[rownames(hip34resX) %in% "APOA1",]
EnhancedVolcano(hip34resX,
                lab = rownames(hip34resX),
                xlim=c(-2 ,2),
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

###STEP 5: RUN TimeSeq on Data
###Timeseq Hippocampus with LFC .5849 (abs 1.5 change) and no Stage 1
con_hip3 <- c(rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
hip3timeseq <- data.frame(sampleID = 1:19, group = c(2, 2, 2, 2,3, 3, 3,3,3,4, 4, 4, 4,4,5, 5, 5, 5,5),
                          timepoint = con_hip3)
gf <- data.frame(chr = c(rep('chr1', 15296), rep('chr2', 2000), rep('chr4', 2000)),
                 start = rep(100, 19296),
                 end = rep(500, 19296),
                 id = row.names(cs_hip.tpm.tsv$abundance))
hip3tca <- TCA(design = hip3timeseq, counts = round(DESeq2::counts(dds_hip_all_no1, normalized=TRUE)), gf)
hip3tca <- DBanalysis(hip3tca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)
hip3tca <- timecourseTable(hip3tca, value = "expression",  lib.norm = FALSE, filter = TRUE,abs.fold = 0.5)
hip3_t <- tcTable(hip3tca)
#19k to 13934 to 1300
length(hip3_t)
#figure out how many clusters
# using timeclust2 with clusGap
clusGap(hip3_t,
        FUNcluster = timeclust2,
        algo = 'cm',
        K.max = 20,
        B = 20)
##11 clusters suspected
chip3tca <- timeclust(hip3tca, algo = "cm", k = 11, standardize = TRUE)
XXXchip3_px <-timeclustplot(chip3tca, value = "z-score(PRKM)", cols = 2)
hip3cxxx<-clustResults(chip3tca)

#Can see where ou significant genes fall
memb_matrix <- hip3cxxx@membership
mem_col <- as.data.frame(colnames(memb_matrix)[max.col(memb_matrix,ties.method="first")])
rownames(mem_col) <- rownames(memb_matrix)
mem_col_subset <-  mem_col[rownames(hip34resSig) %in% rownames(mem_col), ]
mem_col[rownames(mem_col) %in% rownames(hip34resSig), ]
cluster_hip34resSiglog <- subset(mem_col, rownames(mem_col) %in% rownames(hip34resSig))
colnames(cluster_hip34resSiglog) <- "membership"
cluster_hip34resSiglog$membership <- as.factor(cluster_hip34resSiglog$membership)
cluster_hip34resSiglog$membership <- factor(cluster_hip34resSiglog$membership, levels = c("1","2","3","4","5","6","7","8","9","10","11"))
str(cluster_hip34resSiglog$membership)
ggplot(cluster_hip34resSiglog, aes(x=membership))+
  geom_bar()+
  theme_bw()
#Make the figure look better
dribble <- cluster_hip34resSiglog %>%
  add_column(UpDown = NA)
dribble <- as.data.frame(dribble)
for(i in 1:length(cluster_hip34resSiglog$membership)) {
  xcc <- rownames(dribble)
  xcc <- xcc[i]
  dcc5 <-hip34resSig[rownames(hip34resSig) %in% xcc, ]
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
dribble$membership <- factor(dribble$membership, levels = c("1","2","3","4","5","6","7","8","9","10","11"))
ggplot(dribble, aes(x=UpDown))+
  geom_bar(aes(fill = membership), width = 0.5)+ 
  scale_y_continuous(expand = c(0,0))+
  theme_bw()

## STEP 6: Plot DAVID results for both cortex and hippocampus
hipcorDAVID_workable2 <- read_delim("../data/Hippocampus/DifferentialExp/combined_cortex_hippDAVID.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
hipcorDAVID_workable2 <- subset(hipcorDAVID_workable2, Category == "KEGG_PATHWAY")
hipcorDAVID_workable3 <- subset(hipcorDAVID_workable2, PValue < .05)
delme <- subset(hipcorDAVID_workable3, Organ == "Hippocampus")
ggplot(hipcorDAVID_workable2, aes(Term,Count,,fill = Direction))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0), breaks = seq(1,20,1))+
  ylab("Gene Ammount")+
  theme_bw()
ggplot(data = hipcorDAVID_workable3, aes(x = Count,  y = reorder(Term,Count),color = Direction, size = 1-PValue, shape = Organ)) + 
  geom_point() +
  #coord_flip() +
  #scale_color_gradient(low = "red", high = "blue") +
  theme_bw()

#STEP 7: Hippocampus WGCNA
#WGCNA PART1
options(stringsAsFactors = FALSE)
hipp_wg <- DESeq2::counts(dds_hip_all_noout, normalized=TRUE)
hipp_wgt = as.data.frame(t(hipp_wg))
hipp_gsg = goodSamplesGenes(hipp_wgt, verbose = 3)
hipp_gsg$allOK
if (!hipp_gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!hipp_gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(hipp_wgt)[!hipp_gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(hipp_wgt)[!hipp_gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  hipp_wgt = hipp_wgt[hipp_gsg$goodSamples, hipp_gsg$goodGenes]
}
# check to see if there are any obvious outliers
hipp_wgt_sampleTree = hclust(dist(hipp_wgt), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(hipp_wgt_sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#### load in the trait data
control <- c(rep(1,4),rep(0,12),rep(1,5))
shrinkage <- c(rep(0,4),rep(1,8),rep(0,9))
growth <- c(rep(0,12),rep(1,4),rep(0,5))
hipp_trait_2 <- as.data.frame(cbind(control,shrinkage,growth))
rownames(hipp_trait_2) <- rownames(hipp_wgt_datExpr)
#see how traits cluster
# Re-cluster samples
hipp_wgt_sampleTree2 = hclust(dist(hipp_wgt_datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_2 = numbers2colors(hipp_trait_2, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(hipp_wgt_sampleTree2, traitColors_2,
                    groupLabels = names(hipp_trait_2),
                    main = "Sample dendrogram and trait heatmap")
save(hipp_wgt_datExpr, hipp_trait_2, file = "../data/Hippocampus/WGCNA/Hipp_2-01-dataInput.RData")

#WGCNA PART2
workingDir = "../data/Hippocampus/WGCNA/";
setwd(workingDir);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Hipp_2-01-dataInput.RData");
#Note above changes for each set
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(hipp_wgt_datExpr, powerVector = powers, verbose = 5)
#FOR 3 AND 4 ABOVE IT hipp_wgt
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
##NOW NEED TO PICK POWER THRESHOLD
#9 for removed outliers as it crosses .8 threshold and according to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html we are close
softPower = 9;
adjacency = adjacency(hipp_wgt_datExpr, power = softPower);
#Creating TOM  Topological Overlap Matrix,
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
#Clustering using TOM
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicMods)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(hipp_wgt_datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(hipp_wgt_datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#Plot the new merged set
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "../data/Hippocampus/WGCNA/Hipp_1-02-networkConstruction-stepByStep.RData")

#WGCNA PART 3
workingDir = "./WGCNA"
setwd(workingDir)
library(WGCNA)
options(stringsAsFactors = FALSE)
lnames_02 = load(file = "/Users/bill/WGCNA/Hipp_2-01-dataInput.RData")
lnames_02 = load(file = "/Users/bill/WGCNA/Hipp_2-02-networkConstruction-stepByStep.RData")
nGenes = ncol(hipp_wgt_datExpr)
nSamples = nrow(hipp_wgt_datExpr);
MEs0 = moduleEigengenes(hipp_wgt_datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, hipp_trait_2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(hipp_trait_2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = rev(blueWhiteRed(50)),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               bg.lab.y = NULL,
               main = paste("Module-trait relationships"))
growth_trait = as.data.frame(hipp_trait_2$growth)
names(growth_trait) = "growth_trait"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(hipp_wgt_datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(hipp_wgt_datExpr, stage4_trait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(growth_trait), sep="");
names(GSPvalue) = paste("p.GS.", names(growth_trait), sep="");
module = "grey"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for stage4",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# Create the starting data frame
geneInfo0 = data.frame(Genes = names(hipp_wgt_datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
geneInfo0$Genes
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, growth_trait, use = "p")));
modOrder
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.growth_trait));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "../data/Hippocampus/WGCNA/hipp_02_geneInfo.csv")


#Create cytoscape object
modules = c("grey");
# Select module probes
probes = names(hipp_wgt_datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

#make WGCNA DAVID enrichment plot for figure 5C
wgcna_david_grey <-read_delim("/Users/bill/Dehnels_Seasonal_RNAseq/data/Hippocampus/WGCNA/KEGG_grey.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
wgcna_david_grey <- subset(wgcna_david_grey, Category == "KEGG_PATHWAY")
wgcna_david_grey <- subset(wgcna_david_grey, PValue < .05)
ggplot(wgcna_david_grey, aes(reorder(Term,-Count),Count,,fill = PValue))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()



