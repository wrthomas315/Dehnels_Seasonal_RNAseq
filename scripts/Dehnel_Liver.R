## Here is the R code for processing the temporal transcriptomics of Dehnel's phenomenon liver
#Reads were preprocessed on Noctillio server, where they were trimmed using fastp and aligned/quantified using Kalisto
#First I want to get the quantification files onto R and convert transcript abundance to gene abundance

###First set up all your libraries
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

#now import quant files previously created using kallisto
cs_liv_samples <- read.table("../data/Liver/SampleList", header = T)
cs_liv_samples
cs_liv_files <-file.path("../data/Liver/TranscriptAbundances/", cs_liv_samples$Sample_name, "abundance.tsv")
names(cs_liv_files) <- paste0("sample_", cs_liv_samples$Sample_name)
all(file.exists(cs_liv_files))

#write to tables and can save with hashed out code below, but note will overwrite
cs_liv.count.tsv <- tximport(cs_liv_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_liv.tpm.tsv <- tximport(cs_liv_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
#write.table(cs_liv.tpm.tsv$abundance, "../data/Liver/GeneCounts/GeneCounts.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

## STEP2: Graph Liver mass for panel A
SampleData_Rreadable <- read_csv("~/Dehnels_Seasonal_RNAseq/data/Liver/SampleData_Rreadable.csv")
livermass <- SampleData_Rreadable$Liver_Mass
bodymass <- SampleData_Rreadable$Body_Mass
stagemass <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
livermean <- rep(mean(livermass,na.rm=TRUE))
livermass_matrix <- as.data.frame(cbind(livermass,bodymass))
livermass_matrix$norm <- (livermass_matrix$livermass/livermass_matrix$bodymass)
livermass_matrix$rel_mass <- c(rep(median(livermass_matrix$norm[1:5],na.rm=TRUE),24))
indiv_livmat <- as.data.frame(cbind(livermass_matrix,stagemass))

#Stage2vs4 shows significant size change
t.test(indiv_livmat$livermass[15:19],indiv_livmat$livermass[6:9],)
p.adjust(0.000248, method = "bonferroni", n = 4 )
#plot
livmassplot <-ggplot(indiv_livmat, aes(x = stagemass, y = norm)) +
  geom_boxplot(size=.5)+
  xlab("Stage")+
  ylab("Normalized Liver Mass ( g / g BodyMass)")+
  theme_bw()
ggsave("../data/Liver/LiverMass.png", livmassplot)

#Without normalization
ggplot(livermass_df,aes(x=liverstage,y=livermass))+
  geom_boxplot()+
  scale_y_continuous(name="LiverMass", limits=c(0.5,1.4))+
  theme_bw()

## STEP3: Normalization and QC checks
#Begin by setting up design matrix
colnames(cs_liv.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
liv_stages <- factor(c(cs_liv_samples$Run))
liv_organs <- factor(c(cs_liv_samples$Condition))
liv_full1 <- factor(c(cs_liv_samples$Run))
liv_stages_organ_frame <-cbind(as.data.frame(liv_stages),as.data.frame(liv_organs),as.data.frame(liv_full1))

#check to make sure RNA Integrity isn't wacky for alignment purposes
vis_liv_samples <- SampleData_Rreadable
ggplot(vis_liv_samples,aes(x=Liver_RIN,y=Liver_Reads_prefilt, color= Stage))+
  geom_point()+
  theme_bw()
ggplot(vis_liv_samples,aes(x=Liver_RIN,y=Liver_Reads_postfilt, color= Stage))+
  geom_point()+
  theme_bw()
ggplot(vis_liv_samples,aes(x=Liver_RIN,y=Liver_Pseudo_aligned_Percent, color= Stage))+
  geom_point()+
  theme_bw()
#minimal effect on sequencing output, ~3-5% on mapping percent, but spread across stages so not nervous about diferential effect

#Now to normalize our reads
dds_liv_all <- DESeqDataSetFromMatrix(round(cs_liv.count.tsv$counts), DataFrame(liv_stages_organ_frame), ~ liv_full1)
mcols(dds_liv_all) <- cbind(mcols(dds_liv_all), row.names(cs_liv.count.tsv$counts))
rownames(dds_liv_all) <- row.names(cs_liv.count.tsv$counts)
dds_liv_all <- DESeq(dds_liv_all)
#And look at pcas and heatmaps of counts
vst_dds_liv_all <- vst(dds_liv_all)
pcaData_liv_all<- plotPCA(vst_dds_liv_all,intgroup=c("liv_stages","liv_organs"), ntop=300, returnData=TRUE)
ggplot(pcaData_liv_all, aes(x = PC1, y = PC2, color = factor(liv_stages))) +
  geom_point(size=2)+
  theme_bw()
###LiverHeatmap###
livsampleDists <- dist(t(assay(vst_dds_liv_all)))
livsampleDistMatrix <- as.matrix(livsampleDists)
colnames(livsampleDistMatrix) <- NULL
#make the heatmap
pheatmap(livsampleDistMatrix, clustering_distance_rows=livsampleDists,
         clustering_distance_cols = livsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
#make a heatmap of genes now
livtopVarGenes <- head( order( rowVars( assay(vst_dds_liv_all) ), decreasing=TRUE ), 15 )
heatmap.2( assay(vst_dds_liv_all)[ livtopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

## STEP 4 : Differential Expression
#LivStage2vs4
lfc <- 1.58
liv24res <- results(dds_liv_all, contrast = c("liv_full1","Stage4","Stage2"))
liv24resX <-  results(dds_liv_all, contrast = c("liv_full1","Stage4","Stage2"))
liv24resSig <- subset(liv24res,liv24res$padj<.05)
liv24resSigLog <- subset(liv24resSig,abs(liv24resSig$log2FoldChange)>=lfc)
liv24up <- subset(liv24resSigLog,(liv24resSigLog$log2FoldChange)>=0)
liv24down <- subset(liv24resSigLog,(liv24resSigLog$log2FoldChange)<=0)
#And plot with a volcano plot
#First change any read to a triangle if off the plot (>abs(6))
for (i in 1:length(liv24resX$padj)) {
  if  (liv24resX$padj[i]<1e-15 & !is.na (liv24resX$padj[i])) {
    liv24resX$padj[i] <- 1e-15
  }
  if (liv24resX$log2FoldChange[i]>6 & !is.na (liv24resX$log2FoldChange[i])) {
    liv24resX$log2FoldChange[i] <- 6
  }
  if (liv24resX$log2FoldChange[i]< -6 & !is.na (liv24resX$log2FoldChange[i])) {
    liv24resX$log2FoldChange[i] <- -6
  }
}
#create custom key-value pairs for different cell-types
#this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(liv24resX$log2FoldChange) == 6, 17,
  ifelse(liv24resX$padj==1e-15, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
#And also change color for high and low effect, and not significant
keyvals <- ifelse(
  liv24resX$padj > 0.05, 'grey',
  ifelse(liv24resX$log2FoldChange <= -1.58, 'red',
         ifelse(liv24resX$log2FoldChange >= 1.58, 'blue',
                ifelse(liv24resX$log2FoldChange >= 0, 'lightblue',
                       'pink'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'pink'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'lightblue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
#Plot with genes identified in analysis of interet
EnhancedVolcano(liv24resX,
                lab = rownames(liv24resX),
                xlim=c(-6 ,6),
                ylim=c(0,15),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = c('PPARA','FOXO3','PPARD','RXRG','FOXO1','RXRA','PER1'),
                pCutoff = .05,
                FCcutoff = 1.58,
                colCustom = keyvals,
                legendPosition = 'right',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.75)

## STEP 6 : Liver Time Seq analysis
#Timeseq cortex with LFC .5849 (*1.5 change)
liv_TS <- c(rep("Stg1",5),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
liv_TS_frame <- data.frame(sampleID = 1:24, group = c(1, 1, 1,1,1,2, 2, 2,2, 3, 3, 3,3,3,4, 4, 4, 4,4,5, 5, 5, 5,5),
                               timepoint = liv_TSredo)
liv_TSgf <- data.frame(chr = c(rep('chr1', 15296), rep('chr2', 2000), rep('chr4', 2000)),
                           start = rep(100, 19296),
                           end = rep(500, 19296),
                           id = row.names(cs_liv.tpm.tsv$abundance))
library(TCseq)
livtca <- TCA(design = liv_TS_frame, counts = round(DESeq2::counts(dds_liv_all, normalized=TRUE)), liv_TSgf)
livtca <- DBanalysis(livtca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)
livtca <- timecourseTable(livtca, value = "expression",  lib.norm = FALSE, filter = TRUE,abs.fold = 0.5)
liv_x <- tcTable(livtca)
#figure out how many clusters
# using timeclust2 with clusGap
clusGap(liv_x,
        FUNcluster = timeclust2,
        algo = 'cm',
        K.max = 20,
        B = 20)
##14 clusters for liver
#but them into these 14 clusters and plot
lliv_x <- timeclust(livtca, algo = "cm", k = 14, standardize = TRUE)
lliv_px2 <-timeclustplot(lliv_x, value = "z-score(PRKM)", cols = 2,cl.color = "gray50",membership.color= hcl.colors(30, "Viridis",rev = TRUE))
#and save the results
lliv_cxxx<-clustResults(lliv_redo)

#### where do the significant genes reside from Stage 2 vs Stage 4 DESeq
memb_matrix <- lliv_cxxx@membership
mem_col <- as.data.frame(colnames(memb_matrix)[max.col(memb_matrix,ties.method="first")])
rownames(mem_col) <- rownames(memb_matrix)
mem_col_subset <-  mem_col[rownames(liv24resSigLog) %in% rownames(mem_col), ]
mem_col[rownames(mem_col) %in% rownames(liv24resSigLog), ]
cluster_liv24resSigLog <- subset(mem_col, rownames(mem_col) %in% rownames(liv24resSigLog))
colnames(cluster_liv24resSigLog) <- "membership"
cluster_liv24resSigLog$membership <- as.factor(cluster_liv24resSigLog$membership)
cluster_liv24resSigLog$membership <- factor(cluster_liv24resSigLog$membership, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14"))
length(which(cluster_liv24resSigLog$membership == "8"))
ggplot(cluster_liv24resSigLog, aes(x=membership))+
  geom_bar()+
  theme_bw()
#Most in 8 and 12, make graph look better belw
dribble <- cluster_liv24resSigLog %>%
  add_column(UpDown = NA)
dribble <- as.data.frame(dribble)
dribble
for(i in 1:length(cluster_liv24resSigLog$membership)) {
  xcc <- rownames(dribble)
  xcc <- xcc[i]
  dcc5 <-liv24resSigLog[rownames(liv24resSigLog) %in% xcc, ]
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
dribble$membership <- factor(dribble$membership, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14"))
colorblind14  <- c("#5378DE","#71BFB6","#5014E2","#823EC3","#F8AD09","#C28455","#882255","#AA4499","#CC6677","#DDCC77","#88CCEE","#44AA99","#117733","#332288")
ggplot(dribble, aes(x=UpDown))+
  geom_bar(aes(fill = membership), width = 0.5)+ 
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual( values = colorblind14)+
  theme_bw()

#@ STEP 7: Make Figure 2C of cycling TFs
fig4liv_mat <- DESeq2::counts(dds_liv_all, normalized=TRUE)
fig4_foxo1 <-fig4liv_mat[rownames(fig4liv_mat) %in% "FOXO1", ]
fig4_foxo3 <-fig4liv_mat[rownames(fig4liv_mat) %in% "FOXO3", ]
fig4_ppara <-fig4liv_mat[rownames(fig4liv_mat) %in% "PPARA", ]
fig4_pparb <-fig4liv_mat[rownames(fig4liv_mat) %in% "PPARD", ]
fig4_rxra <-fig4liv_mat[rownames(fig4liv_mat) %in% "RXRA", ]
fig4_rxrg <-fig4liv_mat[rownames(fig4liv_mat) %in% "RXRG", ]
fig4_per1 <-fig4liv_mat[rownames(fig4liv_mat) %in% "PER1", ]
splliv_stage  <- rep(c(c(rep("Stg1",5),rep("Stg2",4),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))),7)
livspl <- rbind(fig4_foxo1,fig4_foxo3,fig4_ppara,fig4_pparb,fig4_rxra,fig4_rxrg,fig4_per1)
#get zcore for each gene
z <- matrix(0,nrow(livspl),ncol(livspl))
for (i in 1:nrow(livspl)) {
  for (j in 1:ncol(livspl)) {
    z[i,j] <- (livspl[i,j]-mean(livspl[i,]))/sd(livspl[i,])
  }
}
zz  <- matrix(0,nrow(z),5)
for (i in 1:nrow(z)) {
  zz[i,1] <- sum(z[i,1],z[i,2],z[i,3],z[i,4],z[i,5])/5
  zz[i,2] <- sum(z[i,6],z[i,7],z[i,8],z[i,9])/4
  zz[i,3] <- sum(z[i,10],z[i,11],z[i,12],z[i,13],z[i,14])/5
  zz[i,4] <- sum(z[i,15],z[i,16],z[i,17],z[i,18],z[i,19])/5
  zz[i,5] <- sum(z[i,20],z[i,21],z[i,22],z[i,23],z[i,24])/5
}
real_livspl <- as.matrix(rbind(as.matrix(zz[1,]),as.matrix(zz[2,]),as.matrix(zz[3,]),as.matrix(zz[4,]),as.matrix(zz[5,]),as.matrix(zz[6,]),as.matrix(zz[7,])))
close <-as.matrix(c(rep("FOXO1",5),rep("FOXO3",5),rep("PPARA",5),rep("PPARB",5),rep("RXRA",5),rep("RXRG",5),rep("PER1",5)))
closer <- as.matrix(c(rep("DHE",5),rep("DLE",5),rep("ULE",15),rep("UHE",5),rep("DHE",5)))
vclose<-as.matrix(c(rep("Liver",35)))
soclose<- rep(c("Stg1","Stg2","Stg3","Stg4","Stg5"),7)
realreal <-  cbind(as.data.frame(real_livspl),as.data.frame(close),as.data.frame(closer),as.data.frame(soclose),as.data.frame(vclose))
colnames(realreal) <- c("Expr","Gene","Function","Stage","Tissue")
#Figure 2A plot
ggplot(realreal) +
  geom_point(aes(x = Stage, y =Expr,colour=Function,shape=Gene),size=5)+
  geom_line(aes(x = Stage, y =Expr,colour=Function,group=Gene),size=1)+
  #geom_line(aes(x = Stage, y =Expr,group=Gene,shape=Gene,color=Gene))+
  scale_color_manual(values=c("red","pink","blue","lightblue"))+
  scale_shape_manual(values=c(0,2,1,15,16,17,18))+
  theme_bw()

#STEP8: WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)
liv_wg <- DESeq2::counts(dds_liv_all, normalized=TRUE)
liv_wgt = as.data.frame(t(liv_wg))
liv_gsg = goodSamplesGenes(liv_wgt, verbose = 3)
liv_gsg$allOK
if (!liv_gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!liv_gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(liv_wgt)[!liv_gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(liv_wgt)[!liv_gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  liv_wgt = liv_wgt[liv_gsg$goodSamples, liv_gsg$goodGenes]
}
# check to see if there are any obvious outliers
liv_wgt_sampleTree = hclust(dist(liv_wgt), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(liv_wgt_sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Remove outliers, each stage trait
#remove outliers
# Plot a line to show the cut
abline(h = 800000, col = "red");
# Determine cluster under the line
liv_wgt_clust = cutreeStatic(liv_wgt_sampleTree, cutHeight = 800000, minSize = 10)
liv_wgt_clust
# clust 1 contains the samples we want to keep.
liv_wgt_keepSamples = (liv_wgt_clust==1)
liv_wgt_datExpr = liv_wgt[liv_wgt_keepSamples, ]
liv_wgt_nGenes = ncol(liv_wgt_datExpr)
liv_wgt_nSamples = nrow(liv_wgt_datExpr) 

#Now put in traits
control <- c(rep(1,5),rep(0,14),rep(1,5))
livershrinkage <- c(rep(0,5),rep(1,4),rep(0,15))
livergrowth <- c(rep(0,9),rep(1,10),rep(0,5))
DPshrinkage <- c(rep(0,5),rep(1,9),rep(0,10))
DPgrowth <- c(rep(0,14),rep(1,5),rep(0,5))
liv_trait <- as.data.frame(cbind(control,livershrinkage,livergrowth,DPshrinkage,DPgrowth))
rownames(liv_trait) <- rownames(liv_wgt_datExpr)
liv_wgt_sampleTree2 = hclust(dist(liv_wgt_datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(liv_trait, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(liv_wgt_sampleTree2, traitColors,
                    groupLabels = names(liv_trait),
                    main = "Sample dendrogram and trait heatmap")
save(liv_wgt_datExpr, liv_trait, file = "../data/Liver/WGCNA/Liv-01-dataInput.RData")

#WGCNA Part2
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
# Load the data saved in the first part, or save as new name if conducting in same script
lnames = load(file = "../data/Liver/WGCNA/Liv-01-dataInput.RData");
#Note above changes for each set
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(liv_wgt_datExpr, powerVector = powers, verbose = 5)
#FOR 3 AND 4 ABOVE IT liv_wgt
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
adjacency = adjacency(liv_wgt_datExpr, power = softPower);
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
MEList = moduleEigengenes(liv_wgt_datExpr, colors = dynamicColors)
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
merge = mergeCloseModules(liv_wgt_datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
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
save(MEs, moduleLabels, moduleColors, geneTree, file = "../data/Liver/WGCNA/liv_1-02-networkConstruction-stepByStep.RData")

#WGCNA Part3
#load names again if running in different window
options(stringsAsFactors = FALSE)
lnames_01 = load(file = "/Users/bill/CompleteShrew/PaperData/Liver/WGCNA/Liv-01-dataInput.RData")
lnames_01 = load(file = "/Users/bill/CompleteShrew/PaperData/Liver/WGCNA/Liv-02-networkConstruction-stepByStep.RData")
nGenes = ncol(liv_wgt_datExpr)
nSamples = nrow(liv_wgt_datExpr);
MEs0 = moduleEigengenes(liv_wgt_datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
length(MEs[1,])
moduleTraitCor = cor(MEs, liv_trait, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(liv_trait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#Figure 3A
library(viridis)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(liv_trait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = rev(blueWhiteRed(50)),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               bg.lab.y = NULL,
               main = paste("Module-trait relationships"))
viridis(50, alpha = 1, begin = 0, end = 1, direction = 1, option = "C")
#
livershrinkage = as.data.frame(liv_trait$livershrinkage)
names(livershrinkage) = "livershrinkage"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(liv_wgt_datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(liv_wgt_datExpr, livershrinkage, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(livershrinkage), sep="");
names(GSPvalue) = paste("p.GS.", names(livershrinkage), sep="");
module = "mediumpurple3"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for livershrinkage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# Create the starting data frame
geneInfo0 = data.frame(Genes = names(liv_wgt_datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
geneInfo0$Genes
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, livershrinkage, use = "p")));
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
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.livershrinkage));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "../data/Liver/WGCNA/liv_geneInfo.csv")
# Select modules
modules = c("cyan");
# Select module probes
probes = names(liv_wgt_datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-Liv_thresh", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-Liv_thresh", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);


#STEP 9: Gene set enrichment plot for Supplemental
#Note all DAVID analyses were done using their online tool
cyan_DAVID <- read_delim("../data/Liver/WGCNA/cyan_workable.tsv")
cyan_DAVID$Percent <- as.numeric(cyan_DAVID$Percent)
ggplot(data = cyan_DAVID, aes(x = Count,  y = reorder(Term,Count),color = Benjamini, size = Percent)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw()
