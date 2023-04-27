# Dehnels_Seasonal_RNAseq
Scripts and code to reproduce RNAseq analysis for looking at changes in expression through Dehnel's phenomenon; specifically in the cortex, hippocampus, and liver.

###Goals and strategy

The objective of this project is to evaluate how expression changes 
throughout the size plasticity of Dehnel’s phenomenon in *Sorex araneus*. 
We will be analyzing three different tissue types; the liver, hippocampus, 
and cortex. These regions will help to understand both size and metabolic 
changes underpinning Dehnel’s phenomenon. First, we will have to access 
the quality of our RNA-seq data, filter low quality reads and trim 
adapters, map to the transcriptome and quantify abundance, followed by 
normalization. Then we will 1) analyze differential expression between 
stages of Dehnel’s phenomenon using DESeq2, 2) characterize temporal 
patterns in expression using TCSeq, and 3) build gene correlation networks 
and identify correlation between network structure and traits. Throughout 
the analysis, we will look at resultant genes and test whether they enrich 
KEGG pathways using DAVID Functional Enrichment Tools.

###Data

RNA-seq analyses require alignment to a reference and quantification of 
reads. The genome and original unfiltered reads can be downloaded as 
described below. However, these steps could be skipped when reproducing, 
as count data has been saved in ./data/<tissue>/GeneCounts.

The reference (sorAra2; GCF_000181275.1) can be download from straight 
from NCBI, or using the code below.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/181/275/GCF_000181275.2_SorAra2.0/GCF_000181275.2_SorAra2.0_genomic.gtf.gz ./data/
gunzip ./gtf/GCF_000181275.2_SorAra2.0_genomic.gtf.gz
```

RNA-seq data from this project can also be found on NCBI Sequencing Read 
Archive, <insert project accession>. The list of samples and associated 
accession numbers can be found in the data folder. These can be downloaded 
manually, or using the getter.sh script.

###Quality control, filtering, trimming

