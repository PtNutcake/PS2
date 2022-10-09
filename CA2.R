library(rtracklayer)
library(AnnotationHub)
library(GenomicRanges)
library(biomaRt)

# import the human transcripts for GRCh38.
ahub <- AnnotationHub()
ahub_human <- query(ahub,c("Granges","Homo sapiens","GRCh38"))
gtf <- ahub_human[["AH105317"]]

# Create a Granges object for the promoters of all protein-coding transcripts
protein_coding <- subset(gtf, gtf$transcript_biotype == "protein_coding")

# defined for this problem set as 1500 bp upstream of the TSS and 500 bp downstream of the TSS.
TSS <- subset(protein_coding, protein_coding$type == "five_prime_utr")
promoter <- promoters(TSS, upstream = 1500, downstream = 500)
promoter

# Create a Granges object of all CpG islands for the human genome.
ahub_CpG <- query(ahub,c("Granges","Homo sapiens","CpG Islands"))
ahub_CpG
ahub_hg19CpG <- ahub_CpG[["AH5086"]]

# Transform hg19 version into hg38 version
#BiocManager::install('liftOver')
library(liftOver)
chain <- import.chain("hg19ToHg38.over.chain")
ahub_hg38CpG <- liftOver(ahub_hg19CpG, chain)
ahub_hg38CpG
ahub_hg38CpG <- unlist(ahub_hg38CpG)

# Calculate the fraction of CpG island annotations that overlap a promoter.
seqlevelsStyle(promoter) <- 'UCSC'
findOverlaps(promoter,ahub_hg38CpG)
overlap <- 28378/148358
overlap

# Plot the length distributions for CpG islands that do and do not overlap a promoter
non_overlap <- 1-overlap
label <- c("overlap","non overlap")
barplot(c(overlap, non_overlap), main = "Length Distribution", names.arg = label,ylab = "Fraction")
