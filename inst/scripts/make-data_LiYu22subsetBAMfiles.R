## Here we describe how we generated the subset of RNA-seq BAM files containing
## different levels of genomic DNA (gDNA) contamination, from the RNA-seq data
## published by:
##
## Li, X., Zhang, P., and Yu. Y. Gene expressed at low levels raise false
## discovery rates in RNA samples contaminated with genomic DNA. BMC Genomics,
## 23:554, 2022.
##
## and deposited at https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA007961
##
## The subset of the data generated here corresponds to BAM files containing
## about 100,000 alignments sampled uniformly at random for the RNA-seq
## experiments produced from total RNA libraries mixed with different
## concentrations of gDNA, concretely 0% (no contamination), 1% and
## 10%; see Fig. 2 from Li et al. (2022).
##
## 1. We downloaded the paired-end 2x50 FASTQ files from the previous URL for
## the samples with gDNA concentrations 0%, 1% and 10%, which have the
## following accession identifiers:
##
## ACCESSION  gDNA
## HRR589632    0%
## HRR589633    0%
## HRR589634    0%
## HRR589626    1%
## HRR589627    1%
## HRR589628    1%
## HRR589623   10%
## HRR589624   10%
## HRR589625   10%
##
## 2. We downloaded the hg38 version of the human genome that includes human
## decoy sequences from hs38d1 (GCA_000786075.2) in addition to the sequences
## of the chromosomes, mitochondrial genome, unlocalized scaffolds, and
## unplaced scaffolds:
##
## $ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
## $ gzip -d GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
## $ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai
##
## 3. We downloaded the primary assembly GENCODE annotations of the same
## version as the current release (3.16) of TxDb.Hsapiens.UCSC.hg38.knownGene
## (GENCODE v41):
##
## $ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz
## 
## 4. We built a hg38 genome index for STAR as follows:
##
## $ STAR --runThreadN 10 \
##        --runMode genomeGenerate \
##        --genomeDir $(pwd)/ \
##        --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
##        --sjdbGTFfile gencode.v41.primary_assembly.annotation.gtf \
##        --sjdbOverhang 49
##
## 5. We aligned the FASTA files to the previous STAR index as follows:
##
## $ STAR --genomeDir $GENOMEINDEX \
##       --readFilesIn $RAWDATADIR/$FASTQL $RAWDATADIR/$FASTQR \
##       --readFilesCommand zcat --runThreadN $NTHREADS \
##       --peOverlapNbasesMin 10 --outSAMtype BAM Unsorted \
##       --outReadsUnmapped Fastx --outFileNamePrefix BAM/$PREFIX
##
## where $GENOMEINDEX, $RAWDATADIR, $FASTQL, $FASTQR, $PREFIX and $NTHREADS are
## bash variables storing the STAR genome index, the directory with the raw
## FASTQ files, the FASTQ file for the left/first/R1 reads, the FASTQ file for
## the right/last/R2 reads, the prefix for the FASTQ files and the number of
## execution threads to use by STAR, respectively.
##
## 6. We processed with Picard the initial BAM files output by STAR as follows:
##
## $ mkdir -p BAM/sortedbyqname
##
## $ java -jar picard.jar AddOrReplaceReadGroups \
##        -INPUT BAM/${PREFIX}Aligned.out.bam \
##        -OUTPUT BAM/sortedbyqname/$PREFIX.bam \
##        -SORT_ORDER queryname \
##        -RGID $PREFIX \
##        -RGSM $PREFIX \
##        -RGLB $PREFIX \
##        -RGPL Illumina \
##        -RGPU RUN1
##
## $ mkdir -p BAM/markduplicates
##
## $ java -jar picard.jar MarkDuplicates \
##        -INPUT BAM/sortedbyqname/$PREFIX.bam \
##        -OUTPUT BAM/markduplicates/$PREFIX.bam \
##        -METRICS_FILE BAM/markduplicates/$PREFIX.metrics.txt
##
## $ java -jar picard.jar SortSam \
##        -INPUT BAM/markduplicates/$PREFIX.bam \
##        -OUTPUT BAM/$PREFIX.bam \
##        -SORT_ORDER coordinate \
##        -CREATE_INDEX true
##
## 7. We sampled about 100,000 alignments from the resulting BAM files in two
## steps. First we used the following R script that first sample the read
## identifiers (query names -QNAME- in the BAM file) from the BAM files:
##
library(GenomicRanges)
library(GenomicFiles)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(BiocParallel)

bams <- list.files(".", pattern="*.bam$", full.names=TRUE)
dat <- read.csv("phenoData.csv")
dat <- dat[dat$LIBPREP == "RibosomalDepletion", ]
mask <- gsub(".bam", "", basename(bams)) %in% dat$ACCESSION[dat$gDNA %in% c(0, 1, 10)]
bfl <- BamFileList(bams[mask], yieldSize=1000000L)
sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                       isProperPair=TRUE,
                       isSecondaryAlignment=FALSE,
                       isDuplicate=FALSE,
                       isNotPassingQualityControls=FALSE)
yield <- function(x) readGAlignmentPairs(x, param=ScanBamParam(flag=sbflags, what="qname"))
map <- identity
n <- 100000
cat("Sampling...\n")
gals <- bplapply(bfl, reduceByYield, YIELD=yield, MAP=map,
                 REDUCE=REDUCEsampler(n, TRUE), parallel=FALSE,
                 BPPARAM=MulticoreParam(workers=9))
cat("Processing...\n")
qnames <- lapply(gals, function(x) unique(mcols(first(x))$qname))
cat("Writing...\n")
for (bam in names(qnames)) {
    txt <- gsub("bam", "txt", bam)
    cat(sprintf("exporting to %s\n", txt))
    write.table(data.frame(QNAME=qnames[[bam]]), txt, row.names=FALSE,
                col.names=FALSE, quote=FALSE)
}
## where the CSV file 'phenoData.csv' is derived from the metadata file
## HRA001834.xlsx available at https://ngdc.cncb.ac.cn/gsa-human/browse/HRA001834
## Finally, we use the following shell script to subset the BAM files with
## SAMtools that names the resulting files using the last two digits from the
## accession and the gDNA contamination level:
##
## $ for i in *.txt ; do {
##     p=${i%%.txt}
##     gDNA=`cat ../RawData/phenoData.csv | grep RibosomalDepletion | sed 's/\"//g;s/,/\t/g' | cut -f 1,3 | egrep '[[:space:]]0$|[[:space:]]1$|[[:space:]]10$' | grep $p | cut -f 2`
##     p2=`echo $p | sed 's/HRR5896/s/'`gDNA$gDNA
##     echo $p
##     samtools view -N $i -h -b -o $p2.bam ../Align/BAMsRD/$p.bam
##     samtools sort $p2.bam -o $p2.sorted.bam
##     mv $p2.sorted.bam $p2.bam
##     samtools index $p2.bam
##     mv $p2.bam.bai $p2.bai
## } done
