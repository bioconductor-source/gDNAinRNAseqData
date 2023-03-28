.LiYu22subsetBAMfilesFromUrl <- function(baseUrl) {

  bamfiles <- .getSubDirs(baseUrl)
  bamfiles <- bamfiles[grep(".bam$", bamfiles)]
  gDNAlevel <- as.integer(gsub("s[0-9]+gDNA", "", gsub(".bam", "", bamfiles)))
  accessions <- paste0("HRR5896", gsub("s", "", gsub("gDNA[0-9]+.bam", "", bamfiles)))
  ord <- order(gDNAlevel)
  bamfiles <- bamfiles[ord]
  gDNAlevel <- gDNAlevel[ord]
  accessions <- accessions[ord]
  metadf <- data.frame(title=sprintf("RNA-seq data BAM file subset of %s contaminated with %d%% gDNA", accessions, gDNAlevel),
                       species="Homo sapiens",
                       taxonomyId=9606,
                       genome="hg38",
                       sourceUrl="https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA007961",
                       sourceVersion="2022-08-05",
                       description=sprintf("BAM file subset of %s derived from the RNA-seq data by Li et al., BMC Genomics, 23:554, 2022, contaminated with %d%% gDNA", accessions, gDNAlevel),
                       rDataPath=paste(bamfiles, gsub(".bam", ".bai", bamfiles), sep=":"))
  metadf
}

makeMetadata_LiYu22subsetBAMfiles <- function()
{
  biocver <- as.character(BiocManager::version())
  baseUrl <- "https://functionalgenomics.upf.edu/experimenthub/gdnainrnaseqdata/LiYu22subsetBAMfiles/"
  meta <- .LiYu22subsetBAMfilesFromUrl(baseUrl)
  n <- nrow(meta)
  data.frame(
    Title=meta$title,
    Description=meta$description,
    BiocVersion=rep(biocver, n),
    Genome=meta$genome,
    SourceType=rep("BAM", n),
    SourceUrl=meta$sourceUrl,
    SourceVersion=meta$sourceVersion,
    Species=meta$species,
    TaxonomyId=meta$taxonomyId,
    Coordinate_1_based=rep(TRUE, n),
    DataProvider=rep("NGDC", n),
    Maintainer=rep("Robert Castelo <robert.castelo@upf.edu>", n),
    RDataClass=rep("BamFile", n),
    DispatchClass=rep("BamFile", n),
    Location_Prefix=rep(baseUrl, n),
    RDataPath=meta$rDataPath,
    Tags=rep(paste(c("RNA-seq", "gDNA", "BAM", "gDNAinRNAseqData"), collapse=":"), n))
}

library(XML)
library(RCurl)
library(GenomeInfoDb)

source("../../R/utils.R") ## for .getSubDirs()

metadata <- makeMetadata_LiYu22subsetBAMfiles()
write.csv(metadata, file="../extdata/metadata_LiYu22subsetBAMfiles.csv", row.names=FALSE)
