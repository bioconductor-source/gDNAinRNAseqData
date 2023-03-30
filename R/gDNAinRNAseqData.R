#' RNA-seq data with different levels of gDNA contamination
#'
#' This package provides access to RNA-seq BAM files containing different
#' levels of genomic DNA (gDNA) contamination.\cr\cr
#'
#' Currently, this package allows one to download a subset of the data
#' published in:
#'
#' Li, X., Zhang, P., and Yu. Y. Gene expressed at low levels raise false
#' discovery rates in RNA samples contaminated with genomic DNA. BMC Genomics,
#' 23:554, 2022.
#'
#' The subset of the data acccessible through this package corresponds to
#' BAM files containing about 100,000 alignments sampled uniformly at random
#' for the RNA-seq experiments produced from total RNA libraries mixed with
#' different concentrations of gDNA, concretely 0\% (no contamination), 1\% and
#' 10\%; see Fig. 2 from Li et al. (2022).
#'
#' @references
#' Li, X., Zhang, P., and Yu. Y. Gene expressed at low levels raise false
#' discovery rates in RNA samples contaminated with genomic DNA. BMC Genomics,
#' 23:554, 2022.
#'
#' @docType package
#' @name gDNAinRNAseqData-package
#' @aliases gDNAinRNAseqData-package
#' @aliases gDNAinRNAseqData
#' @keywords package
NULL

#' @describeIn gDNAinRNAseqData-package downloads the BAM files from the
#' RNA-seq data through the ExperimentHub, and returns the path in the
#' filesystem where the BAM files are stored.
#'
#' @param path (Default=`tempdir()`) Filesystem path where to store the
#' BAM files.
#'
#' @param offline (Default=`FALSE`) If there is no internet connection, but
#' the data has been previously downloaded, setting `offline=TRUE` allows one
#' to retrive the data from the ExperimentHub cache.
#'
#' @return `LiYu22subsetBAMfiles()` returns a string character vector of
#' filesystem paths to the downloaded BAM files.
#'
#' @examples
#'
#' ## for LiYu2subsetBAMfiles()
#' bamfiles <- LiYu22subsetBAMfiles()
#' bamfiles
#'
#' @importFrom BiocGenerics subset path basename
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom Rsamtools index
#' 
#' @export
LiYu22subsetBAMfiles <- function(path=tempdir(), offline=FALSE) {
    eh <- NULL
    tryCatch({
        eh <- ExperimentHub(localHub=offline)
    }, error = function(e) {
        message("Could not load the ExperimentHub: ", e)
        return(character(0))
    })
    if (is.null(eh)) {
        message("Could not connect to the ExperimentHub.")
        return(character(0))
    }

    preparerclass <- rdataclass <- NULL ## to avoid NOTE at R CMD check
    ehres <- subset(eh, preparerclass=="gDNAinRNAseqData" & rdataclass=="BamFile")
    bamfiles <- basename(gsub(", .+$", "", ehres$rdatapath))
    baifiles <- gsub(".bam", ".bai", bamfiles)
    for (i in seq_along(names(ehres))) {
        ehbam <- ehres[[i]]
        file.copy(path(ehbam), file.path(path, bamfiles[i]))
        file.copy(index(ehbam), file.path(path, baifiles[i]))
    }

    file.path(path, bamfiles)
}

#' @describeIn gDNAinRNAseqData-package retrieves phenotypic data from the
#' BAM files downloaded with `LiYu22subsetBAMfiles()`.
#'
#' @param bamfiles full filesystem paths to where the BAM files were downloaded
#' with `LiYu22subsetBAMfiles()`.
#'
#' @return `LiYu22phenoData()` returns a `data.frame` object with the gDNA
#' contamination levels for the BAM files specified in the `bamfiles` parameter,
#' according to the publication by Li et al. (2022).
#'
#' @examples
#'
#' ## for LiYu22phenoData()
#' bamfiles <- LiYu22subsetBAMfiles()
#' LiYu22phenoData(bamfiles)
#'
#' @importFrom BiocGenerics basename
#'
#' @export
LiYu22phenoData <- function(bamfiles) {
    fexists <- vapply(bamfiles, FUN=file.exists, FUN.VALUE=logical(1))
    if (any(!fexists))
        stop(sprintf("File(s) %s cannot be found.",
                     paste(bamfiles[!fexists], collapse=", ")))
    gDNAlevel <- as.integer(gsub("s[0-9]+gDNA", "",
                                 gsub(".bam", "", basename(bamfiles))))
    data.frame(gDNA=gDNAlevel, row.names=gsub(".bam", "", basename(bamfiles)))
}
