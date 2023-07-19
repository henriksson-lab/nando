

# https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html

################################################################################
################ Pando: Bug fix ################################################
################################################################################

# There is a bug in pandos initiate_grn.Seurat that may scramble enhancer-gene linkages.
# This is a fixed version

fixed.initiate_grn.Seurat <- function(
    object,
    regions = NULL,
    peak_assay = 'peaks',
    rna_assay = 'RNA',
    exclude_exons = TRUE
){
  
  gene_annot <- Signac::Annotation(object[[peak_assay]])  
  # Through error if NULL
  if (is.null(gene_annot)){
    stop('Please provide a gene annotation for the ChromatinAssay.')
  }
  peak_ranges <- StringToGRanges(rownames(GetAssay(object, assay=peak_assay))) 
  
  # Find candidate ranges by intersecting the supplied regions with peaks
  # Per default take all peaks as candidate regions
  if (!is.null(regions)){
    cand_olaps <- IRanges::findOverlaps(regions, peak_ranges)
    cand_ranges <- IRanges::pintersect(
      peak_ranges[subjectHits(cand_olaps)],
      regions[queryHits(cand_olaps)]
    )
  } else {
    cand_ranges <- peak_ranges
  }
  
  # Exclude exons because they are usually conserved
  if (exclude_exons){
    exon_ranges <- gene_annot[gene_annot$type=='exon', ]
    names(exon_ranges@ranges) <- NULL
    exon_ranges <- IRanges::intersect(exon_ranges, exon_ranges)
    exon_ranges <- GenomicRanges::GRanges(
      seqnames = exon_ranges@seqnames,
      ranges = exon_ranges@ranges
    )
    cand_ranges <- IRanges::setdiff(cand_ranges, exon_ranges, ignore.strand=TRUE)
  }
  
  # Match candidate ranges to peaks  -- an ugly hack added here to avoid duplication
  peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)
  #peak_matches <- subjectHits(peak_overlaps)
  peak_overlaps <- as.data.frame(peak_overlaps)
  peak_overlaps <- peak_overlaps[!duplicated(peak_overlaps$queryHits),]
  peak_matches <- peak_overlaps$subjectHits
  
  regions_obj <- new(
    Class = 'Regions',
    ranges = cand_ranges,
    peaks = peak_matches,
    motifs = NULL
  )
  
  params <- list(
    peak_assay = peak_assay,
    rna_assay = rna_assay,
    exclude_exons = exclude_exons
  )
  
  grn_obj <- new(
    Class = 'RegulatoryNetwork',
    regions = regions_obj,
    params = params
  )
  
  object <- as(object, 'SeuratPlus')
  object@grn <- grn_obj
  return(object)
}



# 
# 
# 
# ##Helper function: foreach applied over a named list
# foreach_namedlist <- function(thelist, myfunc){
#   outlist <- foreach (x = thelist, .final = function(x) setNames(x, names(thelist)),  .verbose = F) %do% {
#     #print(x)
#     myfunc(x)
#   }
#   outlist
# }


