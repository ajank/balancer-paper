# Idea:
# -----
# Look up bi-allelic frequency in regions where CNVs have been called
# CNV calls:
#     CNVnator
# Bi-allelic freuquency:
#     SNV calls from freebayes -> AO and RO fields
#     SNV calls from pileup + GATK -> AD field
# Algorithmic idea:
#     For each CNV call, look up region + surrounding in SNV VCF (tabix)
#     and plot bi-allelic freq. for both samples. Mark CNV prediction
#     with lines. If possible include total read-depth signal, or even
#     GC-normalized read depth signal
#
# State:
# ------
# It is apperently quite difficult to parse the SNV VCF file. Problems
# are multi-allelic sites. I have given up in the interests of more
# urgent analyses.
#

getAllelicReadCounts <- function(snv_tabix, gregion, genome="dm6") {
  
  ### Helper function
  # strsplit by the 'text' by the character 'splitter' 
  # and select k-th field
  selectField <- function(text, splitter=":", k=1) { 
    unlist(strsplit(text, splitter))[k] 
  }
  
  ### Helper function
  # Filter only bi-allelic SNPs from the VCF file
  filter.isSnp <- function(x) {
    refSnp <- nchar(ref(x)) == 1L
    a <- alt(x)
    altSnp <- elementLengths(a) == 1L
    ai <- unlist(a[altSnp])    # all length 1, so unlisting is 1:1 map
    altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
    refSnp & altSnp
  }
  ### End of helper functions
  
  ### Filter the big VCF down to bi-all. SNPs in the region
  #
  require("VariantAnnotation")
  tmp.vcf <- tempfile()
  filterVcf(snv_tabix,
            genome = genome,
            destination = tmp.vcf,
            param = ScanVcfParam(which=gregion),
            filters = FilterRules(list(isSnp=filter.isSnp)))
  
  ### Read the filtered VCF and get Allelic read depth information
  vcf_text = readLines(tmp.vcf)
  x = read.table(textConnection(vcf_text[!grepl("^##", vcf_text)]), 
                 stringsAsFactors=F, 
                 header=T, 
                 comment.char="")
  # clean up file
  file.remove(tmp.vcf)
  
  # if there is at least one entry and one sample
  if (dim(x)[1]>0 && length(colnames(x))>=10) 
  {
    samples = colnames(x)[10:length(colnames(x))]
    
    ### AD field - first try
    if (grepl(":?AD:?", x[1,9])) {
      cat("to do: Using AD field in VCF\n");
      return (data.frame())
    }
    
    ### RO and AO fields
    else if (grepl(":?RO:?", x[1,9]) && grepl(":?AO:?", x[1,9]))
    {
      message("Using RO/AO fields in VCF\n");
      pos_RO = match("RO", unlist(strsplit(x[1,9], ":", fixed=T)))
      pos_AO = match("AO", unlist(strsplit(x[1,9], ":", fixed=T)))
      data = FALSE
      for (s in samples) {
        df = data.frame(chrom  = x$X.CHROM,
                        pos    = x$POS,
                        ref    = as.vector(sapply(x[,s], FUN=selectField, k=pos_RO)),
                        alt    = as.vector(sapply(x[,s], FUN=selectField, k=pos_AO)),
                        sample = s,
                        stringsAsFactors = F)
        df = df[!grepl(".", df$ref, fixed=T) & !grepl(".", df$alt, fixed=T), ]
        df$ref = as.numeric(df$ref)
        df$alt = as.numeric(df$alt)
        if (typeof(data) == typeof(F)) data = df
        else data = rbind(data, df)
      }
      return (data)
    }
    
    ### catch all
    else
    {
      warnings("No suitable tag (AD, RO,...) found!\n");
      data.frame()
    }

  } else { 
    # if there is no sample or less than 1 entry
    data.frame()
  } # if (length(colnames(x))>=10)
}# function