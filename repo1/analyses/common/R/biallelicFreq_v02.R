require("plyr") # join
require("VariantAnnotation")
require("rtracklayer")
require("ggplot2")
require("assertive")

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

getAllelicReadCounts <- function(snv_tabix, gregion, genome="dm6") 
{
  # Sanity check
  is2(gregion, "GRanges")
  is2(snv_tabix, "TabixFile")
  
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
    altSnp <- elementNROWS(a) == 1L
    ai <- unlist(a[altSnp])    # all length 1, so unlisting is 1:1 map
    altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
    refSnp & altSnp
  }
  ### End of helper functions
  
  ### Filter the big VCF down to bi-all. SNPs in the region
  #
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
      # Return empty data.frame bit with colnames
      data.frame(chrom = NA, pos = NA, ref = NA, alt = NA, sample = NA)[-1,]
    }
  } else { 
    # if there is no sample or less than 1 entry
    data.frame(chrom = NA, pos = NA, ref = NA, alt = NA, sample = NA)[-1,]
  } # if (length(colnames(x))>=10)
}# function




# Take coverage tracks for several samples and plot them together
svview_get_coverage <- function(region, ...)
{
  # Get arguments and check sanity
  is2(region, "GRanges")
  tracks <- list(...)
  lapply(tracks, function(x) is2(x, "TabixFile"))
  
  # Get regions from TabixFiles
  D = NULL
  for(i in 1:length(tracks))
  {
    data = import.bedGraph(tracks[[i]], which = region)
    D = rbind(D, data.frame(x = start(data), 
                            len = width(data), 
                            y = data$score,
                            sample = names(tracks)[i]))
  }
  D
}

svview_plot_coverage <- function(region, ...)
{
  D = svview_get_coverage(region, ...)
  ggplot(D) + 
    aes(x=x, y=y, radius=len, col=sample) + 
    geom_spoke(size=2,angle=0, alpha=0.8) +
    ylab(NULL)
}


# takes the raw mapping bedGraph file
svview_get_mappability <- function(region, bedGraph, normalize = 1) 
{
  # Get arguments and check sanity
  is2(region, "GRanges")
  is2(bedGraph, "TabixFile")
  assert_is_a_number(normalize)
  
  raw = import.bedGraph(bedGraph, which = region)
  data.frame(pos = start(raw), 
                   width = width(raw), 
                   score = mcols(raw)$score / normalize)
}


svview_plot_mappability <- function(region, bedGraph, normalize = 1) 
{
  df = svview_get_mappability(region, bedGraph, normalize)
  ggplot(df) + 
    aes(x=pos, y=score, radius=width,col="a") + 
    geom_spoke(size=2,angle=0) +
    ylab(NULL) +
	scale_color_manual(values = c(a = "black"), labels = c(a = "uniqueness score"))
}


svview_biallelicfreq <- function(region, snp_vcf, mapp_bedGraph=NULL, uniqueValue)
{
  # Get arguments and check sanity
  is2(region, "GRanges")
  is2(snp_vcf, "TabixFile")

  data_snps = getAllelicReadCounts(snp_vcf, region)
  plot_snps = ggplot(data_snps)
  
  if (!is.null(mapp_bedGraph)) {
    is2(mapp_bedGraph, "TabixFile")
    
    mapp = svview_get_mappability(region, mapp_bedGraph)
    mapp = data.frame(score = as.numeric(Rle(values=mapp$score, lengths=mapp$width)),
                      pos   = mapp$pos[1] : (mapp$pos[1] + sum(mapp$width) -1) )
    if(missing(uniqueValue))
        uniqueValue = max(mapp$score)
    data_snps = join(data_snps, mapp, by="pos")
    data_snps$uniqueValue = factor(data_snps$score == uniqueValue, 
                                  levels=c(T,F), 
                                  labels=c("good SNP", "in bad region"))
    plot_snps = ggplot(data_snps) + 
      aes(x=pos, y=alt/(ref+alt), col=sample, shape=sample, alpha = uniqueValue) +
      scale_alpha_discrete(range = c("good SNP" = 1, "in bad region" = 0.2), name = NULL)
  }
  else 
  {
    plot_snps = ggplot(data_snps) + 
      aes(x=pos, y=alt/(ref+alt), col=sample, shape=sample)
  }
  
  plot_snps + 
    geom_point(size=1) + 
    ylab(NULL) + 
    scale_y_continuous(breaks=c(0,0.25,0.33,0.5,0.66,0.75,1))
}
