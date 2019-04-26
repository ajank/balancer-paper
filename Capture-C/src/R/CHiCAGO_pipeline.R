options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 6)
  stop("Usage:  Rscript CHiCAGO_pipeline.R dataset RepN1...RepNk alleleSuffix designDir inputDir outputDir")

dataset <- options[1]
replString <- options[2]
alleleSuffix <- options[3]
designDir <- options[4]
inputDir <- options[5]
outputDir <- options[6]

replVector <- paste0("Rep", Filter(function(s) s != "", strsplit(replString, "Rep")[[1]]))
message("Using replicates: ", paste(replVector, collapse = " "))

require(Chicago)
options(warn = 1)

#
#  modified chicagoPipeline
#

unzipFifo <- function(inputFile)
{
  fifo <- system("mktemp", intern = T)
  system(paste("zcat", inputFile, ">", fifo))
  return(fifo)
}

myChicago <- function(repl, distFunParams = NULL, providedDistCurve = FALSE, plotAndExport = TRUE)
{
  message('\n\n*** myChicago(dataset = "', dataset, '", alleleSuffix = "', alleleSuffix, '", designDir = "', designDir, '", inputDir = "', inputDir, '", outputDir = "', outputDir, '", repl = "', repl, '", distFunParams = "', distFunParams, '", providedDistCurve = "', providedDistCurve, '", plotAndExport = "', plotAndExport, '")\n')

  inputFiles <- paste0(inputDir, "/", dataset, "_", repl, alleleSuffix, ".chinput.gz")
  inputFilesUncompressed <- sapply(inputFiles, unzipFifo)
  outPrefix <- paste0(dataset, "_", paste0(repl, collapse = ""), alleleSuffix)

  settingsFile <- paste0(designDir, "/dm6_DpnII.settingsFile")
  cd <- setExperiment(designDir = designDir, settingsFile = settingsFile)

  if (providedDistCurve)
  {
    message("\n*** Using modified settings (estimated using fitDistCurve.R)...\n")
    modifySettingsFile <- paste0(outputDir, "/", outPrefix, ".settings")
    system(paste("cat", modifySettingsFile))
    cd <- modifySettings(cd, settingsFile = modifySettingsFile)
  }

  if (length(inputFiles) > 1)
    cd <- readAndMerge(files = inputFilesUncompressed, cd = cd)
  else
    cd <- readSample(file = inputFilesUncompressed, cd = cd)
  unlink(inputFilesUncompressed)

  # modified chicagoPipeline below

  message("\n*** Running normaliseBaits...\n")
  cd <- normaliseBaits(cd)

  message("\n*** Running normaliseOtherEnds...\n")
  suppressWarnings(dir.create(paste0(outputDir, "/diag_plots")))
  cd <- normaliseOtherEnds(cd, outfile = paste0(outputDir, "/diag_plots/", outPrefix, "_oeNorm.pdf"))

  message("\n*** Running estimateTechnicalNoise...\n")
  cd <- estimateTechnicalNoise(cd, outfile = paste0(outputDir, "/diag_plots/", outPrefix, "_techNoise.pdf"))

  if (is.null(distFunParams))
  {
    message("\n*** Running estimateDistFun...\n")
    ### Note that f is saved in cd@params
    cd <- estimateDistFun(cd, outfile = paste0(outputDir, "/diag_plots/", outPrefix, "_distFun.pdf"))
  }
  else
  {
    message("\n*** Using provided distFunParams...\n")
    cd@params$distFunParams <- distFunParams
  }

  ### Note that f is saved as cd@params$f and
  ### subset is saved as cd@settings$brownianNoise.subset
  message("\n*** Running estimateBrownianComponent...\n")
  cd <- estimateBrownianComponent(cd)

  message("\n*** Running getPvals...\n")
  cd <- getPvals(cd)

  message("\n*** Running getScores...\n")
  cd <- getScores(cd)

  message("\n\nSaving the Chicago object...\n")
  suppressWarnings(dir.create(outputDir))
  saveRDS(cd, paste0(outputDir, "/", outPrefix, ".Rds"))

  logfile <- file(paste0(outputDir, "/", outPrefix, "_params.txt"))
  sink(logfile)
  cat("#  CHiCAGO pipeline settings (chicagoData@settings):\n")
  for (s in names(cd@settings))
    cat(paste(s, cd@settings[[s]], sep="\t"), "\n")
  sink(NULL)
  close(logfile)

  if (plotAndExport)
  {
    message("\n\nPlotting examples...\n")
    suppressWarnings(dir.create(paste0(outputDir, "/examples")))
    baits <- plotBaits(cd, outfile = paste0(outputDir, "/examples/", outPrefix, "_proxExamples.pdf"))

    message("\n\nExporting peak lists...\n")
    exportResults(cd, outfileprefix = paste0(outputDir, "/", outPrefix))
  }

  return(cd)
}

#
#  main loop
#

message("\n\n*** To estimate the weights for p-value weighting procedure, we have to consider the replicates separately...\n")
for (repl in replVector)
  myChicago(repl = repl, plotAndExport = F)

message("\n\n*** Do the estimation, using external script fitDistCurve.R...\n")
# # save the four tuning parameters to a .settings file
# # https://www.bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html#using-different-weights
inputs <- paste(paste0(outputDir, "/", dataset, "_", replVector, alleleSuffix, ".Rds"), collapse = ",")
suppressWarnings(dir.create(outputDir))
system(paste0("Rscript /g/furlong/jankowsk/chicago/chicagoTools/fitDistCurve.R --largeBinSize 200000 --threshold ' -5' --inputs ", inputs,
  " ", outputDir, "/", dataset, "_", paste0(replVector, collapse = ""), alleleSuffix))

message("\n\n*** Process the full set of captured regions...\n")
myChicago(repl = replVector, providedDistCurve = T)
