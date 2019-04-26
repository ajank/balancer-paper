args <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
dataset <- args[2]
bin_size <- as.integer(args[3])
stopifnot(!is.na(bin_size))

require(data.table)
require(ggplot2)
require(locfit)
require(MASS)

source("src/R/functions_HiCExplorer.R")

id <- paste0(genome, "_", dataset, "_filtered_", bin_size)
fout <- paste0("data/distance_decay/decay_", id, ".rda")

# http://www.win-vector.com/blog/2014/05/trimming-the-fat-from-glm-models-in-r/
stripGlmLR <- function(cm) {
  cm$y <- c()
  cm$model <- c()

  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()  
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()

  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()
  attr(cm$terms, ".Environment") <- c()
  attr(cm$formula, ".Environment") <- c()

  cm
}

stripLmLR <- function(cm) {
  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$model <- c()

  attr(cm$terms, ".Environment") <- c()

  cm
}

stripLocfit <- function(cm) {
  cm$locfit$trans <- NULL
  cm$locfit$frame <- NULL

  attr(cm$terms, ".Environment") <- c()
  environment(cm$trans) <- emptyenv()
  cm$frame <- c()

  cm
}

do_fit <- function(dataset, genome, bin_size, raw)
{
  sample_and_fit <- function(dt)
  {
    fit <- list()
    dt_sel <- dt[distance > 0, ]
    fit$sample <- dt_sel[sample(nrow(dt_sel), min(nrow(dt_sel), 50000L)), ]

    # deals properly with zeros, but not particularly useful
    print(system.time(fit$locfit_nolog <- stripLocfit(locfit(value ~ lp(log(distance), nn = 0.5), data = dt_sel))))

    if (raw)
    {
      # the most reasonable is the Poisson fit, non-negative-defined
      print(system.time(try(fit$fit.glm.poisson <- stripGlmLR(glm(value ~ log(distance), family = poisson(link = "log"), data = dt_sel)))))
    }
    else
    {
      # these two fits do not deal properly with zeros
      print(system.time(fit$fit.lm <- stripLmLR(lm(log(value) ~ log(distance), data = dt_sel[value > 0, ]))))
      # print(system.time(fit$locfit <- stripLocfit(locfit(log(value) ~ lp(log(distance), nn = 0.5), data = dt_sel[value > 0, ]))))

      # Gaussian is not really appropriate
      print(system.time(try(fit$fit.glm.gaussian <- stripGlmLR(glm(value ~ log(distance), family = gaussian(link = "log"), data = dt_sel, start = fit$fit.lm$coefficients)))))
    }
  
    # negative binomial fit takes ~10x longer than Poisson or Gaussian
    # print(system.time(try(fit$fit.glm.nb <- stripGlmLR(glm.nb(value ~ log(distance), data = dt_sel)))))
    return(fit)
  }

  stopifnot(!is.na(bin_size))
  fin <- paste0("data/hicexplorer/txt/", id, ifelse(raw, "", "_corrected"), ".txt.gz")

  fm <- fread_map(fin, genome = genome, bin_size = bin_size)
  t <- fm$map
  chrom_map <- fm$chrom_map

  # if (raw)
  # {
    # account for the fact that values on the diagonal are NOT double-counted
    # stopifnot(diag(t) %% 2L == 0L)
    # diag(t) <- diag(t) %/% 2L
  # }

  bins_removed <- apply(t, 2, sum) == 0 & !raw
  t[bins_removed, ] <- NA
  t[, bins_removed] <- NA
  print(summary(as.vector(t)))

  dt <- NULL
  fit <- list()
  for (chrom in grep("Het$", chrom_map$chrom, value = T, invert = T))
  {
    print(chrom)
    m <- match(chrom, chrom_map$chrom)
    int <- chrom_map$map_start[m]:chrom_map$map_end[m]

    dt_chrom <- as.data.table(melt(t[int, int]))
    dt_chrom$chrom <- factor(chrom, chrom_map$chrom)
    dt_chrom$distance <- abs(dt_chrom$Var2 - dt_chrom$Var1) * bin_size
    dt_chrom$Var1 <- NULL
    dt_chrom$Var2 <- NULL
    dt_chrom <- dt_chrom[!is.na(value), ]

    fit[[chrom]] <- sample_and_fit(dt_chrom)
    fit[[chrom]]$mean_cis <- mean(dt_chrom$value)
    fit[[chrom]]$mean_trans <- NA
    fit[[chrom]]$bins <- length(int)
    fit[[chrom]]$bins_after_thresholding <- sum(!bins_removed[int])
    fit[[chrom]]$count <- sum(t[int, int], na.rm = T)

    dt <- rbind(dt, dt_chrom)
  }

  process_fit <- function(fit, id, pattern)
  {
    print(id)
    dt_id <- dt[grepl(pattern, chrom), ]

    d <- nrow(t)
    bins_id <- rep(F, d)
    for (i in seq_len(nrow(chrom_map)))
      if (grepl(pattern, chrom_map$chrom[i]))
        bins_id[chrom_map$map_start[i]:chrom_map$map_end[i]] <- T

    bins_cis <- matrix(F, d, d)
    bins_trans <- matrix(F, d, d)
    bins_trans[bins_id, ] <- T
    bins_trans[, bins_id] <- T
    for (i in seq_len(nrow(chrom_map)))
      if (grepl(pattern, chrom_map$chrom[i]))
      {
        bins_cis[chrom_map$map_start[i]:chrom_map$map_end[i], chrom_map$map_start[i]:chrom_map$map_end[i]] <- T
        bins_trans[chrom_map$map_start[i]:chrom_map$map_end[i], chrom_map$map_start[i]:chrom_map$map_end[i]] <- F
      }

    fit[[id]] <- sample_and_fit(dt_id)
    fit[[id]]$mean_cis <- mean(dt_id$value)
    fit[[id]]$mean_trans <- mean(t[bins_trans], na.rm = T)
    fit[[id]]$bins <- sum(bins_id)
    fit[[id]]$bins_after_thresholding <- sum(!bins_removed[bins_id])
    fit[[id]]$count <- sum(t[bins_id, bins_id], na.rm = T)

    # print(nrow(dt))
    # print(nrow(dt_id))
    # print(fit[[id]]$mean_cis)
    # print(fit[[id]]$mean_trans)
    # print(fit[[id]]$bins)
    # print(fit[[id]]$bins_after_thresholding)
    # print(sum(bins_cis))
    # print(sum(bins_trans))
    # print(mean(t[bins_cis], na.rm = T))
    stopifnot(abs(fit[[id]]$mean_cis - mean(t[bins_cis], na.rm = T)) < 1e-10)

    return(fit)
  }

  fit <- process_fit(fit, "genome", "^chr[23][LR]$|^chr[4XY]$")
  fit <- process_fit(fit, "chr23", "^chr[23][LR]$")
  fit <- process_fit(fit, "chr23X", "^chr[23][LR]$|^chr[X]$")
  fit <- process_fit(fit, "chr234X", "^chr[23][LR]$|^chr[4X]$")

  return(fit)
}

fit_raw <- do_fit(dataset, genome, bin_size, raw = T)
save(fit_raw, file = fout)

# fit_norm <- do_fit(dataset, genome, bin_size, raw = F)
# save(fit_raw, fit_norm, file = fout)
