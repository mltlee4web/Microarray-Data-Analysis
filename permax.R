# By R Gray, DFCI
# Copyright (C) 2000, 2002 Robert Gray
# Modified by Peizheng Chen, 2026

### permax()
###
### Computes row-wise two-group permutation tests for microarray or similar
### high-dimensional data. 
###
### Arguments
###
### data:
###   Numeric matrix or data frame with markers/genes in rows and samples in
###   columns. Marker identifiers used for matching and reporting should be in
###   dimnames(data)[[1]].
###
### ig1:
###   Column numbers belonging to group 1.
###
### nperm:
###   nperm <= 0 computes the complete permutation distribution.
###   nperm > 0 evaluates nperm random permutations.
###
### logs:
###   If TRUE, nonpositive values are first replaced by 1 and the natural
###   logarithm is applied. Summary statistics and non-rank test statistics
###   are then computed from the logged data.
###
### ranks:
###   If TRUE, observations within each row are replaced by ranks before the
###   permutation statistic is calculated. This produces a generalized
###   Wilcoxon rank-sum test. For a stratified analysis, ranking is performed
###   separately within each stratum.
###
### min.np:
###   Retains only rows having at least min.np observations greater than the
###   overall minimum value of the analyzed data. The count is taken across
###   the group 1 and group 2 columns.
###
### ig2:
###   Optional column numbers belonging to group 2. If NULL, all columns not
###   listed in ig1 are assigned to group 2. If supplied, columns outside ig1
###   and ig2 are excluded, and min.np is evaluated only using ig1 and ig2.
###
### WHseed:
###   Initial Wichmann-Hill random-number seed, supplied as a vector of three
###   integers. If omitted for a random-permutation analysis, a seed is
###   generated using runif(). It is not needed when the complete permutation
###   distribution is calculated.
###
### cluster:
###   Optional vector of cluster or stratum identifiers, with one identifier
###   per sample column. When permute.cluster=FALSE, treatment assignments are
###   permuted within clusters while preserving the original number of group 1
###   observations in each cluster. When permute.cluster=TRUE, entire clusters
###   are permuted between groups.
###
### stratify:
###   If TRUE, computes a stratified statistic by combining within-cluster
###   contributions. A cluster vector must be supplied. stratify=TRUE cannot
###   be combined with permute.cluster=TRUE or signed.rank=TRUE.
###
### weights:
###   Optional numeric vector containing one weight per cluster for a
###   stratified analysis. If NULL, equal initial weights of 1/number of
###   clusters are used. The weights are internally adjusted for the group
###   allocation within each informative cluster.
###
### nl, nr:
###   Numbers of markers used to define the lower- and upper-tail critical
###   values. nl selects the lower-tail cutoff from the ordered observed
###   statistics, and nr selects the corresponding upper-tail cutoff. These
###   cutoffs are also used to summarize the permutation distribution of the
###   number of tail-positive markers.
###
### permute.cluster:
###   If TRUE, whole clusters rather than individual samples are permuted.
###   Every cluster must belong entirely to one treatment group. This option
###   cannot be combined with stratify=TRUE.
###
### signed.rank:
###   If TRUE, performs a paired Wilcoxon signed-rank analysis. Every cluster
###   must contain exactly two observations, normally one from each group.
###   This option requires cluster and cannot be combined with stratify=TRUE.
###
### expord:
###   If TRUE, calculates the expected ordered test statistics under the
###   permutation distribution. They are returned in attr(result, "expected").
###
### Value
###
### Returns an object of class c("permax", "data.frame"), ordered from the
### smallest to the largest observed test statistic, with the following
### columns:
###
### stat:
###   Observed test statistic.
###
### pind.lower:
###   Individual, unadjusted lower-tail permutation p-value: the proportion
###   of permutations with a statistic less than or equal to the observed
###   statistic for that row.
###
### pind.upper:
###   Individual, unadjusted upper-tail permutation p-value: the proportion
###   of permutations with a statistic greater than or equal to the observed
###   statistic for that row.
###
### p.lower:
###   Multiplicity-adjusted, step-down lower-tail permutation p-value.
###
### p.upper:
###   Multiplicity-adjusted, step-down upper-tail permutation p-value.
###
### m1, m2:
###   Arithmetic means for groups 1 and 2. If logs=TRUE, these are means of
###   the logged observations.
###
### s1, s2:
###   Sample standard deviations for groups 1 and 2. If logs=TRUE, these are
###   standard deviations of the logged observations.
###
### np1, np2:
###   Numbers of observations in groups 1 and 2 that are greater than the
###   overall minimum value of the analyzed data.
###
### mdiff:
###   Difference between group means. If logs=TRUE, this is the difference
###   between the two geometric means: exp(m1) - exp(m2).
###
### mrat:
###   Ratio of group means. If logs=TRUE, this is the ratio of geometric
###   means: exp(m1 - m2).
###
### Attributes
###
### attr(result, "dist"):
###   Summary of the null permutation distributions for the numbers of
###   lower- and upper-tail positive markers. It contains nl, prop.nl,
###   prop.1l, ave.l, nr, prop.nr, prop.1r, and ave.r.
###
### attr(result, "expected"):
###   Expected order statistics under the permutation distribution. Present
###   only when expord=TRUE.
###
### attr(result, "seed.start"), attr(result, "seed.end"):
###   Starting and ending random-number seeds. Present only when nperm > 0.
###
### attr(result, "call"):
###   The matched call to permax().

permax <- function(data, ig1, nperm=0, logs=TRUE, ranks=FALSE, min.np=1,
                   ig2=NULL, WHseed=NULL, cluster=NULL, stratify=FALSE,
                   weights=NULL, nl=50, nr=50, permute.cluster=FALSE,
                   signed.rank=FALSE, expord=FALSE) {
  cl <- match.call()
  data <- as.matrix(data)
  if (logs) {
    tmp <- data <= 0
    if (any(tmp)) {
      data[tmp] <- 1
    }
    data <- log(data)
  }
  
  if (!is.null(ig2)) {
    # Remove unused columns and adjust ig1
    i2 <- ig1i <- rep(FALSE, ncol(data))
    ig1i[ig1] <- TRUE
    i2[c(ig1, ig2)] <- TRUE
    data <- data[, i2]
    if (!is.null(cluster)) {
      cluster <- cluster[i2]
    }
    ig1 <- seq_len(ncol(data))[ig1i[i2]]
  }
  
  dmin <- min(data)
  n1 <- length(ig1)
  n2 <- ncol(data) - n1
  
  # Descriptive statistics are calculated before filtering and before the
  # permutation-specific transformations.
  d1 <- data[, ig1, drop=FALSE]
  if (n1 > 1L) {
    m1 <- c(d1 %*% rep(1 / n1, n1))
    s1 <- sqrt((d1 - m1)^2 %*% rep(1 / (n1 - 1L), n1))
    d1[d1 <= dmin] <- 0
    d1[d1 > dmin] <- 1
    npos1 <- d1 %*% rep(1, n1)
  } else {
    m1 <- c(d1)
    s1 <- rep(0, length(d1))
    npos1 <- ifelse(d1 > dmin, 1, 0)
  }
  
  d1 <- data[, -ig1, drop=FALSE]
  m2 <- c(d1 %*% rep(1 / n2, n2))
  s2 <- if (n2 > 1L) {
    sqrt((d1 - m2)^2 %*% rep(1 / (n2 - 1L), n2))
  } else {
    rep(0, nrow(d1))
  }
  d1[d1 <= dmin] <- 0
  d1[d1 > dmin] <- 1
  npos2 <- d1 %*% rep(1, n2)
  
  sub <- npos1 + npos2 >= min.np
  data <- data[sub, ]
  n <- nrow(data)
  
  if (!is.null(cluster)) {
    trt <- rep(2, n1 + n2)
    trt[ig1] <- 1
    mclust <- table(cluster)
    nclust <- length(mclust)
    mct1 <- table(cluster, trt)[, 1]
    
    if (permute.cluster) {
      if (stratify) {
        stop("permute.cluster and stratify cannot both == TRUE")
      }
      if (any(mclust != mct1 & mct1 != 0)) {
        stop("clusters cannot contain both groups when permute.cluster==TRUE")
      }
      ipc <- 1L
      o <- order(trt, cluster)
    } else {
      ipc <- 0L
      o <- order(cluster, trt)
    }
    
    data <- data[, o]
    trt <- trt[o]
    ig1 <- seq_len(ncol(data))[trt == 1]
    
    if (stratify) {
      istrt <- 1L
      if (is.null(weights)) {
        weights <- rep(1, nclust) / nclust
      } else if (length(weights) != nclust) {
        stop(paste("weights must have length", format(nclust)))
      }
    } else {
      istrt <- 0L
    }
    
    if (signed.rank) {
      if (max(mclust) != 2L || min(mclust) != 2L) {
        stop("signed.rank requires paired data")
      }
      if (stratify) {
        stop("stratify and signed.rank cannot both = TRUE")
      }
      
      d1 <- data[, ig1] - data[, -ig1]
      d2 <- t(apply(abs(d1), 1, rank))
      d2[d1 == 0] <- 0
      d2 <- ifelse(d1 < 0, -d2, d2)
      data[, ig1] <- d2
      data[, -ig1] <- -d2
      irnk <- 2L
    }
  } else {
    data <- cbind(data[, ig1], data[, -ig1])
    ig1 <- seq_len(n1)
    nclust <- 1L
    mclust <- n1 + n2
    mct1 <- n1
    istrt <- 0L
    ipc <- 0L
  }
  
  if (!signed.rank) {
    if (ranks) {
      if (istrt == 1L) {
        t1 <- c(0, cumsum(mclust))
        for (i in seq_len(nclust)) {
          ii <- (t1[i] + 1L):t1[i + 1L]
          data[, ii] <- t(apply(data[, ii, drop=FALSE], 1, rank))
        }
      } else {
        data <- t(apply(data, 1, rank))
      }
      irnk <- 1L
    } else {
      irnk <- 0L
    }
  }
  
  if (nperm > 0) {
    if (is.null(WHseed)) {
      WHseed <- floor(30000 * runif(3)) + 1
    }
  } else {
    WHseed <- c(0, 0, 0)
    
    if (ipc == 1L) {
      nct1 <- sum(as.numeric(mct1 > 0))
      nn <- exp(sum(log(2:nclust)) - sum(log(2:nct1)) - sum(log(2:(nclust - nct1))))
    } else {
      nn <- 0
      for (i in seq_len(nclust)) {
        if (mclust[i] > mct1[i] && mct1[i] > 0) {
          nn <- nn + sum(log(1:mclust[i])) - sum(log(1:mct1[i])) - sum(log(1:(mclust[i] - mct1[i])))
        }
      }
      nn <- exp(nn)
    }
    
    cat("statistics will be computed for all", format(nn), "combinations\n")
  }
  
  Z <- ptnstd(data=data, ng1=n1, nclust=nclust, mclust=mclust,mct1=mct1, ig1=ig1, 
              irnk=irnk, stratified=istrt, weights=weights, permute_cluster=ipc)
  
  Z2 <- Z$stat
  weights <- Z$weights
  data <- matrix(Z$d, nrow=n, dimnames=dimnames(data))
  
  o <- order(Z2)
  Z2 <- Z2[o]
  
  if (nl > n) {
    nl <- round(n / 2)
  }
  if (nr > n) {
    nr <- round(n / 2)
  }
  
  crit <- c(Z2[nl], Z2[n - nr + 1L])
  data <- data[o, ]
  
  use_t_statistic <- !ranks && !stratify && !signed.rank
  
  permutation <- ptn(data=data, ng1=n1, stat=Z2, nperm=nperm, seed=WHseed,
                     nclust=nclust, mclust=mclust, mct1=mct1,ig1=c(ig1, rep(0, n2)),
                     irnk=irnk, stratified=istrt, weights=weights, nlr=c(nl, nr),
                     crit=crit, permute_cluster=ipc, expord=expord,
                     expected_statistic=if (use_t_statistic) "t" else "sum")
  
  if (nperm > 0) {
    endseed <- permutation$ix
  }
  
  if (use_t_statistic) {
    Z2 <- tst2(cbind(data[, ig1], data[, -ig1]), ng1=n1, ng2=n2)
  }
  
  dist <- c(nl, permutation$dist[1:3],
            nr, permutation$dist[4:6])
  names(dist) <- c("nl", "prop.nl", "prop.1l", "ave.l",
                   "nr", "prop.nr", "prop.1r", "ave.r")
  
  Z <- data.frame(stat=Z2,
                  pind.lower=permutation$pind.lower / permutation$nperm,
                  pind.upper=permutation$pind.upper / permutation$nperm,
                  p.lower=permutation$p.lower / permutation$nperm,
                  p.upper=permutation$p.upper / permutation$nperm)
  
  m1 <- m1[sub]
  m2 <- m2[sub]
  
  if (logs) {
    d1 <- data.frame(m1=m1, m2=m2, s1=s1[sub], s2=s2[sub],
                     np1=npos1[sub], np2=npos2[sub],
                     mdiff=exp(m1) - exp(m2), mrat=exp(m1 - m2))
  } else {
    d1 <- data.frame(m1=m1, m2=m2, s1=s1[sub], s2=s2[sub],
                     np1=npos1[sub], np2=npos2[sub],
                     mdiff=m1 - m2, mrat=m1 / m2)
  }
  
  Z <- cbind(Z, d1[o, ])
  row.names(Z) <- dimnames(data)[[1]]
  class(Z) <- c("permax", "data.frame")
  attr(Z, "dist") <- dist
  attr(Z, "call") <- cl
  
  if (expord) {
    attr(Z, "expected") <- permutation$expected
  }
  
  if (nperm > 0) {
    attr(Z, "seed.start") <- WHseed
    attr(Z, "seed.end") <- endseed
  }
  
  Z
}

ptnstd <- function(data, ng1, nclust, mclust, mct1, ig1,
                   irnk, stratified, weights, permute_cluster) {
  data <- as.matrix(data)
  ng1 <- as.integer(ng1)
  nclust <- as.integer(nclust)
  mclust <- as.integer(mclust)
  mct1 <- as.integer(mct1)
  ig1 <- as.integer(ig1)
  irnk <- as.integer(irnk)
  stratified <- isTRUE(as.logical(stratified))
  invisible(permute_cluster)
  
  if (is.null(weights)) {
    weights_out <- numeric(0)
  } else {
    weights_out <- as.numeric(weights)
  }
  
  standardize_block <- function(columns, center_only) {
    if (length(columns) == 0L) {
      return(invisible(NULL))
    }
    
    block <- data[, columns, drop=FALSE]
    
    for (i in seq_len(nrow(block))) {
      values <- block[i, ]
      value_mean <- sum(values) / length(values)
      
      if (center_only) {
        block[i, ] <- values - value_mean
      } else {
        value_ss <- sum((values - value_mean)^2)
        value_sd <- sqrt(value_ss / (length(values) - 1L))
        
        # stdmv leaves a zero-variance row unchanged.
        if (isTRUE(value_sd > 0)) {
          block[i, ] <- (values - value_mean) / value_sd
        }
      }
    }
    
    data[, columns] <<- block
    invisible(NULL)
  }
  
  if (stratified) {
    if (length(weights_out) != nclust) {
      stop("weights must have one value per cluster for a stratified test")
    }
    
    column_offset <- 0L
    
    for (j in seq_len(nclust)) {
      columns <- column_offset + seq_len(mclust[j])
      
      if (mct1[j] > 0L && mct1[j] < mclust[j]) {
        weights_out[j] <- weights_out[j] * mclust[j] /
          (mct1[j] * (mclust[j] - mct1[j]))
      } else {
        weights_out[j] <- 0
      }
      
      if (irnk == 1L) {
        standardize_block(columns, center_only=TRUE)
      } else if (irnk != 2L) {
        standardize_block(columns, center_only=FALSE)
      }
      
      column_offset <- column_offset + mclust[j]
    }
  } else {
    if (irnk == 1L) {
      standardize_block(seq_len(ncol(data)), center_only=TRUE)
    } else if (irnk != 2L) {
      standardize_block(seq_len(ncol(data)), center_only=FALSE)
    }
  }
  
  stat <- tsum(data=data,
               ig1=ig1[seq_len(ng1)],
               stratified=stratified,
               mclust=mclust,
               mct1=mct1,
               weights=weights_out)
  
  list(d=data, stat=stat, weights=weights_out)
}

# The data matrix is assumed to have already been standardized by ptnstd().
tsum <- function(data, ig1, stratified, mclust, mct1, weights) {
  if (!stratified) {
    return(rowSums(data[, ig1, drop=FALSE]))
  }
  
  ans <- numeric(nrow(data))
  group_offset <- 0L
  
  for (j in seq_along(mclust)) {
    number_group1 <- mct1[j]
    
    if (number_group1 > 0L && number_group1 < mclust[j]) {
      selected <- group_offset + seq_len(number_group1)
      ans <- ans + weights[j] * rowSums(data[, ig1[selected], drop=FALSE])
    }
    
    group_offset <- group_offset + number_group1
  }
  
  ans
}

tst2 <- function(data, ng1, ng2=ncol(data) - ng1) {
  data <- as.matrix(data)
  ng1 <- as.integer(ng1)
  ng2 <- as.integer(ng2)
  ans <- numeric(nrow(data))
  
  group1_columns <- seq_len(ng1)
  group2_columns <- ng1 + seq_len(ng2)
  
  for (j in seq_len(nrow(data))) {
    group1 <- data[j, group1_columns]
    group2 <- data[j, group2_columns]
    
    mean1 <- sum(group1) / ng1
    mean2 <- sum(group2) / ng2
    ss1 <- sum((group1 - mean1)^2)
    ss2 <- sum((group2 - mean2)^2)
    
    if (ss1 == 0 && ss2 == 0) {
      ans[j] <- 0
    } else {
      ans[j] <- (mean1 - mean2) / sqrt((1 / ng1 + 1 / ng2) * (ss1 + ss2) / (ng1 + ng2 - 2L))
    }
  }
  
  ans
}

ptn <- function(data, ng1, stat, nperm, seed, nclust, mclust, mct1, ig1,
                irnk, stratified, weights,nlr, crit, permute_cluster,
                expord=FALSE, expected_statistic=c("sum", "t")) {
  data <- as.matrix(data)
  stat <- as.numeric(stat)
  crit <- as.numeric(crit)
  ng1 <- as.integer(ng1)
  nperm_requested <- as.integer(nperm)
  seed_state <- as.integer(seed)
  nclust <- as.integer(nclust)
  mclust <- as.integer(mclust)
  mct1 <- as.integer(mct1)
  current_group1 <- as.integer(ig1[seq_len(ng1)])
  irnk <- as.integer(irnk)
  stratified <- isTRUE(as.logical(stratified))
  permute_cluster <- isTRUE(as.logical(permute_cluster))
  nlr <- as.integer(nlr)
  expected_statistic <- match.arg(expected_statistic)
  invisible(irnk)
  
  number_genes <- nrow(data)
  number_samples <- ncol(data)
  exhaustive <- nperm_requested <= 0L
  
  lower_individual <- if (exhaustive) {
    rep(1, number_genes)
  } else {
    numeric(number_genes)
  }
  upper_individual <- if (exhaustive) {
    rep(1, number_genes)
  } else {
    numeric(number_genes)
  }
  lower_multiple <- if (exhaustive) {
    rep(1, number_genes)
  } else {
    numeric(number_genes)
  }
  upper_multiple <- if (exhaustive) {
    rep(1, number_genes)
  } else {
    numeric(number_genes)
  }
  
  if (exhaustive) {
    dist_counts <- c(1, 1, nlr[1L], 1, 1, nlr[2L])
    number_evaluated <- 1L
  } else {
    dist_counts <- numeric(6L)
    number_evaluated <- 0L
  }
  
  expected_sum <- if (expord) numeric(number_genes) else NULL
  
  if (exhaustive && expord) {
    if (expected_statistic == "sum") {
      observed_for_expected <- stat
    } else {
      observed_group2 <- setdiff(seq_len(number_samples), current_group1)
      observed_for_expected <- tst2(cbind(data[, current_group1, drop=FALSE],
                                          data[, observed_group2, drop=FALSE]),
                                    ng1=length(current_group1), 
                                    ng2=length(observed_group2))
    }
    expected_sum <- expected_sum + sort(observed_for_expected)
  }
  
  process_grouping <- function(group1) {
    permutation_stat <- tsum(
      data=data,
      ig1=group1,
      stratified=stratified,
      mclust=mclust,
      mct1=mct1,
      weights=weights
    )
    
    lower_individual <<- lower_individual + (permutation_stat <= stat)
    upper_individual <<- upper_individual + (permutation_stat >= stat)
    
    lower_hits <- sum(permutation_stat <= crit[1L])
    upper_hits <- sum(permutation_stat >= crit[2L])
    
    running_maximum <- stat[1L] - 1
    for (i in seq_len(number_genes)) {
      running_maximum <- max(permutation_stat[i], running_maximum)
      if (running_maximum >= stat[i]) {
        upper_multiple[i] <<- upper_multiple[i] + 1
      }
    }
    
    running_minimum <- stat[number_genes] + 1
    for (i in seq.int(number_genes, 1L, by=-1L)) {
      running_minimum <- min(permutation_stat[i], running_minimum)
      if (running_minimum <= stat[i]) {
        lower_multiple[i] <<- lower_multiple[i] + 1
      }
    }
    
    dist_counts[1L] <<- dist_counts[1L] + (lower_hits >= nlr[1L])
    dist_counts[2L] <<- dist_counts[2L] + (lower_hits >= 1L)
    dist_counts[3L] <<- dist_counts[3L] + lower_hits
    dist_counts[4L] <<- dist_counts[4L] + (upper_hits >= nlr[2L])
    dist_counts[5L] <<- dist_counts[5L] + (upper_hits >= 1L)
    dist_counts[6L] <<- dist_counts[6L] + upper_hits
    
    if (expord) {
      if (expected_statistic == "sum") {
        expected_values <- permutation_stat
      } else {
        group2 <- setdiff(seq_len(number_samples), group1)
        expected_values <- tst2(cbind(data[, group1, drop=FALSE],
                                      data[, group2, drop=FALSE]),
                                ng1=length(group1), 
                                ng2=length(group2))
      }
      expected_sum <<- expected_sum + sort(expected_values)
    }
    
    number_evaluated <<- number_evaluated + 1L
    invisible(NULL)
  }
  
  if (exhaustive) {
    if (permute_cluster) {
      selected_clusters <- which(mct1 > 0L)
      
      repeat {
        number_selected <- length(selected_clusters)
        combination_exhausted <- FALSE
        
        if (number_selected == 0L) {
          combination_exhausted <- TRUE
        } else if (selected_clusters[number_selected] < nclust) {
          selected_clusters[number_selected] <-
            selected_clusters[number_selected] + 1L
        } else {
          advanced_position <- 0L
          
          if (number_selected > 1L) {
            for (j in seq.int(number_selected - 1L, 1L, by=-1L)) {
              if (selected_clusters[j] < nclust - number_selected + j) {
                advanced_position <- j
                break
              }
            }
          }
          
          if (advanced_position > 0L) {
            selected_clusters[advanced_position] <-
              selected_clusters[advanced_position] + 1L
            if (advanced_position < number_selected) {
              for (k in seq.int(advanced_position + 1L, number_selected)) {
                selected_clusters[k] <- selected_clusters[k - 1L] + 1L
              }
            }
          } else {
            selected_clusters <- seq_len(number_selected)
            combination_exhausted <- TRUE
          }
        }
        
        cluster_starts <- c(0L, cumsum(mclust))[selected_clusters] + 1L
        current_group1 <- integer(0)
        for (j in seq_along(selected_clusters)) {
          current_group1 <- c(
            current_group1,
            cluster_starts[j] + seq_len(mclust[selected_clusters[j]]) - 1L
          )
        }
        
        if (combination_exhausted) {
          break
        }
        
        process_grouping(current_group1)
      }
    } else {
      repeat {
        data_offset <- 0L
        group_offset <- 0L
        grouping_available <- FALSE
        
        # Advance the product of the within-cluster combination spaces.  The
        # first cluster changes fastest.
        for (j in seq_len(nclust)) {
          number_selected <- mct1[j]
          positions <- if (number_selected > 0L) {
            group_offset + seq_len(number_selected)
          } else {
            integer(0)
          }
          
          cluster_index <- current_group1[positions]
          cluster_maximum <- data_offset + mclust[j]
          cluster_minimum <- data_offset + 1L
          cluster_exhausted <- FALSE
          
          if (number_selected == 0L) {
            cluster_exhausted <- TRUE
          } else if (cluster_index[number_selected] < cluster_maximum) {
            cluster_index[number_selected] <-
              cluster_index[number_selected] + 1L
          } else {
            advanced_position <- 0L
            
            if (number_selected > 1L) {
              for (k in seq.int(number_selected - 1L, 1L, by=-1L)) {
                if (cluster_index[k] <
                    cluster_maximum - number_selected + k) {
                  advanced_position <- k
                  break
                }
              }
            }
            
            if (advanced_position > 0L) {
              cluster_index[advanced_position] <-
                cluster_index[advanced_position] + 1L
              if (advanced_position < number_selected) {
                for (k in seq.int(advanced_position + 1L, number_selected)) {
                  cluster_index[k] <- cluster_index[k - 1L] + 1L
                }
              }
            } else {
              cluster_index <- cluster_minimum + seq_len(number_selected) - 1L
              cluster_exhausted <- TRUE
            }
          }
          
          if (number_selected > 0L) {
            current_group1[positions] <- cluster_index
          }
          
          if (!cluster_exhausted) {
            grouping_available <- TRUE
            break
          }
          
          data_offset <- data_offset + mclust[j]
          group_offset <- group_offset + number_selected
        }
        
        if (!grouping_available) {
          break
        }
        
        process_grouping(current_group1)
      }
    }
  } else {
    if (nperm_requested > 0L) {
      for (sample_number in seq_len(nperm_requested)) {
        if (permute_cluster) {
          number_group1_clusters <- sum(mct1 > 0L)
          available_clusters <- seq_len(nclust)
          selected_clusters <- integer(number_group1_clusters)
          
          if (number_group1_clusters > 0L) {
            for (i in seq_len(number_group1_clusters)) {
              seed_state <- (c(171L, 172L, 170L) * seed_state) %% c(30269L, 30307L, 30323L)
              uniform <- sum(seed_state / c(30269, 30307, 30323)) %% 1
              remaining <- nclust - i + 1L
              position <- as.integer(uniform * remaining + 1)
              selected_clusters[i] <- available_clusters[position]
              available_clusters[position] <- available_clusters[remaining]
            }
          }
          
          selected_clusters <- sort(selected_clusters)
          cluster_starts <- c(0L, cumsum(mclust))[selected_clusters] + 1L
          current_group1 <- integer(0)
          for (j in seq_along(selected_clusters)) {
            current_group1 <- c(
              current_group1,
              cluster_starts[j] +
                seq_len(mclust[selected_clusters[j]]) - 1L
            )
          }
        } else {
          current_group1 <- integer(0)
          data_offset <- 0L
          
          for (j in seq_len(nclust)) {
            number_selected <- mct1[j]
            available_columns <- seq_len(mclust[j])
            selected_columns <- integer(number_selected)
            
            if (number_selected > 0L) {
              for (i in seq_len(number_selected)) {
                seed_state <- (c(171L, 172L, 170L) * seed_state) %% c(30269L, 30307L, 30323L)
                uniform <- sum(seed_state / c(30269, 30307, 30323)) %% 1
                remaining <- mclust[j] - i + 1L
                position <- as.integer(uniform * remaining + 1)
                selected_columns[i] <- available_columns[position]
                available_columns[position] <- available_columns[remaining]
              }
            }
            
            current_group1 <- c(current_group1, selected_columns + data_offset)
            data_offset <- data_offset + mclust[j]
          }
        }
        
        process_grouping(current_group1)
      }
    }
  }
  
  # Enforce the same monotonicity constraints as the final loop.
  if (number_genes >= 2L) {
    for (i in seq.int(2L, number_genes)) {
      lower_multiple[i] <- max(lower_multiple[i - 1L], lower_multiple[i])
      reverse_index <- number_genes - i + 1L
      upper_multiple[reverse_index] <- max(upper_multiple[reverse_index], upper_multiple[reverse_index + 1L])
    }
  }
  
  dist <- dist_counts / number_evaluated
  expected <- if (expord) expected_sum / number_evaluated else NULL
  
  list(pind.lower=lower_individual,
       pind.upper=upper_individual,
       p.lower=lower_multiple,
       p.upper=upper_multiple,
       nperm=number_evaluated,
       ix=seed_state,
       dist=dist,
       expected=expected)
}

plot.expord <- function(x, del=0, ...) {
  # 1. Pull the attribute of expected order
  expected_vals <- attr(x, 'expected')
  
  # 2. Check if expected stats exist
  if (is.null(expected_vals)) {
    stop("This permax object does not contain expected order statistics. \nPlease re-run permax() with the argument 'expord=TRUE'.")
  }
  
  # 3. Sort the Observed Statistics
  obs <- sort(x$stat, decreasing = FALSE)
  exp <- expected_vals 
  
  # 4. Create the Plot
  plot(exp, obs,
       xlab = "Expected Order Statistics",
       ylab = "Observed Statistics",
       main = "Observed vs. Expected Statistics",
       pch = 1,      
       ...)
  
  # 5. Add the Reference Line (Diagonal y = x)
  abline(0, 1, col = "black")
  
  # 6. Add 'del' lines 
  if (del > 0) {
    abline(del, 1, col = "blue", lty = 2)  
    abline(-del, 1, col = "blue", lty = 2) 
  }
}
