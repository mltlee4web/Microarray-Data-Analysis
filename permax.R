### =======================================================================
### FUNCTION: permax
### Performs permutation-based differential expression analysis.
### Handles exact and Monte Carlo permutations, with support for 
### clustered (dependent) and stratified (blocked) experimental designs.
### =======================================================================
### 
### INPUT:
### data      : Data matrix; genes/markers in rows, samples in columns. 
###             Gene codes used for matching should be in dimnames(data)[[1]].
### ig1       : Numeric vector of column indices belonging to Group 1.
### ig2       : (Optional) Numeric vector of column indices for Group 2. 
###             If missing, all columns not in ig1 are assumed to be Group 2.
### nperm     : Number of permutations. 
###             <= 0: compute the exact, full permutation distribution.
###             > 0:  compute nperm random Monte Carlo samples.
### logs      : Logical. If TRUE, summary stats are computed from logs of data, 
###             and logs are used in the t-statistics (if ranks=FALSE).
### ranks     : Logical. If TRUE, ranks are used in the test statistics 
###             (yielding a Wilcoxon/Rank-Sum test).
### min.np    : Integer. Data is subset to only include rows with at least 
###             'min.np' values > min(data) across the columns in ig1 and ig2.
### WHseed    : Initial random number seed (integer vector). If missing, generated 
###             from runif(). Not needed if exact permutations are calculated.
### cluster   : Vector defining cluster memberships for samples. If specified, 
###             samples in the same cluster are permuted together to preserve 
###             within-cluster dependence (e.g., repeated measures on same mouse).
### signed.rank : Logical. Used only if ranks=TRUE. If FALSE (default), performs 
###               Wilcoxon Rank-Sum test. If TRUE, performs Wilcoxon Signed-Rank test.
### stratify  : Vector defining strata/blocks. Permutations are strictly restricted 
###             to occur WITHIN each stratum. Can be combined with 'cluster' 
###             for nested hierarchical designs.
### permute.cluster : Logical. If TRUE and 'cluster' is defined, permutes whole clusters.
### nl        : Integer. Cutoff for the lower-tail False Discovery Rate (FDR). 
###             Uses the nl'th most negative observed statistic as the critical threshold.
### nr        : Integer. Cutoff for the upper-tail False Discovery Rate (FDR).
###             Uses the nr'th most positive observed statistic as the critical threshold.
### expord    : Logical. If TRUE, calculates the Expected Order Statistics across 
###             permutations (used for generating Q-Q plots).
###
### OUTPUT:
### Returns a data.frame (class 'permax') with rows corresponding to surviving genes:
###   stat    : The standardized test statistic (t-stat or rank-sum).
###   pind    : Individual permutation p-values (2-sided).
###   p2      : 2-sided Family-Wise Error Rate (FWER) using dist of max overall rows.
###   p.lower : 1-sided p-value for lower expression in Group 1.
###   p.upper : 1-sided p-value for higher expression in Group 1.
###   nml     : Count of permutations where this row was the most significant for p.lower.
###   nmr     : Count of permutations where this row was the most significant for p.upper.
###   m1, m2  : Means of groups 1 and 2 (geometric means if logs=TRUE).
###   s1, s2  : Standard deviations of groups 1 and 2 (of logs if logs=TRUE).
###   np1, np2: Number of positive signals (> minimum value) in groups 1 and 2.
###   mdiff   : Difference of means (difference of geometric means if logs=TRUE).
###   mrat    : Ratio of means (ratio of geometric means if logs=TRUE).
###
### ATTRIBUTES ATTACHED TO OUTPUT (Access via attr(res, "attribute_name")):
###   expected   : Vector of Expected Order Statistics for QQ-plots (if expord=TRUE).
###   seed.start : The random seed used at the start of the permutations.
###   seed.end   : The state of the random seed at the end of the permutations.
###   fdr.lower  : Estimated Positive False Discovery Rate (pFDR) for the bottom 'nl' genes.
###                (Expected false positives / Actual observed lower hits)
###   fdr.upper  : Estimated Positive False Discovery Rate (pFDR) for the top 'nr' genes.
###                (Expected false positives / Actual observed upper hits)
### =======================================================================

permax <- function(data, ig1, nperm=0, logs=TRUE, ranks=FALSE, min.np=1,
                   ig2, WHseed=NULL, cluster=NULL, signed.rank=FALSE, stratify=NULL, 
                   expord=FALSE, permute.cluster=TRUE, nl=25, nr=25, ...) {
  
  # DATA CLEANING AND REORDERING
  data <- as.matrix(data)
  if (logs) {
    data[data<=0] <- 1
    data <- log(data)
  }
  dmin <- min(data)
  
  if (missing(ig2)) {
    # Splits data into [Group 1 | Everything Else] based on ig1
    data_reordered <- cbind(data[,ig1], data[,-ig1])
    # Reorder cluster/stratify to match the new column order
    if(!is.null(cluster)) cluster <- c(cluster[ig1], cluster[-ig1])
    if(!is.null(stratify)) stratify <- c(stratify[ig1], stratify[-ig1])
    data <- data_reordered
  } else {
    data_reordered <- cbind(data[,ig1], data[,ig2])
    if(!is.null(cluster)) cluster <- c(cluster[ig1], cluster[ig2])
    if(!is.null(stratify)) stratify <- c(stratify[ig1], stratify[ig2])
    data <- data_reordered
  }
  
  # GET TRUE POSITIVE SIGNALS
  n1 <- length(ig1)
  n2 <- ncol(data)-n1
  ig1 <- 1:n1
  
  # Number of positive signals in the first group
  d1 <- data[,ig1,drop=FALSE]
  m1 <- c(d1 %*% rep(1/n1, n1))
  s1 <- if(n1>1) c(sqrt((d1-m1)^2 %*% rep(1/(n1-1), n1))) else rep(0,nrow(d1))
  d1[d1<=dmin] <- 0
  d1[d1>dmin] <- 1
  npos1 <- d1 %*% rep(1,n1)
  
  # Number of positive signals in the second group
  d1 <- data[,-ig1,drop=FALSE]
  m2 <- c(d1 %*% rep(1/n2,n2))
  s2 <- if(n2>1) c(sqrt((d1-m2)^2 %*% rep(1/(n2-1),n2))) else rep(0,nrow(d1))
  d1[d1<=dmin] <- 0
  d1[d1>dmin] <- 1
  npos2 <- d1 %*% rep(1,n2)
  
  # Filter true positive signals
  sub <- npos1+npos2 >= min.np
  data <- cbind(data[sub,ig1],data[sub,-ig1])
  if (ranks) {
    if (signed.rank && !is.null(cluster)) {
      # Wilcoxon signed-rank
      # 1. Identify the paired columns using the cluster IDs
      unique_clusters <- unique(cluster)
      col_ig1 <- numeric(length(unique_clusters))
      col_ig2 <- numeric(length(unique_clusters))
      
      for(i in seq_along(unique_clusters)) {
        c_id <- unique_clusters[i]
        cols <- which(cluster == c_id)
        # Because data was reordered above, Group 1 is columns 1 to n1
        col_ig1[i] <- cols[cols <= n1]
        col_ig2[i] <- cols[cols > n1]
      }
      
      # 2. Calculate the differences (Group 1 - Group 2)
      diffs <- data[, col_ig1, drop=FALSE] - data[, col_ig2, drop=FALSE]
      
      # 3. Calculate Signed Ranks across the differences
      signed_ranks <- t(apply(diffs, 1, function(x) sign(x) * rank(abs(x))))
      
      # 4. Create the Mirror Matrix for the permutation engine!
      data[, col_ig1] <- signed_ranks
      data[, col_ig2] <- -signed_ranks
      
    }
    else {
      # Wilcoxon rank-sum
      data <- t(apply(data,1,rank))
    }
  }
  
  # Calculate number of permutation (by observation/cluster) and specify the WHseed
  if (nperm <= 0) {
    if (!is.null(stratify)) {
      # 1. STRATIFIED MATH
      nn <- 1
      for (s in unique(stratify)) {
        idx <- which(stratify == s)
        if (!is.null(cluster) && permute.cluster) {
          mk <- length(unique(cluster[idx]))
          n1k <- length(unique(cluster[intersect(ig1, idx)]))
        } else {
          mk <- length(idx)
          n1k <- length(intersect(ig1, idx))
        }
        nn <- nn * exp(lchoose(mk, n1k))
      }
    } else if (!is.null(cluster) && !permute.cluster) {
      # 2. PERMUTE WITHIN CLUSTERS (Paired Design)
      nn <- 1
      for (c_id in unique(cluster)) {
        idx <- which(cluster == c_id)
        n1k <- length(intersect(ig1, idx))
        nn <- nn * exp(lchoose(length(idx), n1k))
      }
    } else if (!is.null(cluster) && permute.cluster) {
      # 3. CLUSTERED MATH (Permute whole clusters)
      n_clust_total <- length(unique(cluster))
      n_clust_g1 <- length(unique(cluster[ig1]))
      nn <- exp(lchoose(n_clust_total, n_clust_g1)) 
    } else {
      # 4. STANDARD MATH
      nn <- exp(lchoose(n1+n2, n1))
    }
    
    cat('statistics will be computed for all', format(round(nn)), 'groupings\n')
    nperm <- round(nn)
    WHseed <- c(0,0,0)
    
  } else if (is.null(WHseed)) {
    WHseed <- floor(30000*runif(3)) + 1
  }
  
  # PERFORM PERMUTATION TEST
  Z_list <- perm_test(data = data, ig1 = ig1, nperm = nperm, seed = WHseed, 
                      cluster = cluster, stratify = stratify, expord = expord, 
                      permute_cluster = permute.cluster)
  
  # COMPUTE ACTUAL FDR RATES
  crit_lower <- sort(Z_list$stat)[nl]
  crit_upper <- sort(Z_list$stat, decreasing = TRUE)[nr]
  obs_hits_lower <- sum(Z_list$stat <= crit_lower)
  obs_hits_upper <- sum(Z_list$stat >= crit_upper)
  fdr_lower <- if (obs_hits_lower > 0) mean(Z_list$fdr.lower.counts) / obs_hits_lower else NA
  fdr_upper <- if (obs_hits_upper > 0) mean(Z_list$fdr.upper.counts) / obs_hits_upper else NA
  
  if (nperm>0) endseed <- Z_list$ix
  Z <- data.frame(stat=Z_list$stat,pind=Z_list$pind/nperm,p2=Z_list$p2/nperm,
                  p.lower=Z_list$p.lower/nperm,p.upper=Z_list$p.upper/nperm,
                  nml=Z_list$nml,nmr=Z_list$nmr)
  
  # PROCESS OUTPUT TABLE
  # Append descriptive statistics
  if (logs){
    Z <- cbind(Z,m1=m1[sub],m2=m2[sub],s1=s1[sub],s2=s2[sub],np1=npos1[sub],
               np2=npos2[sub],mdiff=exp(m1[sub])-exp(m2[sub]),mrat=exp(m1[sub]-m2[sub]))
  } else {
    Z <- cbind(Z,m1=m1[sub],m2=m2[sub],s1=s1[sub],s2=s2[sub],np1=npos1[sub],
               np2=npos2[sub],mdiff=m1[sub]-m2[sub],mrat=m1[sub]/m2[sub])
  }
  
  row.names(Z) <- dimnames(data)[[1]]
  class(Z) <- c('permax','data.frame')
  
  if (expord) attr(Z, 'expected') <- Z_list$expord
  if (nperm > 0) {
    attr(Z, 'seed.start') <- WHseed
    attr(Z, 'seed.end') <- endseed
  }
  attr(Z, 'fdr.lower') <- fdr_lower
  attr(Z, 'fdr.upper') <- fdr_upper
  
  return(Z)
}

perm_test <- function(data, ig1, nperm, seed, 
                      cluster=NULL, stratify=NULL, expord=FALSE,
                      permute_cluster=TRUE, nl=25, nr=25) {
  
  # SETUP
  n_genes <- nrow(data)
  n_samples <- ncol(data)
  n1 <- length(ig1)
  
  if (!is.null(seed)) set.seed(seed[1])
  
  # STANDARDIZE DATA
  data_std <- matrix(0, nrow = n_genes, ncol = n_samples)
  data_std_weighted <- matrix(0, nrow = n_genes, ncol = n_samples) # Optimization matrix
  obs_sum <- numeric(n_genes)
  
  if (!is.null(stratify)) {
    # Standardize WITHIN each stratum and calculate weighted statistic
    unique_strata <- unique(stratify)
    K <- length(unique_strata) # Total number of strata
    
    for (s in unique_strata) {
      idx <- which(stratify == s)
      stratum_data <- data[, idx, drop=FALSE]
      
      # Standardize within stratum
      stratum_means <- rowMeans(stratum_data)
      stratum_sds <- apply(stratum_data, 1, sd)
      stratum_sds[stratum_sds == 0] <- 1
      data_std[, idx] <- (stratum_data - stratum_means) / stratum_sds
      
      # Stratified Math: (1/K) * [m_k / (n1k * n2k)] * sum(Group 1 in stratum)
      ig1_in_stratum <- intersect(ig1, idx)
      n1k <- length(ig1_in_stratum)
      mk <- length(idx)
      n2k <- mk - n1k
      
      if (n1k > 0 && n2k > 0) {
        stratum_multiplier <- (1/K) * (mk / (n1k * n2k))
        sum_group1 <- rowSums(data_std[, ig1_in_stratum, drop=FALSE])
        obs_sum <- obs_sum + (stratum_multiplier * sum_group1)
        
        data_std_weighted[, idx] <- data_std[, idx] * stratum_multiplier # Optimization matrix
      }
    }
  } else {
    # Standardize GLOBALLY
    row_means <- rowMeans(data)
    row_sds <- apply(data, 1, sd)
    row_sds[row_sds == 0] <- 1 
    data_std <- (data - row_means) / row_sds
    
    obs_sum <- rowSums(data_std[, ig1, drop=FALSE])
    data_std_weighted <- data_std
  }
  
  # FDR cutoff
  crit_lower <- sort(obs_sum)[nl]
  crit_upper <- sort(obs_sum, decreasing = TRUE)[nr]
  perm_lower_counts <- numeric(nperm)
  perm_upper_counts <- numeric(nperm)
  
  # STRATIFY and CLUSTER SETUP
  if (!is.null(stratify)) {
    exact_nperm <- 1
    for (s in unique_strata) {
      idx <- which(stratify == s)
      
      if (!is.null(cluster) && permute_cluster) {
        # STRATIFIED and CLUSTERED: Count unique CLUSTERS within this stratum
        mk <- length(unique(cluster[idx]))
        n1k <- length(unique(cluster[intersect(ig1, idx)]))
      } else {
        # STRATIFIED only: Count raw columns
        mk <- length(idx)
        n1k <- length(intersect(ig1, idx))
      }
      exact_nperm <- exact_nperm * exp(lchoose(mk, n1k))
    }
    
  } else if (!is.null(cluster) && !permute_cluster) {
    # PAIRED design: 2^N
    exact_nperm <- 1 
    for (c_id in unique(cluster)) {
      idx <- which(cluster == c_id)
      exact_nperm <- exact_nperm * exp(lchoose(length(idx), length(intersect(ig1, idx))))
    }
  } else if (!is.null(cluster) && permute_cluster) {
    # CLUSTERED ONLY
    exact_nperm <- exp(lchoose(length(unique(cluster)), length(unique(cluster[ig1]))))
  } else {
    # STANDARD
    exact_nperm <- exp(lchoose(n_samples, n1))
  }
  
  is_exact <- (round(nperm) == round(exact_nperm))

  # Initialize counters
  count_pind <- numeric(n_genes)
  count_p2 <- numeric(n_genes)
  count_lower <- numeric(n_genes)
  count_upper <- numeric(n_genes)
  count_nml <- numeric(n_genes)
  count_nmr <- numeric(n_genes)
  max_stats_dist <- numeric(nperm)
  if (expord) sum_ordered_stats <- numeric(n_genes) else sum_ordered_stats <- NULL
  
  # GENERATE PERMUTATION INDICES
  if (is_exact) {
    if (!is.null(stratify)) {
      # 1. STRATIFIED EXACT (Nested or Pure)
      stratum_combos_list <- list()
      
      for (s in unique(stratify)) {
        idx <- which(stratify == s)
        
        if (!is.null(cluster) && permute_cluster) {
          # STRATIFIED and CLUSTERED: Clusters within Stratum
          stratum_clusters <- unique(cluster[idx])
          g1_clusters <- unique(cluster[intersect(ig1, idx)])
          n1k <- length(g1_clusters)
          stratum_combos_list[[as.character(s)]] <- combn(stratum_clusters, n1k, simplify = FALSE)
        } else {
          # STRATIFIED only: Raw columns within Stratum
          n1k <- length(intersect(ig1, idx))
          stratum_combos_list[[as.character(s)]] <- combn(idx, n1k, simplify = FALSE)
        }
      }
      # Cross all strata to get every valid full-sample permutation
      grid_indices <- expand.grid(lapply(stratum_combos_list, seq_along))
      all_combos <- apply(grid_indices, 1, function(row_idx) {
        unlist(mapply(function(lst, i) lst[[i]], stratum_combos_list, row_idx, SIMPLIFY = FALSE))
      })
      
    } else if (!is.null(cluster) && !permute_cluster) {
      cluster_combos_list <- list()
      for (c_id in unique(cluster)) {
        idx <- which(cluster == c_id)
        n1k <- length(intersect(ig1, idx))
        cluster_combos_list[[as.character(c_id)]] <- combn(idx, n1k, simplify = FALSE)
      }
      grid_indices <- expand.grid(lapply(cluster_combos_list, seq_along))
      all_combos <- apply(grid_indices, 1, function(row_idx) {
        unlist(mapply(function(lst, i) lst[[i]], cluster_combos_list, row_idx, SIMPLIFY = FALSE))
      })
    } else if (!is.null(cluster) && permute_cluster) {
      # 2. CLUSTERED EXACT (No Strata)
      all_combos <- combn(unique(cluster), length(unique(cluster[ig1])))
      
    } else {
      # 3. STANDARD EXACT
      all_combos <- combn(1:n_samples, n1)
    }
  }
  
  # PERMUTATION LOOP
  for (i in 1:nperm) {
    
    # STEP 1: GET PERMUTATION INDICES
    if (is_exact) {
      # Deterministic: USE the pre-calculated combination
      curr_perm_idx <- all_combos[, i]
      
      # If clustered, convert cluster IDs back to actual column indices
      if (!is.null(cluster) && permute_cluster) {
        curr_perm_idx <- which(cluster %in% curr_perm_idx)
      }
    } else {
      # Random: Build a restricted shuffle on the fly
      curr_perm_idx <- integer(0)
      
      if (!is.null(stratify)) {
        # 1. STRATIFIED RANDOM
        for (s in unique_strata) {
          idx <- which(stratify == s)
          
          if (!is.null(cluster) && permute_cluster) {
            # STRATIFIED and CLUSTERED: Sample clusters within this stratum
            stratum_clusters <- unique(cluster[idx])
            n1k <- length(unique(cluster[intersect(ig1, idx)]))
            
            if (length(stratum_clusters) == 1) {
              sampled_clusters <- stratum_clusters
            } else {
              sampled_clusters <- sample(stratum_clusters, n1k)
            }
            curr_perm_idx <- c(curr_perm_idx, idx[cluster[idx] %in% sampled_clusters])
            
          } else {
            # STRATIFIED only: Sample columns within this stratum
            n1k <- length(intersect(ig1, idx))
            
            if (length(idx) == 1) {
              sampled_cols <- idx
            } else {
              sampled_cols <- sample(idx, n1k)
            }
            curr_perm_idx <- c(curr_perm_idx, sampled_cols)
          }
        }
      } else if (!is.null(cluster) && !permute_cluster) {
          # Random sampler for Paired Design
          for (c_id in unique(cluster)) {
            idx <- which(cluster == c_id)
            n1k <- length(intersect(ig1, idx))
            sampled_cols <- if(length(idx)==1) idx else sample(idx, n1k)
            curr_perm_idx <- c(curr_perm_idx, sampled_cols)
          }
        } else if (!is.null(cluster) && permute_cluster) {
        # 2. CLUSTERED RANDOM (No Strata)
        n_c1 <- length(unique(cluster[ig1]))
        unique_clusters <- unique(cluster)
        
        if (length(unique_clusters) == 1) {
          sampled_clusters <- unique_clusters
        } else {
          sampled_clusters <- sample(unique_clusters, n_c1)
        }
        curr_perm_idx <- which(cluster %in% sampled_clusters)
        
      } else {
        # 3. STANDARD RANDOM
        if (n_samples == 1) {
          curr_perm_idx <- 1
        } else {
          curr_perm_idx <- sample(1:n_samples, n1)
        }
      }
    }
    
    # STEP 2: CALCULATE FAKE STATISTIC
    perm_sum <- rowSums(data_std_weighted[, curr_perm_idx, drop=FALSE])
    perm_lower_counts[i] <- sum(perm_sum <= crit_lower)
    perm_upper_counts[i] <- sum(perm_sum >= crit_upper)
    
    # === STEP 3: TALLY THE SCOREBOARDS ===
    abs_obs <- abs(obs_sum)
    abs_perm <- abs(perm_sum)
    
    # Individual P-value tallies (two-sided)
    count_pind <- count_pind + (abs_perm >= abs_obs)
    
    # Track the maximum statistic for the Westfall-Young FWER correction
    max_perm_stat <- max(abs_perm)
    max_stats_dist[i] <- max_perm_stat
    
    # One-sided tallies
    count_lower <- count_lower + (perm_sum <= obs_sum)
    count_upper <- count_upper + (perm_sum >= obs_sum)
    
    # Track which genes achieved the absolute maximum score across the whole matrix
    is_max <- (abs_perm == max_perm_stat)
    count_nml <- count_nml + (is_max & (perm_sum < 0)) 
    count_nmr <- count_nmr + (is_max & (perm_sum > 0))
    
    # Build Expected Order QQ-Plot Data
    if (expord) sum_ordered_stats <- sum_ordered_stats + sort(perm_sum, decreasing = FALSE)
  }
  
  # CALCULATE FAMILY-WISE ERROR RATE
  for (g in 1:n_genes) {
    count_p2[g] <- sum(max_stats_dist >= abs(obs_sum[g]))
  }
  if (expord) {final_expord <- sum_ordered_stats / nperm}
    else {final_expord <- NULL}
  
  return(list(
    stat = obs_sum, pind = count_pind, p2 = count_p2,
    p.lower = count_lower, p.upper = count_upper,
    nml = count_nml, nmr = count_nmr, ix = seed, expord = final_expord,
    fdr.lower.counts = perm_lower_counts, fdr.upper.counts = perm_upper_counts
  ))
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
