#' event_burstiness
#'
#' @param <ts> character (or logicial) vector
#' @param <min_iet> minimum interevent sequence spacing
#' @keywords burstiness
#' @description This function returns the burstiness coefficient for each of the
#' unique codes in a series. The series must be event-based; see
#' time_burstiness() to find the burstiness of codes in a time-series.
#' @export
#' @references Kim, E.-K., & Jo, H.-H. (2016). Measuring burstiness for finite
#' event sequences. Physical Review E, 94(3), 032311.
#' https://doi.org/10.1103/PhysRevE.94.032311
#' @examples
#' event_burstiness(guastello, 1)
#' event_burstiness(guastello)
#'
event_burstiness <- function(ts,
                             min_iet = NULL) {
  # Takes a sequence of logical or character entries, and returns
  #  the Kim & Jo (2016) burstiness for each of the unique entries

  if(is.character(ts) == FALSE &         # Validate sequence
     is.logical(ts) == FALSE) {
    cat("Invalid sequence type passed to event_burstiness().\n")
    cat("Did you mean to use time_burstiness()?\n")
    return(NULL)
  }

  if(length(ts) > 0) {                   # If ts has entries
    sort(unique(ts)) ->                  # Unique codes
      code_vec
  } else {
    cat("Empty sequence.\n")
    return(NULL)
  }

  if(length(code_vec) > 0) {             # If there are valid codes
    rep(0.0, length(code_vec)) ->        # Pre-allocation
      B
  } else {
    cat("No valid codes found.\n")
    return(NULL)
  }

  if(is.null(min_iet) == TRUE) {         # Use Kim & Jo, eqn. 22
    # No minimum inter-event time.
    for(index in 1:length(code_vec)) {   # For each code
      which(ts == code_vec[index]) ->    # Get the positions of the current
        code_pos_vec                     #  code.
      diff(code_pos_vec, 1) ->           # Interevent times
        gaps_vec
      mean(gaps_vec, na.rm = TRUE) ->    # mean(interevent times)
        mean_iet
      sd(gaps_vec, na.rm = TRUE) ->      # sd(interevent times)
        sd_iet
      sd_iet / mean_iet ->               # Coefficient of variation
        r
      length(code_pos_vec) ->            # Number of interevent times
        n
      (r * sqrt(n+1) - sqrt(n-1)) /      # Burstiness (K & J, eq. 22)
        (r * (sqrt(n+1) - 2) + sqrt(n-1)) ->
        B[index]
    }
    code_vec -> names(B)                 # Identify each B
    return(B)
  }

  if(is.null(min_iet) == FALSE) {      # Minimum inter-event time
    if(min_iet %% 1 == 0) {            # Must be an integer!
      if(min_iet >= 1) {               # Must be a positive integer!
        min_iet / length(ts) ->        # For K & J, eq. 28
          y_tilde
        for(index in 1:length(code_vec)) {   # For each code
          which(ts == code_vec[index]) ->    # Get the positions of the current
            code_pos_vec                     #  code.
          diff(code_pos_vec, 1) ->           # Interevent times
            gaps_vec
          mean(gaps_vec, na.rm = TRUE) ->    # mean(interevent times)
            mean_iet
          sd(gaps_vec, na.rm = TRUE) ->      # sd(interevent times)
            sd_iet
          sd_iet / mean_iet ->               # Coefficient of variation
            r
          length(code_pos_vec) ->            # Number of interevent times
            n
          ((n - 2) * (r * sqrt(n+1) - (1 - n * y_tilde) * sqrt(n-1))) /
            (r * (n * sqrt(n+1) - 2*(n-1)) +
               (1 - n * y_tilde) * sqrt(n-1) * (n - 2 * sqrt(n+1))) ->
            B[index]                         # Burstiness (K & J, eq. 28)
        }
        code_vec -> names(B)
        return(B)
      }
      cat("Invalid minimum inter-event time encountered!/n") # Non-positive
      cat("Minimum inter-event time must be NULL or a positive integer./n")
      cat(paste(min_iet, "is invalid./n"))
      return(NULL)
    }
    cat("Invalid minimum inter-event time encountered!/n") # Non-integer
    cat("Minimum inter-event time must be NULL or a positive integer./n")
    cat(paste(min_iet, "is invalid./n"))
    return(NULL)
  }

  cat("Invalid data/n")            # Should never get here
  cat("Did you mean to use time_burstiness()?/n")
  return(NULL)
}


#' time_burstiness
#' @param <ts> numeric vector
#' @param <min_iet> minimum interevent sequence spacing
#' @keywords burstiness
#' @description This function returns the burstiness coefficient for a series
#' of events. The series must be time-based; see event_burstiness() to find
#' the burstiness of codes in a event-series.
#' @export
#' @references Kim, E.-K., & Jo, H.-H. (2016). Measuring burstiness for finite
#' event sequences. Physical Review E, 94(3), 032311.
#' https://doi.org/10.1103/PhysRevE.94.032311
#' @examples
#' event_burstiness(guastello, 1)
#' event_burstiness(guastello)
#' time_burstiness()
#'
time_burstiness <- function(times,
                            min_iet = NULL) {
  # Takes a vector of times and computes a burstiness coefficient.
  # Uses Kim & Jo (2016) burstiness calculation.
  if(is.numeric(ts) == FALSE) {         # Validate sequence
    cat("Invalid sequence type passed to time_burstiness().\n")
    cat("Did you mean to use event_burstiness()?\n")
    return(NULL)
  }

  if(length(ts) == 0) {                  # If ts has no entries
    cat("Empty sequence.\n")
    return(NULL)
  }

  if(is.null(min_iet) == TRUE) {         # Use Kim & Jo, eqn. 22
    # No minimum inter-event time.
    diff(ts, 1) ->                       # Interevent times
      iet_vec
    mean(iet_vec, na.rm = TRUE) ->       # mean(interevent times)
      mean_iet
    sd(iet_vec, na.rm = TRUE) ->         # sd(interevent times)
      sd_iet
    sd_iet / mean_iet ->                 # Coefficient of variation
      r
    length(code_pos_vec) ->              # Number of interevent times
      n
    (r * sqrt(n+1) - sqrt(n-1)) /        # Burstiness (K & J, eq. 22)
      (r * (sqrt(n+1) - 2) + sqrt(n-1)) ->
      B
    return(B)
  }

  if(is.null(min_iet) == FALSE) {            # Minimum inter-event time
    if(min_iet %% 1 == 0) {                  # Must be an integer!
      if(min_iet >= 1) {                     # Must be a positive integer!
        min_iet / length(ts) ->              # For K & J, eq. 28
          y_tilde
        diff(ts, 1) ->                       # Interevent times
          iet_vec
        mean(iet_vec, na.rm = TRUE) ->       # mean(interevent times)
          mean_iet
        sd(iet_vec, na.rm = TRUE) ->         # sd(interevent times)
          sd_iet
        sd_iet / mean_iet ->                 # Coefficient of variation
          r
        length(code_pos_vec) ->              # Number of interevent times
          n
        ((n - 2) * (r * sqrt(n+1) - (1 - n * y_tilde) * sqrt(n-1))) /
          (r * (n * sqrt(n+1) - 2*(n-1)) +
             (1 - n * y_tilde) * sqrt(n-1) * (n - 2 * sqrt(n+1))) ->
          B                                # Burstiness (K & J, eq. 28)
        return(B)
      }
      cat("Invalid minimum inter-event time encountered!/n") # Non-positive
      cat("Minimum inter-event time must be NULL or a positive integer./n")
      cat(paste(min_iet, "is invalid./n"))
      return(NULL)
    }
    cat("Invalid minimum inter-event time encountered!/n") # Non-integer
    cat("Minimum inter-event time must be NULL or a positive integer./n")
    cat(paste(min_iet, "is invalid./n"))
    return(NULL)
  }

  cat("Invalid data/n")            # Should never get here
  cat("Did you mean to use time_burstiness()?/n")
  return(NULL)
}

#' orbde
#'
#' @param <data_seq> character or integer time series
#' @keywords orbital decomposition
#' @description This function returns the result of orbital decomposition.
#' All of the metrics associated with orbital decomposition are returned in a
#' matrix with named columns.
#' @references Guastello, S. J. (n.d.). Orbital Decomposition: Identification
#' of Dynamical Patterns in Categorical Data. In S. J. Guastello & R. A. M.
#' Gregson (Eds.), Nonlinear Dynamical Systems Analysis for the Behavioral
#' Sciences Using Real Data (pp. 499–516). CRC Press.
#' Guastello, S. J., Hyde, T., & Odak, M. (1998). Symbolic dynamic patterns of
#' verbal exchange in a creative problem solving group. Nonlinear Dynamics,
#' Psychology, and Life Sciences, 2(1), 35–58.

#' @export
#' @examples
#' orbde(guastello)
#'
orbde <- function(data_seq) {
  # Returns a table  with the final orbital decomposition analysis:
  #  - String length, C
  #  - Trace of the transition matrix for the string length C, M^C
  #  - Topological entropy, H_T (this is better than H_S for these purposes;
  #    use minimum of H_T)
  #  - Shannon entropy, H_S
  #  - D_L - dimensionality via Lyapunov
  #  - chi-square
  #  - degrees of freedom, df
  #  - p-value, p
  #  - Number of possible sub-strengths of length C, N* = length + 1 - C
  #  - phi-square (phi-sq = chi-sq/N*)
  # This all follows Guastello, Chapter 21, in Guastello & Gregson, 2011,
  #  _Nonlinear Dynamical Systems Analysis for the Behavioral Sciences Using
  #  Real Data_.
  #
  # 52 codes should be sufficient for any application; it is the limit on the
  #  ORBDE software as well.
  require(data.table)  # Necessary to allow parallelizing to work
  # data.table::shift() is used here.
  require(glue)

  # First, recast everything into single character codes:
  unique(data_seq) -> uni_seq
  length(uni_seq) -> n_codes
  if(n_codes > 52) {
    cat("Cannot do more than 52 codes\n") # See below.
    return(NULL)
  }
  c(LETTERS, letters)[1:n_codes] ->
    uni_rep                         # c(LETTERS, letters) has 52 elements
  uni_seq -> names(uni_rep)
  uni_rep[data_seq] -> data_seq

  # For now, remove all missing data from the sequence. In the case where
  #  the missing data is due to uncoded data, this will cause problems,
  #  but so be it.
  if(any(is.na(data_seq)) == TRUE) {
    data_seq[-which(is.na(data_seq))] -> data_seq
  }

  unique(data_seq) -> uni_seq
  length(uni_seq) -> n_codes

  table(data_seq) -> freqs
  # May want this for the previous line:
  # table(data_seq, useNA = "always") -> freqs

  glue_collapse(data_seq) -> coll_seq
  # Begin processing data: For C = 1
  1 -> C
  length(data_seq) -> seq_len -> N_star

  shift(data_seq, 1) -> shifted_seq
  length(which(data_seq == shifted_seq)) -> recurs

  if(recurs > 0) {
    unique(data_seq[data_seq == shifted_seq]) -> repeats




    repeats[!is.na(repeats)] -> repeats




    length(repeats) -> trMC
    if (trMC > 0) {
      log2(trMC) / C -> H_T
    } else {
      -Inf -> H_T
    }
    exp(H_T) -> D_L
  } else {
    0 -> trMC
    -Inf -> H_T
    0 -> D_L
  }
  n_codes - 1 -> dof # For the singlets, this is true;

  table(data_seq) -> freqs  # I think this is a repeated command,
  #  but everything breaks if I remove it.
  freqs / seq_len -> p_obs_keep  # This gets used for C > 1 calculations
  rep(seq_len / n_codes,
      n_codes) -> F_exp
  freqs / seq_len -> p_obs_tab
  - 1 * sum(p_obs_tab * log(p_obs_tab)) -> H_S

  sum(freqs * log(freqs / F_exp)) * 2 ->
    chi_sq

  dchisq(chi_sq, dof) -> p
  chi_sq / N_star -> phi_sq

  data.frame("C" = 1,
             "trM" = trMC,
             "Ht" = H_T,
             "Dl" = D_L,
             "chi^2" = chi_sq,
             "df" = dof,
             "N*" = N_star,
             "Phi^2" = phi_sq,
             "Hs" = H_S,
             "p" = p) -> OD_tab

  # Processing for 1 < C < N_star; while trace(C^M) > 0
  seq_len -> N_star
  for(len in 1:(length(data_seq) / 2)) {    # This is the ORBDE default
    # Identify all recurrences of length C
    C + 1 -> C
    N_star - 1 -> N_star

    # I suspect that this next part is the part that takes time
    vector(mode = "character",
           length = N_star) -> current_seq
    C -> end_index
    # The old way: routine is >10 times slower with this:
    #    for (index in 1:N_star) {
    #      paste0(data_seq[index:end_index],
    #             collapse = "") -> current_seq[index]
    #      end_index + 1 -> end_index
    #    }

    # New way: routine is much faster than with the old way
    for(index in 1:N_star) {
      substr(coll_seq, index, end_index) -> current_seq[index]
      end_index + 1 -> end_index
    }

    shift(current_seq, C) -> shifted_seq

    # Can we do a loop here and skip the which?
    # E.g.: if(current_seq == shifted_seq) {...}
    which(current_seq == shifted_seq) -> repeat_nums
    # I doubt it is faster.

    if (length(repeat_nums) > 0) {
      length(repeat_nums) -> recurs

      unique(current_seq[repeat_nums]) -> repeats  # Which codes are repeated?



      repeats[!is.na(repeats)] -> repeats # Just to be sure




      length(repeats) -> trMC                   # How many repeated codes
      #  are there?
      log2(trMC) / C -> H_T                  # Topological entropy
      exp(H_T) -> D_L                        # Lyapunov dimension

      # The number of unique repeated codes:
      table(current_seq) -> repeated_codes -> F_obs_tab
      length(repeated_codes[repeated_codes > 1]) -> dof
      # Shannon entropy. Must do this before collapsing codes:
      sum((F_obs_tab / N_star) *
            (log(N_star/F_obs_tab))) -> H_S

      F_obs_tab[F_obs_tab > 1] -> F_rep_obs_tab  # Repeated codes

      # Frequency expected
      length(F_rep_obs_tab) -> n_rep
      rep(1, n_rep) -> F_exp           # For all the repeats and

      # Is there a faster way here?
      for (i in 1:n_rep) {
        for (j in 1:C) {
          substr(names(F_rep_obs_tab)[i], j, j) -> fn
          p_obs_keep[fn] * F_exp[i] -> F_exp[i]
        }
      }
      F_exp * N_star -> F_exp

      # Guastello eq. 21.6:
      #  chi-squared = 2 * sum( F_obs * ln(F_obs/F_expected) )
      sum(F_rep_obs_tab * log(F_rep_obs_tab / F_exp)) * 2 -> chi_sq
      abs(N_star - sum(F_rep_obs_tab)) -> singles
      # Need to do something for all the non-recurrent sequences, as they
      #  contribute here too.
      if(singles > 0.01 ) {        # 0.01 is capricious, but probably good
        chi_sq + 2 *singles * log(singles / (N_star - sum(F_exp))) -> chi_sq
      }

      # Now for p-value and phi-squared (Guastello eq 21.7):
      dchisq(chi_sq, dof) -> p
      chi_sq / N_star -> phi_sq

      if(!is.na(trMC)) {
        if(trMC != 0) {
          data.frame("C" = C,
                     "trM" = trMC,
                     "Ht" = H_T,
                     "Dl" = D_L,
                     "chi^2" = chi_sq,
                     "df" = dof,
                     "N*" = N_star,
                     "Phi^2" = phi_sq,
                     "Hs" = H_S,
                     "p" = p) -> OD_line
          rbind(OD_tab, OD_line) -> OD_tab
        }
      }
    }
  }

  if(OD_tab$trM[1] == 0) {
    OD_tab[-1,] -> OD_tab
  }

  return(OD_tab)
}

#' IPL_fit
#'
#' @param <ts> character vector
#' @keywords power-law
#' @description This function computes a power-law distribution fit to
#' @export
#' @examples
#' @references Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009).
#' Power-law distributions in empirical data. SIAM Review, 51(4), 661–703.
#' https://doi.org/10.1137/070710111
#' IPL_fit(guastello, C)
#'
IPL_fit <- function(ts,       # data sequence
                    C = 1) {  # subsequences length
                              # (from orbital decomposition, perhaps)
  # Create the subsequences
  require(igraph)
  length(ts) - C + 1 -> len
  if(len > 0) {
    vector(mode = "character", length = len) ->
      sub_seq
    for (seq_start in 1:len) {
      paste0(ts[seq_start:(seq_start + C - 1)],
             collapse = "") ->
        sub_seq[seq_start]
    }

    table(sub_seq) -> freqs
    length(unique(sub_seq)) -> num_codes

    # Compute the entropy; see note following this chunk
    # More

    # Fit the data with a power law
    igraph::fit_power_law(unname(freqs)) ->
      fit_ls
#    fit_ls[[2]] -> Shape
#    fit_ls[[6]] -> Fit

  }
  return(fit_ls)
}

#' makeMarkov
#'
#' This function returns a numeric sequence corresponding to a text
#' file, suitable for RQA use.
#' @param <x> text file name
#' @keywords clean text
#' @export
#' @examples
#' makeMarkov(ts)
#'
makeMarkov <- function(ts) {   # category vector
  # With some help from SO:
  # https://stackoverflow.com/questions/47329234/how-to-build-a-markovs-chain-transition-probability-matrix
  # The time series is either a single thing:
  #  "ATTCAACACATCCAGCCACATGCTCCGAG"
  # or multiple things:
  # c("A", "B", "A", "C")
  if(length(ts) == 1) {   # Needs a split
    ts <- unlist(strsplit(ts, split =""))
  }

  sort(unique(ts)) ->
    ts_unique

  matrix(0,
         ncol = length(ts_unique),
         nrow=length(ts_unique),
         dimnames = list(ts_unique, ts_unique)) ->
    markov

  for (i in 1:(length(ts) - 1)){
    #    ts == ts[i] ->
    #      rw
    #    ts == ts[i + 1] ->
    #      cl
    markov[ts[i], ts[i+1]] + 1 ->
      markov[ts[i], ts[i+1]]
  }
  markov <- markov / rowSums(markov)
  return(markov)
}

#' bootMarkov
#'
#' This function returns a numeric sequence corresponding to a text
#' file, suitable for RQA use.
#' @param <obs_adj_mat> observed adjacency matrix
#' @param <obs_node_freq> observed node frequencies
#' @param <directed> is the network directed (TRUE) or not (FALSE)
#' @param <loops> are there self-loops in the matrix?
#' @param <ci> confidence interval level
#' @param <R> number of replications
#' @keywords clean text
#' @description This function returns a matrix of the same dimension as
#' obs_adj_mat, with values of {-1,0,1}. Values of -1 (+1) indicate that the
#' corresponding entry in obs_adj_mat has an entry that is less (more) than
#' would be expected by the observed node frequencies and outside the
#' bootstrapped confidence interval matrix created from the observed node
#' frequencies.
#' @export
#' @examples
#' bootMarkov()
#'
bootMarkov <- function(obs_adj_mat,       # Observed adjacency matrix
                       obs_node_freq,     # Observed node frequencies
                       directed = TRUE,   # Directed network?
                       loops = TRUE,      # Loops in the matrix?
                       ci = 0.95,         # Confidence interval
                       R = 1e3) {         # Number of replications
  # Bootstrap whether the adj_mat entries fall within the
  #  ci ("black") or not ("blue" for less than, and "green" for
  #  greater than).
  # Updates to do:
  #  Work with non-directed networks
  #  match a degree distribution
  #
  length(obs_node_freq) -> num_nodes
  if(directed == TRUE) {
    sum(obs_adj_mat) -> num_edges      # This should be sum(obs_node_freq) - 1, I think
  }
  rep(1:num_nodes, obs_node_freq) -> node_pop

  matrix(0,
         ncol = num_nodes,
         nrow = num_nodes) -> result_mat
  array(0,
        dim = c(num_nodes,
                num_nodes,
                R)) -> samples_arr
  array(FALSE,
        dim = c(num_nodes,
                num_nodes,
                2)) -> quantile_arr

  for (i in 1:R) {
    replicate(num_edges,
              sample(node_pop,
                     2,
                     replace = loops)) -> edge_samp
    for (j in 1:num_edges) {
      samples_arr[edge_samp[1, j], edge_samp[2, j], i] + 1 ->
        samples_arr[edge_samp[1, j], edge_samp[2, j], i]
    }
  }

  # Now, to get the confidence intervals
  min_quant = (1 - ci) / 2
  max_quant = 1 - min_quant
  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (obs_adj_mat[i, j] < quantile(samples_arr[i, j, ],
                                       min_quant)) {
        -1 -> result_mat[i, j]
      } else {
        if (obs_adj_mat[i, j] > quantile(samples_arr[i, j, ],
                                         max_quant)) {
          1 -> result_mat[i, j]
        }
      }
    }
  }
  dimnames(obs_adj_mat) -> dimnames(result_mat)
  return(result_mat)
}

#' SSG
#'
#' This function returns a numeric sequence corresponding to a text
#' file, suitable for RQA use.
#' @param <x> text file name
#' @keywords clean text
#' @export
#' @examples
#' SSG(ts)
#'
SSG <- function(ts1,              # time series
                ts2,             # time series 2
                times = NULL,    # ts times - not yet implemented
                cats1 = NULL,
                cats2 = NULL) {  # category vector
  require("scales")

  if(length(unique(ts1)) > 10 | length(unique(ts2)) > 10) {
    cat("More than 10 unique codes; cannot be plotted.\nAre data continuous?\n")
    return(NULL)
  }

  if(is.null(cats1) == TRUE) {          # Set up the axes labels
    1:length(unique(ts1)) -> cat_num1
    sort(unique(ts1)) -> names(cat_num1)
  } else {
    1:length(cats1) -> cat_num1
    cats1 -> names(cat_num1)
  }

  if(is.null(cats2) == TRUE) {          # Set up the axes labels
    1:length(unique(ts2)) -> cat_num2
    sort(unique(ts2)) -> names(cat_num2)
  } else {
    1:length(cats2) -> cat_num2
    cats2 -> names(cat_num2)
  }

  cat_num1[ts1] ->
    seq_num1
  cat_num2[ts2] ->
    seq_num2

  jitter(x = seq_num1,
         factor = 1.5,
         amount = NULL) -> jt1
  jitter(seq_num2,
         factor = 1.5,
         amount = NULL) -> jt2

  if(is.null(times) == TRUE) {
    rep(1, length(jt1)) -> pt_sizes
  }
  if(is.null(times) == FALSE) {
    diff(times) -> pt_sizes
    sqrt(pt_sizes * 4 / max(pt_sizes)) ->
      pt_sizes
    if(length(pt_sizes) ==
       (length(jt1) - 1)) {
      c(pt_sizes, 0.1) -> pt_sizes            # Placeholder for the last time
    }
  }
  {
    plot(jt1, jt2,
         pch = 16,
         xaxt = "n",  # No numbers on x-axis
         yaxt = "n",
         xlab = "First time series",
         ylab = "Second time series",
         cex = pt_sizes)
    grid(length(cat_num1), length(cat_num2))
    axis(1,
         at = 1:length(cat_num1),
         #         at = seq(from = 0.5,
         #                  to = length(cat_num1)+0.5,
         #                  by = 1),
         lab = names(cat_num1))
    axis(2,
         at = 1:length(cat_num2),
         #         at = seq(from = 0.5,
         #                  to = length(cat_num1)+0.5,
         #                  by = 1),
         lab = names(cat_num2))
    # xlab, ylab, change number to categories, change boxes,
    #  gridlines
    #    x <- seq(-4, 4, len = 101)
    #    plot(x,sin(x),type="l",xaxt="n", col="red",
    #                                    xlab=expression(paste("Phase Angle ",phi)),
    #                                    ylab=expression("sin "*phi))
    #    axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
    #         lab = expression(-pi, -pi/2, 0, pi/2, pi))
    lines(jt1, jt2,
          col = alpha("blue", 0.60),
          lwd = 0.7)
  }
  return(NULL)
}

#' SSG_metrics
#'
#' This function returns a numeric sequence corresponding to a text
#' file, suitable for RQA use.
#' @param <x> text file name
#' @keywords clean text
#' @export
#' @examples
#' SSG_metrics(ts)
#'
SSG_metrics <- function(ts1,             # time series
                        ts2,             # time series 2
                        times = NULL,    # ts times - not yet implemented
                        cats1 = NULL,
                        cats2 = NULL) {  # category vector

  if(is.null(cats1) == TRUE) {          # Set up the axes labels
    1:length(unique(ts1)) -> cat_num1
    sort(unique(ts1)) -> names(cat_num1)
  } else {
    1:length(cats1) -> cat_num1
    cats1 -> names(cat_num1)
  }

  if(is.null(cats2) == TRUE) {          # Set up the axes labels
    1:length(unique(ts2)) -> cat_num2
    sort(unique(ts2)) -> names(cat_num2)
  } else {
    1:length(cats2) -> cat_num2
    cats2 -> names(cat_num2)
  }

  cat_num1[ts1] -> seq_num1
  cat_num2[ts2] -> seq_num2

  jitter(ts1) -> jt1
  jitter(ts2) -> jt2

  if(is.null(times) == TRUE) {
    plot(jt1, jt2,
         pch = 16)
    lines(jt1, jt2)
  }


}

