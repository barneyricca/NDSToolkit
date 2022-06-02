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
#' @param <x> text file name
#' @keywords clean text
#' @export
#' @examples
#' bootMarkov(ts)
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
