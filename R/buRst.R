#' event_burstiness
#'
#' This function returns the burstiness coefficient for each of the
#' unique codes in a series. The series must be event-based; see
#' time_burstiness() to find the burstiness of codes in a time-series.
#' @param <ts> character (or logicial) vector
#' @param <min_iet> minimum interevent sequence spacing
#' @keywords burstiness
#' @export
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
#'
#' This function returns the burstiness coefficient for a series of events.
#' The series must be time-based; see event_burstiness() to find the
#' burstiness of codes in a event-series.
#' @param <ts> series of numeric times of events
#' @param <min_iet> minimum interevent time
#' @keywords burstiness
#' @export
#' @examples
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
