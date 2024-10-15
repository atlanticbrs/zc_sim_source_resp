### https://github.com/animaltags/tagtools_r/blob/main/R/find_dives.R
### https://github.com/animaltags/tagtools_r/blob/main/R/fir_nodelay.R
# GPL >=3

require(signal)

#' Find time cues for dives
#'
#' This function is used to find the time cues for the start and end of either dives in a depth record or flights in an altitude record.
#' @param p A depth or altitude time series (a sensor data list or  a vector) in meters.
#' @param sampling_rate The sampling rate of the sensor data in Hz (samples per second).
#' @param mindepth The threshold in meters at which to recognize a dive or flight. Dives shallow or flights lower than mindepth will be ignored.
#' @param surface (optional) The threshold in meters at which the animal is presumed to have reached the surface. Default value is 1. A smaller value can be used if the dive/altitude data are very accurate and you need to detect shallow dives/flights.
#' @param findall (optional) When TRUE, forces the algorithm to include incomplete dives at the start and end of the record. Default is FALSE which only recognizes complete dives.
#' @return dives is a data frame with one row for each dive/flight found. The columns of dives are: start (time in seconds of the start of each dive/flight), end (time in seconds of the start of each dive/flight), max (maximum depth/altitude reached in each dive/flight), tmax	(time in seconds at which the animal reaches the max depth/altitude).
#' @export
#' @examples
#' BW <- beaked_whale
#' dives <- find_dives(p = BW$P$data, 
#' sampling_rate = BW$P$sampling_rate, 
#' mindepth = 25, surface = 5, 
#' findall = FALSE)

find_dives <- function(p, mindepth, sampling_rate = NULL, surface = 1, findall = 0) {
  if (nargs() < 2) {
    stop("inputs for p and mindepth are required")
  }
  if (is.list(p)) {
    sampling_rate <- p$sampling_rate
    p <- p$data
    if (is.null(p)) {
      stop("p cannot be an empty vector")
    }
  } else {
    if (nrow(p) == 1) {
      p <- t(p)
    }
    if (is.null(sampling_rate)) {
      stop("sampling_rate is required when p is a vector")
    }
  }
  
  searchlen <- 20 # how far to look in seconds to find actual surfacing
  dpthresh <- 0.25 # vertical velocity threshold for surfacing
  dp_lp <- 0.25 # low-pass filter frequency for vertical velocity
  # find threshold crossings and surface times
  # hack for case where first depth obs is > mindepth
  if (p[1] > mindepth){p[1] <- mindepth - 0.25}
  
  tth <- which(diff(p > mindepth) > 0)
  tsurf <- which(p < surface)
  ton <- 0 * tth
  toff <- ton
  k <- 0
  empty <- integer(0)
  # sort through threshold crossings to find valid dive start and end points
  for (kth in 1:length(tth)) {
    if (all(tth[kth] > toff)) {
      ks0 <- which(tsurf < tth[kth])
      ks1 <- which(tsurf > tth[kth])
      if (findall || ((!identical(ks0, empty)) & (!identical(ks1, empty)))) {
        k <- k + 1
        if (identical(ks0, empty)) {
          ton[k] <- 1
        } else {
          ton[k] <- max(tsurf[ks0])
        }
        if (identical(ks1, empty)) {
          toff[k] <- length(p)
        } else {
          toff[k] <- min(tsurf[ks1])
        }
      }
    }
  }
  # truncate dive list to only dives with starts and stops in the record
  ton <- ton[1:k]
  toff <- toff[1:k]
  # filter vertical velocity to find actual surfacing moments
  n <- round(4 * sampling_rate / dp_lp)
  dp <- fir_nodelay(
    matrix(c(0, diff(p)), ncol = 1) * sampling_rate,
    n, dp_lp / (sampling_rate / 2)
  )
  # for each ton, look back to find last time whale was at the surface
  # for each toff, look forward to find next time whale is at the surface
  dmax <- matrix(0, length(ton), 2)
  for (k in 1:length(ton)) {
    ind <- ton[k] + (-round(searchlen * sampling_rate):0)
    ind <- ind[which(ind > 0)]
    ki <- max(which(dp[ind] < dpthresh))
    if (length(ki) == 0 | is.infinite(ki)) {
      ki <- 1
    }
    ton[k] <- ind[ki]
    ind <- toff[k] + (0:round(searchlen * sampling_rate))
    ind <- ind[which(ind <= length(p))]
    ki <- min(which(dp[ind] > -dpthresh))
    if (length(ki) == 0 | is.infinite(ki)) {
      ki <- 1
    }
    toff[k] <- ind[ki]
    dm <- max(p[ton[k]:toff[k]])
    km <- which.max(p[ton[k]:toff[k]])
    dmax[k, ] <- c(dm, ((ton[k] + km - 1) / sampling_rate))
  }
  # assemble output
  t0 <- cbind(ton, toff)
  t1 <- t0 / sampling_rate
  t2 <- dmax
  dmat <- cbind(t1, t2)
  dmat <- matrix(dmat[stats::complete.cases(dmat)], byrow = FALSE, ncol = 4)
  dives <- data.frame(
    start = dmat[, 1], end = dmat[, 2],
    max = dmat[, 3], tmax = dmat[, 4]
  )
  return(dives)
}

#' Delay-free filtering
#'
#' This function is used to gather a delay-free filtering using a linear-phase (symmetric) FIR filter followed by group delay correction. Delay-free filtering is needed when the relative timing between signals is important e.g., when integrating signals that have been sampled at different rates.
#' @param x The signal to be filtered. It can be multi-channel with a signal in each column, e.g., an acceleration matrix. The number of samples (i.e., the number of rows in x) must be larger than the filter length, n.
#' @param n The length of symmetric FIR filter to use in units of input samples (i.e., samples of x). The length should be at least 4/fc. A longer filter gives a steeper cut-off.
#' @param fc The filter cut-off frequency relative to sampling_rate/2=1. If a single number is given, the filter is a low-pass or high-pass. If fc is a vector with two numbers, the filter is a bandpass filter with lower and upper cut-off frequencies given by fc(1) and fc(2). For a bandpass filter, n should be at least 4/fc(1) or 4/diff(fc) whichever is larger.
#' @param qual An optional qualifier determining if the filter is: "low" for low-pass (the default value if fc has a single number), or "high" for high-pass. Default is "low".
#' @param return_coefs Logical. Return filter coefficients instead of filtered signal? If TRUE, the function will return the FIR filter coefficients instead of the filtered signal. Default is FALSE.
#' @export
#' @return If return_coefs is FALSE (the default), \code{fir_nodelay()} returns the filtered signal (same size as x). If return_coefs is TRUE, returns the vector of filter coefficients only.
#' @note The filter is generated by a call to \code{\link{fir1}}: \code{h <- fir1(n, fc, qual)}.
#' @note h is always an odd length filter even if n is even. This is needed to ensure that the filter is both symmetric and has a group delay which is an integer number of samples.
#' @note The filter has a support of n samples, i.e., it uses n samples from x to compute each sample in y.
#' @note The input samples used are those from n/2 samples before to n/2 samples after the sample number being computed. This means that samples at the start and end of the output vector y need input samples before the start of x and after the end of x. These are faked by reversing the first n/2 samples of x and concatenating them to the start of x. The same trick is used at the end of x. As a result, the first and last n/2 samples in y are untrustworthy. This initial condition problem is true for any filter but the FIR filter used here makes it easy to identify precisely which samples are unreliable.
#' @examples 
#' x <- sin(t(2 * pi * 0.05 * (1:100)) +
#'   t(cos(2 * pi * 0.25 * (1:100))))
#' Y <- fir_nodelay(x = x, n = 30, fc = 0.2, qual = "low")
#' plot(c(1:length(x)), x,
#'   type = "l", col = "grey42",
#'   xlab = "index", ylab = "input x and output y"
#' )
#' lines(c(1:length(Y)), Y, lwd = 2)
#'
fir_nodelay <- function(x, n, fc, qual = "low", return_coefs = FALSE) {
  # input checking
  # ================================================================
  # make sure x is a column vector or matrix
  if (!(sum(class(x) %in% c("matrix", "vector")))) {
    x <- as.matrix(x)
  }
  if (is.vector(x)) x <- as.matrix(x, nrow = length(x))
  
  # in case of multi-channel data, make sure matrix rows are samples and columns are channels
  if (dim(x)[2] > dim(x)[1]) x <- t(x)
  
  # make sure n is even to ensure an integer group delay
  n <- floor(n / 2) * 2
  
  
  # generate fir filter
  # ============================================================
  h <- signal::fir1(n = n, w = fc, type = qual)
  
  if (return_coefs) {
    return(h)
  } else { # carry out filtering
    
    # append fake samples to start and end of x to absorb filter delay
    # (output from these will be removed before returning result to user)
    nofsampling_rate <- floor(n / 2)
    top_pad <- matrix(x[nofsampling_rate:2, ], ncol = ncol(x))
    bot_pad <- matrix(x[(nrow(x) - 1):(nrow(x) - nofsampling_rate), ], ncol = ncol(x))
    x_pad <- rbind(top_pad, x, bot_pad)
    
    # filter the signal
    # ============================================================
    # apply filter to padded signal
    y <- apply(x_pad, MARGIN = 2, FUN = signal::filter, filt = h, nrow = nrow(x_pad))
    
    # account for filter offset (remove padding)
    y <- y[n - 1 + (1:nrow(x)), ]
    
    return(y)
  }
}