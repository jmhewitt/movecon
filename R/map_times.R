#' Discretize observation times to a grid of discrete times
#' 
#' Discretization maps the input times \code{t} to the closest matching time 
#' in the regular grid of discrete times implied by \code{t0}, \code{tf}, and
#' \code{dt}.
#'
#' @param t0 \code{POSIXct} or \code{numeric} for start of time grid
#' @param tf \code{POSIXct} or \code{numeric} for end of time grid
#' @param dt Number of seconds between time steps
#' @param t Vector of times to map to grid
#' 
#' @export
#' 
#' @examples map_times(t0 = 0, tf = 1, dt = .001, t = c(.123, .234, .345))
#' @examples map_times(
#'     t0 = round(Sys.time()), 
#'     tf = round(Sys.time()) + 3600, 
#'     dt = 1, 
#'     t = Sys.time() + c(20, 60, 90)
#'  )
#'
map_times = function(t0, tf, dt, t) {
  
  # convert all inputs to numeric (or seconds)
  t0 = as.numeric(t0)
  tf = as.numeric(tf)
  dt = as.numeric(dt)
  t = as.numeric(t)
  
  # number of timesteps
  nt = (tf - t0)/dt
  
  if(t0 >= tf) {
    stop('t0 must be before tf')
  }
  
  if(nt != round(nt)) {
    stop('t0, tf, and dt must imply a whole number of time steps')
  }
  
  if(any(unique(diff(order(t))) != 1)) {
    stop('Observation times t must be in increasing order')
  }
  
  if(length(unique(t)) != length(t)) {
    stop('Observation times t must be unique')
  }
  
  if(t[1] < t0) {
    stop('t0 must be before first observation time t')
  }
  
  if(tf < tail(t,1)) {
    stop('tf must be after last observation time t')
  }
  
  # discrete time grid
  tseq = seq(from = t0, to = tf, by = dt)
  
  # map time steps by rounding to closest match
  tinds = sapply(t, function(t) which.min(abs(t - tseq)))
  
  if(any(unique(diff(order(tinds))) != 1)) {
    stop('Discretized times must be in increasing order')
  }
  
  if(length(unique(tinds)) != length(tinds)) {
    stop('Discretized times must be unique')
  }
  
  # package results
  list(
    tgridded = tinds, 
    ngrid = nt
  )
}
