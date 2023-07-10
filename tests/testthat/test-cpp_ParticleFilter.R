require(stars)

# load test raster
dat = read_stars(system.file("tif/L7_ETMs.tif", package = "stars"))

# get grid definition
coords = unique(st_coordinates(dat)[, c('x', 'y')])
eastings = unique(coords$x)
northings = unique(coords$y)

# extract covariates s.t. all covariate data is grouped by location
covariates = matrix(
  data = apply(X = dat[['L7_ETMs.tif']], MARGIN = 1:2, FUN = identity), 
  nrow = dim(dat)[3]
)

# add an intercept to the covariates
covariates = rbind(1, covariates)

#
# build constrained domain
#

band1_avg = mean(dat$L7_ETMs.tif[,,1])
linear_constraint = c(-band1_avg, 1, rep(0, nrow(covariates)-2))

statespace_constrained = build_statespace(
  eastings = eastings, northings = northings, covariates = covariates, 
  # exclude locations whose band 1 value is below average
  linear_constraint = linear_constraint
)

# identify locations that will be un/defined
valid_locs = which(dat$L7_ETMs.tif[,,1] >= band1_avg, arr.ind = TRUE)
invalid_locs = which(dat$L7_ETMs.tif[,,1] < band1_avg, arr.ind = TRUE)

#
# test: movement stays within a constrained space
#

# identify locations near the the linear constraint boundary
contours = st_contour(dat, breaks = band1_avg)

# pick a location: geometry[[break]][[contour_ind]][[1]][location, coordinate]
boundary_location = contours$geometry[[1]][[100]][[1]][1,]

# crudely map the coordinate onto the grid
boundary_coord_ind = which.min(colSums((t(coords) - boundary_location)^2))
boundary_coord_inds = c(
  easting_ind = which.min(abs(coords$x[boundary_coord_ind] - eastings)),
  northing_ind = which.min(abs(coords$y[boundary_coord_ind] - northings))
)

# simulate movement
path = Test__Particle_Steps(
  statespace = statespace_constrained, 
  last_movement_direction = 'west', 
  easting_ind = boundary_coord_inds['easting_ind'] - 1, 
  northing_ind = boundary_coord_inds['northing_ind'] - 1, 
  directional_persistence = 0, 
  beta = rep(0, nrow(covariates)),
  delta = .9, 
  nsteps = 1000
)

# exact likelihood for movement
ll_exact = function(directional_persistence) {
  sapply(directional_persistence, function(directional_persistence) {
  sum(sapply(seq_along(path)[-1], function(ind) {
    
    # determine where movement was from
    from_inds = c(
      easting_ind = which.min(abs(path[[ind-1]]$location$easting - eastings)),
      northing_ind = which.min(abs(path[[ind-1]]$location$northing - northings))
    )
    
    # determine uniformized transition rate
    tx_rate = Test__Location_Based_Movement_Transition_Rate(
      statespace = statespace_constrained, 
      last_movement_direction = path[[ind-1]]$last_movement_direction,
      easting_ind = from_inds['easting_ind'] - 1, 
      northing_ind = from_inds['northing_ind'] - 1, 
      beta = rep(0, nrow(covariates))
    )
    
    # determine transition probabilities
    probs = Test__Directional_Transition_Probabilities(
      statespace = statespace_constrained, 
      last_movement_direction = path[[ind-1]]$last_movement_direction, 
      easting_ind = from_inds['easting_ind'] - 1, 
      northing_ind = from_inds['northing_ind'] - 1, 
      directional_persistence = directional_persistence
    )
    
    # associate transition probabilities with movements
    names(probs) = sapply(path[[ind-1]]$to, function(x) x$last_movement_direction)
    
    # extract log-prob for transition
    if(identical(path[[ind]]$location, path[[ind-1]]$location)) {
      # self-transition
      log(tx_rate * .1)
    } else {
      # transition
      log(probs[path[[ind]]$last_movement_direction]) + log(tx_rate * .9)
    }
  }))
  })
}
  
  



directional_persistence_seq = seq(from = -1, to = 1, length.out = 10)

ll_seq_exact = ll_exact(directional_persistence_seq)

plot(directional_persistence_seq, ll_seq_exact)

tick = Sys.time()
ll_seq = sapply(directional_persistence_seq, function(directional_persistence) {
  # approximate likelihood for movement
  list(Test__Particle_Filter_Likelihood(
    eastings = sapply(path, function(x) x$location$easting), 
    northings = sapply(path, function(x) x$location$northing), 
    semi_majors = rep(.1, length(path)),  
    semi_minors = rep(.1, length(path)), 
    orientations = rep(0, length(path)), 
    statespace = statespace_constrained, 
    last_movement_direction = 'west', 
    easting_ind = boundary_coord_inds['easting_ind'] - 1, 
    northing_ind = boundary_coord_inds['northing_ind'] - 1, 
    nparticles = 1e4, 
    directional_persistence = directional_persistence, 
    beta = rep(0, nrow(covariates)), 
    delta = .9
  ))
})
tock = Sys.time()

tock - tick
# Time difference of 1.564818 mins for 10 evaluations
# Time difference of 16.9716 secs


plot(directional_persistence_seq, ll_seq)



lp_seq = exp(ll_seq - log_sum(ll_seq))
lp_seq = lp_seq / sum(lp_seq)

plot(directional_persistence_seq, (lp_seq))

sum(directional_persistence_seq * lp_seq)
