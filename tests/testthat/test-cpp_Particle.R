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

# verify all visited locations satisfy linear constraint
expect_true(
  all(
    linear_constraint %*% sapply(path, function(x) x$location$covariates) >= 0
  )
)

# # view path
# scale = .0001
# path_range = list(
#   xlim = range(sapply(path, function(x) x$location$easting)) * 
#     (1 + c(-scale, scale)),
#   ylim = range(sapply(path, function(x) x$location$northing)) * 
#     (1 + c(-scale, scale))
# )
# # base map
# image(
#   x = dat > band1_avg,
#   band = 1,
#   axes = TRUE,
#   xlim = path_range$xlim,
#   ylim = path_range$ylim
# )
# # path
# for(i in 2:length(path)) {
#   lines(
#     x = c(path[[i-1]]$location$easting, path[[i]]$location$easting),
#     y = c(path[[i-1]]$location$northing, path[[i]]$location$northing),
#     col = 'blue'
#   )
# }

#
# test: movement explores a region
#

# choose a location with a large range for movement
starting_location = c(296000, 9119500)

# crudely map the coordinate onto the grid
starting_coord_ind = which.min(colSums((t(coords) - starting_location)^2))
starting_coord_inds = c(
  easting_ind = which.min(abs(coords$x[starting_coord_ind] - eastings)),
  northing_ind = which.min(abs(coords$y[starting_coord_ind] - northings))
)

# simulate movement
path = Test__Particle_Steps(
  statespace = statespace_constrained, 
  last_movement_direction = 'west', 
  easting_ind = starting_coord_inds['easting_ind'] - 1, 
  northing_ind = starting_coord_inds['northing_ind'] - 1, 
  directional_persistence = 1.25, 
  beta = rep(0, nrow(covariates)),
  delta = .9, 
  nsteps = 1000
)

# verify all visited locations satisfy linear constraint
expect_true(
  all(
    linear_constraint %*% sapply(path, function(x) x$location$covariates) >= 0
  )
)

# # view path
# scale = .0001
# path_range = list(
#   xlim = range(sapply(path, function(x) x$location$easting)) * 
#     (1 + c(-scale, scale)),
#   ylim = range(sapply(path, function(x) x$location$northing)) * 
#     (1 + c(-scale, scale))
# )
# # base map
# image(
#   x = dat > band1_avg,
#   band = 1,
#   axes = TRUE,
#   xlim = path_range$xlim,
#   ylim = path_range$ylim
# )
# # path
# for(i in 2:length(path)) {
#   lines(
#     x = c(path[[i-1]]$location$easting, path[[i]]$location$easting),
#     y = c(path[[i-1]]$location$northing, path[[i]]$location$northing),
#     col = 'blue'
#   )
# }

#
# test: transition distribution is sampled correctly
#

set.seed(2023)

# set movement parameters
last_movement_direction = 'west'
beta = rnorm(n = nrow(covariates), sd = .001)
directional_persistence = 1.25
delta = .9

# identify neighbors of the starting state
state = extract_statespace_state(
  statespace = statespace_constrained,
  last_movement_direction = last_movement_direction, 
  easting_ind = starting_coord_inds['easting_ind'] - 1, 
  northing_ind = starting_coord_inds['northing_ind'] - 1
)

# compute the transition rate from the starting state
tx_rate = Test__Location_Based_Movement_Transition_Rate(
  statespace = statespace_constrained,
  last_movement_direction = last_movement_direction, 
  easting_ind = starting_coord_inds['easting_ind'] - 1, 
  northing_ind = starting_coord_inds['northing_ind'] - 1,
  beta = beta
)

# compute transition probabilities to neighbors
probs = Test__Directional_Transition_Probabilities(
  statespace = statespace_constrained,
  last_movement_direction = last_movement_direction, 
  easting_ind = starting_coord_inds['easting_ind'] - 1, 
  northing_ind = starting_coord_inds['northing_ind'] - 1, 
  directional_persistence = directional_persistence
)
names(probs) = sapply(state$to, function(x) x$last_movement_direction)

# uniformized transition distribution
target_probs = c('none' = 1 - tx_rate * delta, tx_rate * delta * probs)

# sample transition events
move_samples = replicate(
  n = 5e3, 
  expr = {
    # simulate movement
    path = Test__Particle_Steps(
      statespace = statespace_constrained, 
      last_movement_direction = last_movement_direction, 
      easting_ind = starting_coord_inds['easting_ind'] - 1, 
      northing_ind = starting_coord_inds['northing_ind'] - 1, 
      directional_persistence = directional_persistence, 
      beta = beta,
      delta = delta, 
      nsteps = 1
    )
    
    # determine what movement occurred
    if(identical(path[[1]]$location, path[[2]]$location)) {
      'none'
    } else {
      path[[2]]$last_movement_direction
    }
  }
)

# MC estimate of transition distribution
empirical_probs = table(move_samples) / length(move_samples)

# compare sampler to theoretical transition probabilities
expect_equivalent(
  as.numeric(empirical_probs[names(target_probs)]),
  target_probs, 
  tolerance = .015 # relatively high tolerance b/c MC samples are small
)
