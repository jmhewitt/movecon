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

statespace = build_statespace(
  eastings = eastings, northings = northings, covariates = covariates, 
  linear_constraint = rep(0, nrow(covariates))
)

#
# test: no directional persistence implies uniform transition
#

probs = Test__Directional_Transition_Probabilities(
  statespace = statespace, 
  last_movement_direction = 'west', 
  easting_ind = 0, 
  northing_ind = 30, 
  directional_persistence = 0
)

# yields a probability distribution
expect_equal(sum(probs), 1)

# uniform
expect_equal(probs, rep(1/length(probs), length(probs)))

#
# test: strong directional persistence behaves as expected
#

# specify state and associated
test_inds = c(easting_ind = 100, northing_ind = 30)
last_movement_direction = 'west'
opposite_movement_direction = 'east'
orthogonal_movement_directions = c('north', 'south')

state = extract_statespace_state(
  statespace = statespace, 
  last_movement_direction = last_movement_direction,
  easting_ind = test_inds['easting_ind'], 
  northing_ind =test_inds['northing_ind']
)

probs = Test__Directional_Transition_Probabilities(
  statespace = statespace, 
  last_movement_direction = last_movement_direction, 
  easting_ind = test_inds['easting_ind'], 
  northing_ind =test_inds['northing_ind'],
  directional_persistence = 1
)

names(probs) = sapply(state$to, function(x) x$last_movement_direction)

# yields a probability distribution
expect_equal(sum(probs), 1)

# preference along direction of movement/avoid opposite direction
expect_gt(probs[last_movement_direction], probs[opposite_movement_direction])

# no preference in orthogonal direction
expect_equivalent(
  probs[orthogonal_movement_directions[1]],
  probs[orthogonal_movement_directions[2]]
)

#
# test: investigate caching; it appears that caching is slower than computing
#

if(FALSE) {
  
  # repeated evaluation of directional preference computations with caching
  cache_times = microbenchmark::microbenchmark(
    Test__Directional_Transition_Probability_Cache(
      statespace = statespace, 
      last_movement_direction = last_movement_direction, 
      easting_ind = test_inds['easting_ind'], 
      northing_ind =test_inds['northing_ind'],
      directional_persistence = 1, 
      reps = 1000
    ), 
    times = 1e3
  )
  
  # repeated evaluation of directional preference computations without caching
  no_cache_times = microbenchmark::microbenchmark(
    Test__Directional_Transition_Probability_No_Cache(
      statespace = statespace, 
      last_movement_direction = last_movement_direction, 
      easting_ind = test_inds['easting_ind'], 
      northing_ind =test_inds['northing_ind'],
      directional_persistence = 1, 
      reps = 1000
    ),
    times = 1e3
  )
  
  # ratio shows that caching is slower than direct evaluation
  mean(cache_times$time)/mean(no_cache_times$time)
}

#
# test: transition rates are evaluated as expected
#

# random covariates, for demonstration purposes
set.seed(2023)
beta = rnorm(n = nrow(covariates), sd = .01)

# manually compute transition rate
location = extract_statespace_location(
  statespace = statespace, 
  easting_ind = test_inds['easting_ind'], 
  northing_ind =test_inds['northing_ind']
)
tx_rate_manual = exp(sum(location$covariates * beta))
  
# compute transition rate through C++ structures
tx_rate = Test__Location_Based_Movement_Transition_Rate(
  statespace = statespace, 
  last_movement_direction = last_movement_direction, 
  easting_ind = test_inds['easting_ind'], 
  northing_ind =test_inds['northing_ind'],
  beta = beta
)

expect_equal(tx_rate_manual, tx_rate)
