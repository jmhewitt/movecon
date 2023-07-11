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

# build statespace search util
search = build_statespace_search(statespace = statespace_constrained)

#
# find the closest match for a valid location
#

# get grid indices for a valid location
test_ind = c(
  easting_ind = unname(valid_locs[218,'row']), 
  northing_ind = unname(valid_locs[218,'col'])
)

# map coordinates to grid
res = nearest_location_in_domain(
  statespace_search = search,
  easting = eastings[test_ind['easting_ind']],
  northing = northings[test_ind['northing_ind']]
)

# extract states associated with the same mapped location
states = states_at_nearest_location_in_domain(
  statespace_search = search,
  easting = eastings[test_ind['easting_ind']],
  northing = northings[test_ind['northing_ind']]
)

# verify we recover the valid location
expect_equal(
  c(eastings[test_ind['easting_ind']], northings[test_ind['northing_ind']]),
  as.numeric(res[c('easting', 'northing')])
)

# verify the states are all associated with the location
for(ind in seq_along(states)) {
  expect_identical(states[[ind]]$location, res)
}

# verify we get all possible states for this location
expect_identical(
  sapply(states, function(s) s$last_movement_direction),
  c('north', 'south', 'east', 'west')
)

#
# find the closest match for something near a valid location
#

# get grid indices for a valid location
test_ind = c(
  easting_ind = unname(valid_locs[218,'row']), 
  northing_ind = unname(valid_locs[218,'col'])
)

# map perturbed coordinates to grid
res = nearest_location_in_domain(
  statespace_search = search,
  easting = eastings[test_ind['easting_ind']] + .1,
  northing = northings[test_ind['northing_ind']] - .1
)

# verify we recover the valid location
expect_equal(
  c(eastings[test_ind['easting_ind']], northings[test_ind['northing_ind']]),
  as.numeric(res[c('easting', 'northing')])
)

#
# find the closest match for an invalid location
#

# get grid indices for an invalid location
test_ind = c(
  easting_ind = unname(invalid_locs[218,'row']), 
  northing_ind = unname(invalid_locs[218,'col'])
)

# map to closest good coordinates in domain
res = nearest_location_in_domain(
  statespace_search = search,
  easting = eastings[test_ind['easting_ind']],
  northing = northings[test_ind['northing_ind']]
)

expect_failure(
  expect_equal(
    c(eastings[test_ind['easting_ind']], northings[test_ind['northing_ind']]),
    as.numeric(res[c('easting', 'northing')])
  )  
)
