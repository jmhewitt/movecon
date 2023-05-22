require(stars)

# load test raster
dat = read_stars(system.file("tif/L7_ETMs.tif", package = "stars"))

# get grid definition
coords = st_coordinates(dat)
eastings = unique(coords$x)
northings = unique(coords$y)

# extract covariates s.t. all covariate data is grouped by location
covariates = matrix(
  data = apply(X = dat[['L7_ETMs.tif']], MARGIN = 1:2, FUN = identity), 
  nrow = dim(dat)[3]
)

statespace = build_statespace(
  eastings = eastings, northings = northings, covariates = covariates
)

test_ind = c(easting = 70, northing = 70)

location = extract_statespace_location(
  statespace = statespace, 
  easting_ind = test_ind['easting'] - 1, 
  northing_ind = test_ind['northing'] - 1
)

#
# locations are indexed correctly and information can be retrieved
#

expect_equal(
  location$covariates, 
  dat[["L7_ETMs.tif"]][test_ind['easting'], test_ind['northing'], ]
)

expect_equal(location$easting, eastings[test_ind['easting']])
expect_equal(location$northing, northings[test_ind['northing']])


#
# invalid states do not exist in the C++ structures
#

expect_error(
  extract_statespace_state(
    statespace = statespace, 
    last_movement_direction = 'south',
    easting_ind = 0, 
    northing_ind = 0
  )
)


#
# states are indexed correctly and information can be retrieved
#


state = extract_statespace_state(
  statespace = statespace, 
  last_movement_direction = 'west',
  easting_ind = 0, 
  northing_ind = 0
)

state_east = extract_statespace_state(
  statespace = statespace, 
  last_movement_direction = 'east',
  easting_ind = 1, 
  northing_ind = 0
)

state_south = extract_statespace_state(
  statespace = statespace, 
  last_movement_direction = 'south',
  easting_ind = 0, 
  northing_ind = 1
)

state_ids = sapply(state$to, function(x) x$last_movement_direction)

expect_equal(
  state$to[[which(state_ids == 'east')]]$location, 
  state_east$location
)

expect_equal(
  state$to[[which(state_ids == 'south')]]$location, 
  state_south$location
)

#
# a state's last_movement_direction makes sense with the actual location pairs
#

for(loc_ind in length(state$to)) {
  dx = state$to[[loc_ind]]$location$easting - state$location$easting
  dy = state$to[[loc_ind]]$location$northing - state$location$northing
  expect_true(
    switch(
      EXPR = state$to[[loc_ind]]$last_movement_direction,
      'south' = dy < 0,
      'north' = dy > 0,
      'east' = dx > 0,
      'west' = dx < 0
    )
  )
}
