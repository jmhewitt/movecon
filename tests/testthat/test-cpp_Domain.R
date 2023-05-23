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
# test: invalid states do not exist in the C++ structures
#

expect_error(
  extract_statespace_state(
    statespace = statespace, 
    last_movement_direction = 'south',
    easting_ind = 0, 
    northing_ind = 0
  )
)

test_inds = rbind(
  c(easting = 70, northing = 70),
  c(easting = 1, northing = 1)
)

#
# test: locations are indexed correctly and covariates can be retrieved
#

apply(test_inds, 1, function(test_ind) {
  
  location = extract_statespace_location(
    statespace = statespace, 
    easting_ind = test_ind['easting'] - 1, 
    northing_ind = test_ind['northing'] - 1
  )
  
  expect_equal(
    location$covariates[-1], # remove intercept
    dat[["L7_ETMs.tif"]][test_ind['easting'], test_ind['northing'], ]
  )
  
  expect_equal(location$easting, eastings[test_ind['easting']])
  expect_equal(location$northing, northings[test_ind['northing']])
  
  0
})

#
# test: states are indexed correctly and information can be retrieved
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

apply(test_inds, 1, function(test_ind) {
  
  for(last_movement_direction in c('north', 'east', 'south', 'west')) {
    
    # try to extract state
    state = tryCatch(
      extract_statespace_state(
        statespace = statespace, 
        last_movement_direction = last_movement_direction,
        easting_ind = test_ind['easting'] - 1, 
        northing_ind = test_ind['northing'] - 1
      ), 
      error = function(e) { }
    )
    
    # the loop does not guarantee all states will be defined, so skip when 
    # needed since this test goal is not to check state definition, but to check
    # connections
    if(is.null(state))
      next
    
    #
    # a state's last_movement_direction makes sense with the location pairs
    #
    
    for(loc_ind in seq_along(state$to)) {
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
    
  }
  
  0
})

#
# build constrained domain
#

band1_avg = mean(dat$L7_ETMs.tif[,,1])

statespace_constrained = build_statespace(
  eastings = eastings, northings = northings, covariates = covariates, 
  # exclude locations whose band 1 value is below average
  linear_constraint = c(-band1_avg, 1, rep(0, nrow(covariates)-2))
)

# identify locations that will be un/defined
valid_locs = which(dat$L7_ETMs.tif[,,1] >= band1_avg, arr.ind = TRUE)
invalid_locs = which(dat$L7_ETMs.tif[,,1] < band1_avg, arr.ind = TRUE)

#
# test: check locations that should exist
#

test_ind = valid_locs[nrow(valid_locs)/2,]

location = extract_statespace_location(
  statespace = statespace_constrained, 
  easting_ind = test_ind['row'] - 1, 
  northing_ind = test_ind['col'] - 1
)

expect_equal(
  location$covariates[-1], # remove intercept
  dat[["L7_ETMs.tif"]][test_ind['row'], test_ind['col'], ]
)

expect_equal(location$easting, eastings[test_ind['row']])
expect_equal(location$northing, northings[test_ind['col']])

#
# test: check locations that should not exist
#

test_ind = invalid_locs[nrow(invalid_locs)/2,]

expect_error(
  extract_statespace_location(
    statespace = statespace_constrained, 
    easting_ind = test_ind['row'] - 1,
    northing_ind = test_ind['col'] - 1
  )
)

#
# test: check that not all states exist near constraint boundaries
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

# investigate all states at this boundary location
states = lapply(c('north', 'east', 'south', 'west'), function(direction) {
  # return a state if defined, or null otherwise
  tryCatch(
    extract_statespace_state(
      statespace = statespace_constrained, 
      last_movement_direction = direction,
      easting_ind = boundary_coord_inds['easting_ind'] - 1, 
      northing_ind = boundary_coord_inds['northing_ind'] - 1
    ),
    error = function(e) { }
  )
})

# verify not all states are defined
expect_gt(sum(sapply(states, is.null)), 0)
