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

# get grid indices for a valid location
test_ind = c(
  easting_ind = unname(valid_locs[2e4,'row']), 
  northing_ind = unname(valid_locs[2e4,'col'])
)

set.seed(2023)

# sample points on grid
states = sample_gaussian_states(
  statespace_search = search, 
  easting = eastings[test_ind['easting_ind']], 
  northing = northings[test_ind['northing_ind']], 
  semi_major = 1e3, 
  semi_minor = 1e2, 
  orientation = 45, 
  n = 1e4
)

# extract sampled locations
pts = do.call(rbind, lapply(states$states, function(s) {
  unlist(s$location[c('easting', 'northing')])
}))

# extract last movement directions
last_movement_directions = sapply(states$states, function(s) {
  s$last_movement_direction
})

# basic tests: we get many coordinates and movement directions
expect_gt(length(unique(pts[,1])), 1)
expect_gt(length(unique(pts[,2])), 1)
expect_equal(length(unique(last_movement_directions)), 4)

# #
# # explore sampled movement directions
# #
# 
# last_movement_directions = sapply(states$states, function(s) {
#   s$last_movement_direction
# })
# 
# table(last_movement_directions) / length(last_movement_directions)

# #
# # map sampled points
# #
# 
# dat2 = dat
# dat2$L7_ETMs.tif[,,1] = dat$L7_ETMs.tif[,,1] >= band1_avg
# 
# image(dat2)
# 
# points(pts[,1], pts[,2], pch = '.', col = 'green')
# 
# points(
#   eastings[test_ind['easting_ind']],
#   northings[test_ind['northing_ind']],
#   col = 'blue'
# )
