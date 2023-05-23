#
# test: directional covariates are implemented correctly
#

expect_equal(TestDirectionalCovariate('north', 'north'), 1)
expect_equal(TestDirectionalCovariate('north', 'east'), 0)
expect_equal(TestDirectionalCovariate('north', 'south'), -1)
expect_equal(TestDirectionalCovariate('north', 'west'), 0)

expect_equal(TestDirectionalCovariate('east', 'north'), 0)
expect_equal(TestDirectionalCovariate('east', 'east'), 1)
expect_equal(TestDirectionalCovariate('east', 'south'), 0)
expect_equal(TestDirectionalCovariate('east', 'west'), -1)

expect_equal(TestDirectionalCovariate('south', 'north'), -1)
expect_equal(TestDirectionalCovariate('south', 'east'), 0)
expect_equal(TestDirectionalCovariate('south', 'south'), 1)
expect_equal(TestDirectionalCovariate('south', 'west'), 0)

expect_equal(TestDirectionalCovariate('west', 'north'), 0)
expect_equal(TestDirectionalCovariate('west', 'east'), -1)
expect_equal(TestDirectionalCovariate('west', 'south'), 0)
expect_equal(TestDirectionalCovariate('west', 'west'), 1)
