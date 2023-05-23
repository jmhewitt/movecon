#
# test: directional covariates are implemented correctly
#

expect_equal(Test__Directional_Covariate('north', 'north'), 1)
expect_equal(Test__Directional_Covariate('north', 'east'), 0)
expect_equal(Test__Directional_Covariate('north', 'south'), -1)
expect_equal(Test__Directional_Covariate('north', 'west'), 0)

expect_equal(Test__Directional_Covariate('east', 'north'), 0)
expect_equal(Test__Directional_Covariate('east', 'east'), 1)
expect_equal(Test__Directional_Covariate('east', 'south'), 0)
expect_equal(Test__Directional_Covariate('east', 'west'), -1)

expect_equal(Test__Directional_Covariate('south', 'north'), -1)
expect_equal(Test__Directional_Covariate('south', 'east'), 0)
expect_equal(Test__Directional_Covariate('south', 'south'), 1)
expect_equal(Test__Directional_Covariate('south', 'west'), 0)

expect_equal(Test__Directional_Covariate('west', 'north'), 0)
expect_equal(Test__Directional_Covariate('west', 'east'), -1)
expect_equal(Test__Directional_Covariate('west', 'south'), 0)
expect_equal(Test__Directional_Covariate('west', 'west'), 1)
