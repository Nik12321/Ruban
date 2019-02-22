test_that("nuclear function", {
  expect_error(nuclearFunction(10, r = c (10, 13)))
  expect_error(nuclearFunction(10, r = 3.4))
  expect_error(nuclearFunction(10, r = 0))
  expect_error(nuclearFunction(10, r = -1))
  expect_error(nuclearFunction(10, r = "4"))
  expect_error(nuclearFunction(10, r = NULL))
  expect_error(nuclearFunction(10, r = NaN))

  expect_error(nuclearFunction(10, s = c (123, 123)))
  expect_error(nuclearFunction(10, s = -1))
  expect_error(nuclearFunction(10, s = "123"))
  expect_error(nuclearFunction(10, s = NaN))
  expect_error(nuclearFunction(10, s = NULL))

  expect_error(nuclearFunction(z = c (12, 123)))
  expect_error(nuclearFunction(z = NULL))
  expect_error(nuclearFunction(z = "123"))
  expect_error(nuclearFunction(z = NaN))

  expect_warning(nuclearFunction(10, x = -1))
  expect_warning(nuclearFunction(10, x = NaN))
  expect_warning(nuclearFunction(10, x = c (10, 10)))
  expect_warning(nuclearFunction(10, x = NULL))

})

test_that("Test Rouban algorithm for uncorrect values", {
  expect_error(Rouban(x = NaN, deltaX = c(10, 10), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(123, NaN), deltaX = c(10, 10), f = function(x) {
    x[1] + x[2]
    }
                     , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = NULL, deltaX = c(10, 10), f = function(x) {
    x[1] + x[2]
    }
                     , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(123, NULL), deltaX = c(10, 10), f = function(x) {
    x[1] + x[2]
    }
                     , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = 123, deltaX = c(10, 10), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(123, 100, 100), deltaX = c(10, 10), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))

  expect_error(Rouban(x = c(10, 10), deltaX = NaN, f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(123, NaN), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = NULL, f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(123, NULL), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = 123, f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(123, 123, 123), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10)))

  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = NaN, upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, NaN), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = NULL, upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, NULL), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = 123, upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10, -10), upper = c(10, 10)))

  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , upper = NaN, lower = c(-10, -10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , upper = c(123, NaN), lower = c(-10, -10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , upper = NULL, lower = c(-10, -10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , upper = c(123, NULL), lower = c(-10, -10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , upper = 123, lower = c(-10, -10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , upper = c(123, 123, 1230), lower = c(-10, -10)))


  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), n = 0))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), n = NULL))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), n = 23.3))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), n = NaN))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), n = c(1, 2)))

  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), y = -0.1))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), y = NULL))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), y = NaN))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), y = c(1, 2)))

  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), q = 0))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), q = NULL))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), q = 1.2))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), q = NaN))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), q = c(1, 2)))

  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), e = 0))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), e = NULL))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), e = NaN))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, -10), upper = c(10, 10), e = c(1, 2)))

  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(20, -10), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(-10, 20), upper = c(10, 10)))
  expect_error(Rouban(x = c(10, 10), deltaX = c(20, 20), f = function(x) {
    x[1] + x[2]
    }
    , lower = c(20, 20), upper = c(10, 10)))
})
