source(file.path('..','..','R','pipeline.R'))

context('infer_spikes')

test_that('infer_spikes returns fit and spike columns', {
  trace <- rep(0, 10)
  result <- suppressWarnings(infer_spikes(trace, penalty = 1, ar_order = 1))
  expect_true(is.data.frame(result))
  expect_equal(names(result), c('fit', 'spike'))
  expect_equal(nrow(result), length(trace))
})
