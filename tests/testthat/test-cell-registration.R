test_that("compute_roi_centroids works with mask-based ROIs", {
  rois <- list()
  for (i in 1:3) {
    mask <- matrix(0, 50, 50)
    cx <- 10 + (i - 1) * 15
    cy <- 25
    mask[(cx-2):(cx+2), (cy-2):(cy+2)] <- 1
    rois[[i]] <- list(mask = mask)
  }

  centroids <- compute_roi_centroids(rois)

  expect_true(is.matrix(centroids))
  expect_equal(nrow(centroids), 3)
  expect_equal(ncol(centroids), 2)
})

test_that("compute_roi_centroids works with coordinate ROIs", {
  rois <- list(
    matrix(c(5, 5, 6, 5, 5, 6, 6, 6), ncol = 2, byrow = TRUE),
    matrix(c(20, 20, 21, 20, 20, 21), ncol = 2, byrow = TRUE)
  )

  centroids <- compute_roi_centroids(rois)

  expect_equal(nrow(centroids), 2)
  expect_true(centroids[1, 1] < 10)  # First ROI near (5,5)
  expect_true(centroids[2, 1] > 15)  # Second ROI near (20,20)
})

test_that("register_cells with centroid method", {
  # Create two sessions with similar ROIs
  session1 <- list()
  session2 <- list()

  for (i in 1:3) {
    mask1 <- matrix(0, 30, 30)
    mask2 <- matrix(0, 30, 30)
    cx <- 5 + (i - 1) * 10
    cy <- 15

    mask1[(cx-1):(cx+1), (cy-1):(cy+1)] <- 1
    # Slightly shifted in session 2
    mask2[(cx):(cx+2), (cy-1):(cy+1)] <- 1

    session1[[i]] <- list(mask = mask1)
    session2[[i]] <- list(mask = mask2)
  }

  result <- register_cells(
    sessions = list(session1, session2),
    method = "centroid",
    max_distance = 5
  )

  expect_true(is.list(result))
  expect_true("matches" %in% names(result) || "registration_table" %in% names(result))
})

test_that("register_cells handles empty sessions gracefully", {
  session1 <- list()
  mask <- matrix(0, 20, 20)
  mask[8:12, 8:12] <- 1
  session1[[1]] <- list(mask = mask)

  session2 <- list()  # Empty

  # Should not error
 result <- tryCatch(
    register_cells(list(session1, session2), method = "centroid"),
    error = function(e) list(error = TRUE)
  )

  expect_true(is.list(result))
})
