// Image Operations C++ Implementation
// High-performance image processing functions
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <queue>
using namespace Rcpp;
using namespace arma;

//' Connected components labeling (C++ implementation)
//'
//' Fast 8-connected component labeling using union-find.
//'
//' @param binary_image Binary matrix (0/1)
//'
//' @return Integer matrix with component labels
//' @export
// [[Rcpp::export]]
IntegerMatrix connected_components_cpp(IntegerMatrix binary_image) {
  int nrow = binary_image.nrow();
  int ncol = binary_image.ncol();

  IntegerMatrix labels(nrow, ncol);
  std::fill(labels.begin(), labels.end(), 0);

  int current_label = 0;

  // BFS for each unlabeled foreground pixel
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (binary_image(i, j) > 0 && labels(i, j) == 0) {
        current_label++;

        // BFS queue
        std::queue<std::pair<int, int>> q;
        q.push(std::make_pair(i, j));
        labels(i, j) = current_label;

        while (!q.empty()) {
          std::pair<int, int> p = q.front();
          q.pop();
          int pi = p.first;
          int pj = p.second;

          // 8-connectivity neighbors
          for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
              if (di == 0 && dj == 0) continue;

              int ni = pi + di;
              int nj = pj + dj;

              if (ni >= 0 && ni < nrow && nj >= 0 && nj < ncol) {
                if (binary_image(ni, nj) > 0 && labels(ni, nj) == 0) {
                  labels(ni, nj) = current_label;
                  q.push(std::make_pair(ni, nj));
                }
              }
            }
          }
        }
      }
    }
  }

  labels.attr("n_components") = current_label;
  return labels;
}

//' Distance transform (C++ implementation)
//'
//' Compute Euclidean distance transform of binary image.
//'
//' @param binary_image Binary matrix (0 = background, >0 = foreground)
//'
//' @return Numeric matrix with distances to nearest background pixel
//' @export
// [[Rcpp::export]]
NumericMatrix distance_transform_cpp(IntegerMatrix binary_image) {
  int nrow = binary_image.nrow();
  int ncol = binary_image.ncol();

  // Initialize with large values for foreground, 0 for background
  NumericMatrix dist(nrow, ncol);
  double INF = nrow + ncol;

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (binary_image(i, j) > 0) {
        dist(i, j) = INF;
      } else {
        dist(i, j) = 0;
      }
    }
  }

  // Forward pass
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (dist(i, j) > 0) {
        double val = INF;
        if (i > 0) val = std::min(val, dist(i-1, j) + 1.0);
        if (j > 0) val = std::min(val, dist(i, j-1) + 1.0);
        if (i > 0 && j > 0) val = std::min(val, dist(i-1, j-1) + 1.414);
        if (i > 0 && j < ncol-1) val = std::min(val, dist(i-1, j+1) + 1.414);
        dist(i, j) = std::min(dist(i, j), val);
      }
    }
  }

  // Backward pass
  for (int i = nrow - 1; i >= 0; i--) {
    for (int j = ncol - 1; j >= 0; j--) {
      if (dist(i, j) > 0) {
        double val = dist(i, j);
        if (i < nrow-1) val = std::min(val, dist(i+1, j) + 1.0);
        if (j < ncol-1) val = std::min(val, dist(i, j+1) + 1.0);
        if (i < nrow-1 && j < ncol-1) val = std::min(val, dist(i+1, j+1) + 1.414);
        if (i < nrow-1 && j > 0) val = std::min(val, dist(i+1, j-1) + 1.414);
        dist(i, j) = val;
      }
    }
  }

  return dist;
}

//' Gaussian blur (C++ implementation)
//'
//' Fast Gaussian smoothing using separable kernel.
//'
//' @param image Numeric matrix
//' @param sigma Standard deviation of Gaussian kernel
//'
//' @return Smoothed image matrix
//' @export
// [[Rcpp::export]]
NumericMatrix gaussian_blur_cpp(NumericMatrix image, double sigma) {
  int nrow = image.nrow();
  int ncol = image.ncol();

  // Compute kernel size
  int ksize = (int)(6 * sigma);
  if (ksize % 2 == 0) ksize++;
  int half = ksize / 2;

  // Create 1D Gaussian kernel
  NumericVector kernel(ksize);
  double sum = 0;
  for (int i = 0; i < ksize; i++) {
    double x = i - half;
    kernel[i] = exp(-x * x / (2 * sigma * sigma));
    sum += kernel[i];
  }
  for (int i = 0; i < ksize; i++) {
    kernel[i] /= sum;
  }

  // Horizontal pass
  NumericMatrix temp(nrow, ncol);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double val = 0;
      double wsum = 0;
      for (int k = -half; k <= half; k++) {
        int jj = j + k;
        if (jj >= 0 && jj < ncol) {
          val += image(i, jj) * kernel[k + half];
          wsum += kernel[k + half];
        }
      }
      temp(i, j) = val / wsum;
    }
  }

  // Vertical pass
  NumericMatrix result(nrow, ncol);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double val = 0;
      double wsum = 0;
      for (int k = -half; k <= half; k++) {
        int ii = i + k;
        if (ii >= 0 && ii < nrow) {
          val += temp(ii, j) * kernel[k + half];
          wsum += kernel[k + half];
        }
      }
      result(i, j) = val / wsum;
    }
  }

  return result;
}

//' Extract trace from ROI mask (C++ implementation)
//'
//' Fast trace extraction using mask.
//'
//' @param movie 3D array (height x width x frames)
//' @param mask Binary mask matrix
//' @param method "mean", "sum", or "weighted"
//'
//' @return Numeric vector of trace values
//' @export
// [[Rcpp::export]]
NumericVector extract_trace_cpp(NumericVector movie, NumericMatrix mask, String method = "mean") {
  IntegerVector dims = movie.attr("dim");
  int nrow = dims[0];
  int ncol = dims[1];
  int nframes = dims[2];

  // Find mask pixels
  std::vector<int> mask_i, mask_j;
  std::vector<double> mask_w;
  double wsum = 0;

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (mask(i, j) > 0) {
        mask_i.push_back(i);
        mask_j.push_back(j);
        mask_w.push_back(mask(i, j));
        wsum += mask(i, j);
      }
    }
  }

  int npix = mask_i.size();
  NumericVector trace(nframes);

  for (int f = 0; f < nframes; f++) {
    double val = 0;

    for (int p = 0; p < npix; p++) {
      int i = mask_i[p];
      int j = mask_j[p];
      double w = mask_w[p];

      double pix_val = movie[i + j * nrow + f * nrow * ncol];

      if (method == "weighted") {
        val += pix_val * w;
      } else {
        val += pix_val;
      }
    }

    if (method == "mean") {
      trace[f] = val / npix;
    } else if (method == "weighted") {
      trace[f] = val / wsum;
    } else {
      trace[f] = val;
    }
  }

  return trace;
}

//' Local correlation image (C++ implementation)
//'
//' Compute pixel-wise correlation with neighbors for cell detection.
//'
//' @param movie 3D array (height x width x frames)
//' @param radius Neighborhood radius
//'
//' @return Correlation image matrix
//' @export
// [[Rcpp::export]]
NumericMatrix local_correlation_cpp(NumericVector movie, int radius = 1) {
  IntegerVector dims = movie.attr("dim");
  int nrow = dims[0];
  int ncol = dims[1];
  int nframes = dims[2];

  NumericMatrix corr_img(nrow, ncol);

  // Pre-compute mean and variance for each pixel
  mat means(nrow, ncol, fill::zeros);
  mat vars(nrow, ncol, fill::zeros);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double sum = 0, sum2 = 0;
      for (int f = 0; f < nframes; f++) {
        double v = movie[i + j * nrow + f * nrow * ncol];
        sum += v;
        sum2 += v * v;
      }
      means(i, j) = sum / nframes;
      vars(i, j) = sqrt(sum2 / nframes - means(i, j) * means(i, j) + 1e-10);
    }
  }

  // Compute local correlations
  for (int i = radius; i < nrow - radius; i++) {
    for (int j = radius; j < ncol - radius; j++) {
      double corr_sum = 0;
      int count = 0;

      for (int di = -radius; di <= radius; di++) {
        for (int dj = -radius; dj <= radius; dj++) {
          if (di == 0 && dj == 0) continue;

          int ni = i + di;
          int nj = j + dj;

          // Compute correlation between pixels
          double cov = 0;
          for (int f = 0; f < nframes; f++) {
            double v1 = movie[i + j * nrow + f * nrow * ncol] - means(i, j);
            double v2 = movie[ni + nj * nrow + f * nrow * ncol] - means(ni, nj);
            cov += v1 * v2;
          }
          cov /= nframes;

          double corr = cov / (vars(i, j) * vars(ni, nj) + 1e-10);
          corr_sum += corr;
          count++;
        }
      }

      corr_img(i, j) = corr_sum / count;
    }
  }

  return corr_img;
}

//' Median filter (C++ implementation)
//'
//' Fast median filtering for image denoising.
//'
//' @param image Numeric matrix
//' @param radius Filter radius
//'
//' @return Filtered image matrix
//' @export
// [[Rcpp::export]]
NumericMatrix median_filter_cpp(NumericMatrix image, int radius) {
  int nrow = image.nrow();
  int ncol = image.ncol();

  NumericMatrix result(nrow, ncol);

  int window_size = (2 * radius + 1) * (2 * radius + 1);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      std::vector<double> vals;
      vals.reserve(window_size);

      for (int di = -radius; di <= radius; di++) {
        for (int dj = -radius; dj <= radius; dj++) {
          int ni = i + di;
          int nj = j + dj;
          if (ni >= 0 && ni < nrow && nj >= 0 && nj < ncol) {
            vals.push_back(image(ni, nj));
          }
        }
      }

      // Find median
      int n = vals.size();
      std::nth_element(vals.begin(), vals.begin() + n/2, vals.end());
      result(i, j) = vals[n/2];
    }
  }

  return result;
}

//' Threshold image (C++ implementation)
//'
//' @param image Numeric matrix
//' @param threshold Threshold value
//' @param method "binary", "otsu", or "adaptive"
//'
//' @return Binary integer matrix
//' @export
// [[Rcpp::export]]
IntegerMatrix threshold_image_cpp(NumericMatrix image, double threshold, String method = "binary") {
  int nrow = image.nrow();
  int ncol = image.ncol();

  IntegerMatrix result(nrow, ncol);

  if (method == "binary") {
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        result(i, j) = (image(i, j) > threshold) ? 1 : 0;
      }
    }
  } else if (method == "otsu") {
    // Compute histogram
    int nbins = 256;
    NumericVector hist(nbins, 0.0);
    double minval = image(0, 0), maxval = image(0, 0);

    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        if (image(i, j) < minval) minval = image(i, j);
        if (image(i, j) > maxval) maxval = image(i, j);
      }
    }

    double range = maxval - minval + 1e-10;
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        int bin = (int)((image(i, j) - minval) / range * (nbins - 1));
        hist[bin] += 1.0;
      }
    }

    // Normalize histogram
    double total = nrow * ncol;
    for (int i = 0; i < nbins; i++) {
      hist[i] /= total;
    }

    // Find optimal threshold (Otsu's method)
    double best_thresh = 0, best_var = 0;
    double w0 = 0, mu0 = 0, mu_total = 0;

    for (int i = 0; i < nbins; i++) {
      mu_total += i * hist[i];
    }

    for (int t = 0; t < nbins; t++) {
      w0 += hist[t];
      if (w0 == 0) continue;

      double w1 = 1.0 - w0;
      if (w1 == 0) break;

      mu0 += t * hist[t];
      double mu1 = (mu_total - mu0) / w1;
      double mu0_norm = mu0 / w0;

      double var = w0 * w1 * (mu0_norm - mu1) * (mu0_norm - mu1);
      if (var > best_var) {
        best_var = var;
        best_thresh = minval + t * range / nbins;
      }
    }

    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        result(i, j) = (image(i, j) > best_thresh) ? 1 : 0;
      }
    }

    result.attr("threshold") = best_thresh;
  }

  return result;
}
