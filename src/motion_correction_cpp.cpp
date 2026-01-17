// Motion Correction C++ Implementation
// High-performance motion correction using Rcpp
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Shift image by given displacement (C++ implementation)
//'
//' Fast image shifting using bilinear interpolation.
//'
//' @param image Numeric matrix representing the image
//' @param shift_x X displacement (can be fractional)
//' @param shift_y Y displacement (can be fractional)
//'
//' @return Shifted image matrix
//' @export
// [[Rcpp::export]]
NumericMatrix shift_image_cpp(NumericMatrix image, double shift_x, double shift_y) {
  int nrow = image.nrow();
  int ncol = image.ncol();

  NumericMatrix result(nrow, ncol);

  // Integer and fractional parts
  int ix = (int)floor(shift_x);
  int iy = (int)floor(shift_y);
  double fx = shift_x - ix;
  double fy = shift_y - iy;

  // Weights for bilinear interpolation
  double w00 = (1.0 - fx) * (1.0 - fy);
  double w01 = (1.0 - fx) * fy;
  double w10 = fx * (1.0 - fy);
  double w11 = fx * fy;

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      // Source coordinates
      int si = i - iy;
      int sj = j - ix;

      // Check bounds and interpolate
      double val = 0.0;

      if (si >= 0 && si < nrow && sj >= 0 && sj < ncol) {
        val += w00 * image(si, sj);
      }
      if (si >= 0 && si < nrow && sj + 1 >= 0 && sj + 1 < ncol) {
        val += w10 * image(si, sj + 1);
      }
      if (si + 1 >= 0 && si + 1 < nrow && sj >= 0 && sj < ncol) {
        val += w01 * image(si + 1, sj);
      }
      if (si + 1 >= 0 && si + 1 < nrow && sj + 1 >= 0 && sj + 1 < ncol) {
        val += w11 * image(si + 1, sj + 1);
      }

      result(i, j) = val;
    }
  }

  return result;
}

//' Compute phase correlation between two images
//'
//' Fast phase correlation for motion estimation.
//'
//' @param ref Reference image matrix
//' @param target Target image matrix
//' @param max_shift Maximum allowed shift
//'
//' @return NumericVector with c(shift_x, shift_y, correlation)
//' @export
// [[Rcpp::export]]
NumericVector phase_correlation_cpp(NumericMatrix ref, NumericMatrix target, int max_shift) {
  int nrow = ref.nrow();
  int ncol = ref.ncol();

  // Convert to arma matrices for FFT
  mat ref_mat(ref.begin(), nrow, ncol, false);
  mat target_mat(target.begin(), nrow, ncol, false);

  // Remove mean
  ref_mat -= mean(mean(ref_mat));
  target_mat -= mean(mean(target_mat));

  // Compute cross-correlation in spatial domain (simplified)
  // For a full implementation, would use FFT via fftw3
  double best_corr = -1e9;
  int best_dx = 0;
  int best_dy = 0;

  // Search over shift range
  for (int dy = -max_shift; dy <= max_shift; dy++) {
    for (int dx = -max_shift; dx <= max_shift; dx++) {
      double corr = 0.0;
      double n = 0.0;

      // Compute normalized cross-correlation
      for (int i = max_shift; i < nrow - max_shift; i++) {
        for (int j = max_shift; j < ncol - max_shift; j++) {
          int ti = i + dy;
          int tj = j + dx;

          if (ti >= 0 && ti < nrow && tj >= 0 && tj < ncol) {
            corr += ref_mat(i, j) * target_mat(ti, tj);
            n += 1.0;
          }
        }
      }

      if (n > 0) {
        corr /= n;
        if (corr > best_corr) {
          best_corr = corr;
          best_dx = dx;
          best_dy = dy;
        }
      }
    }
  }

  NumericVector result(3);
  result[0] = (double)best_dx;
  result[1] = (double)best_dy;
  result[2] = best_corr;

  return result;
}

//' Motion correct entire movie
//'
//' Apply motion correction to all frames efficiently.
//'
//' @param movie 3D array (height x width x frames)
//' @param template Reference image
//' @param max_shift Maximum shift
//'
//' @return List with corrected movie and shifts
//' @export
// [[Rcpp::export]]
List motion_correct_movie_cpp(NumericVector movie, NumericMatrix tmpl, int max_shift) {
  IntegerVector dims = movie.attr("dim");
  int nrow = dims[0];
  int ncol = dims[1];
  int nframes = dims[2];

  // Output arrays
  NumericVector corrected(nrow * ncol * nframes);
  corrected.attr("dim") = dims;

  NumericMatrix shifts(nframes, 2);

  // Process each frame
  for (int f = 0; f < nframes; f++) {
    // Extract frame
    NumericMatrix frame(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        frame(i, j) = movie[i + j * nrow + f * nrow * ncol];
      }
    }

    // Compute shift
    NumericVector shift_result = phase_correlation_cpp(tmpl, frame, max_shift);
    shifts(f, 0) = shift_result[0];
    shifts(f, 1) = shift_result[1];

    // Apply shift
    NumericMatrix corrected_frame = shift_image_cpp(frame, shift_result[0], shift_result[1]);

    // Store result
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        corrected[i + j * nrow + f * nrow * ncol] = corrected_frame(i, j);
      }
    }
  }

  return List::create(
    Named("corrected") = corrected,
    Named("shifts") = shifts
  );
}

//' Compute template by averaging frames
//'
//' @param movie 3D array
//' @param n_frames Number of frames to average (0 = all)
//'
//' @return Template image matrix
//' @export
// [[Rcpp::export]]
NumericMatrix compute_template_cpp(NumericVector movie, int n_frames = 0) {
  IntegerVector dims = movie.attr("dim");
  int nrow = dims[0];
  int ncol = dims[1];
  int total_frames = dims[2];

  if (n_frames <= 0 || n_frames > total_frames) {
    n_frames = total_frames;
  }

  NumericMatrix tmpl(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double sum = 0.0;
      for (int f = 0; f < n_frames; f++) {
        sum += movie[i + j * nrow + f * nrow * ncol];
      }
      tmpl(i, j) = sum / n_frames;
    }
  }

  return tmpl;
}
