#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
// [[Rcpp::export]]
NumericVector f_y_hat_epnechikov(NumericVector Y, NumericVector y, double h) {
  int n = Y.size(), m = y.size();

  // Center Y to protect the moment expansion below from catastrophic
  // cancellation when Y is on a large scale (income, wages, etc.).
  double center = 0.0;
  for (int i = 0; i < n; i++) center += Y[i];
  center /= n;

  NumericVector Ys = clone(Y) - center;
  Ys = Ys.sort();

  double sqrt5  = std::sqrt(5.0);
  double coef   = 3.0 / (4.0 * sqrt5);
  double inv_h  = 1.0 / h;
  double win    = h * sqrt5;
  double inv5h2 = 1.0 / (5.0 * h * h);

  // Prefix sums of the 1st and 2nd moments of Ys -- O(n) to build,
  // O(1) to query the moments of any contiguous range.
  std::vector<double> cs1(n + 1, 0.0), cs2(n + 1, 0.0);
  const double* ys = Ys.begin();
  for (int i = 0; i < n; i++) {
    cs1[i + 1] = cs1[i] + ys[i];
    cs2[i + 1] = cs2[i] + ys[i] * ys[i];
  }

  NumericVector res(m);
  const double* yend = Ys.end();

  for (int i = 0; i < m; i++) {
    double yi = y[i] - center;
    const double* lo_ptr = std::lower_bound(ys, yend, yi - win);
    const double* hi_ptr = std::upper_bound(lo_ptr, yend, yi + win);
    int lo = (int)(lo_ptr - ys);
    int hi = (int)(hi_ptr - ys);

    double S0 = (double)(hi - lo);
    double S1 = cs1[hi] - cs1[lo];
    double S2 = cs2[hi] - cs2[lo];

    // sum_{p in window} (p - yi)^2 = S2 - 2*yi*S1 + yi^2*S0
    double quad = S2 - 2.0 * yi * S1 + yi * yi * S0;
    double s = S0 - quad * inv5h2;   // sum_{p in window} (1 - u^2/5)

    res[i] = s * coef * inv_h / n;
  }
  return res;
}



// [[Rcpp::export]]
IntegerVector rect_counts_rcpp(NumericVector X_sorted, NumericVector x_eval,
                               NumericVector h_vals) {
  const int m = x_eval.size();
  IntegerVector counts(m);
  const double* xs = X_sorted.begin();
  const double* xend = X_sorted.end();

  for (int j = 0; j < m; j++) {
    const double x = x_eval[j];
    const double h = h_vals[j];
    const double* lo = std::lower_bound(xs, xend, x - h);
    const double* hi = std::upper_bound(lo, xend, x + h);
    counts[j] = (int)(hi - lo);
  }
  return counts;
}

// [[Rcpp::export]]
NumericVector counts_to_density(IntegerVector counts, NumericVector h_vals,
                                int n) {
  const int m = counts.size();
  NumericVector f_hat(m);
  const double inv2n = 1.0 / (2.0 * n);
  for (int j = 0; j < m; j++)
    f_hat[j] = counts[j] * inv2n / h_vals[j];
  return f_hat;
}

// [[Rcpp::export]]
NumericVector boot_core(NumericVector Ys, NumericVector Xs,
                        NumericVector Zs, int B) {
  int n1 = Ys.size(), n2 = Xs.size(), n3 = Zs.size();
  NumericVector results(B);
  std::vector<int> counts_y(n1), counts_z(n3), cdf_y_counts(n1);
  std::vector<double> cdf_z(n3);

  GetRNGstate();
  for (int b = 0; b < B; b++) {
    std::fill(counts_y.begin(), counts_y.end(), 0);
    std::fill(counts_z.begin(), counts_z.end(), 0);
    for (int i = 0; i < n1; i++) counts_y[(int)(unif_rand() * n1)]++;
    for (int i = 0; i < n3; i++) counts_z[(int)(unif_rand() * n3)]++;

    int cumul = 0;
    for (int j = 0; j < n1; j++) {
      cumul += counts_y[j];
      cdf_y_counts[j] = cumul;
    }

    cumul = 0;
    for (int j = 0; j < n3; j++) {
      cumul += counts_z[j];
      cdf_z[j] = (double)cumul / n3;
    }

    double s = 0.0;
    for (int i = 0; i < n2; i++) {
      double xb = Xs[(int)(unif_rand() * n2)];
      int pos = (int)(std::upper_bound(Zs.begin(), Zs.end(), xb) - Zs.begin());
      double u = (pos == 0) ? 0.0 : cdf_z[pos - 1];
      int rank_y = (int)std::ceil(u * n1);
      if (rank_y < 1) rank_y = 1;
      if (rank_y > n1) rank_y = n1;
      int idx_y = (int)(std::lower_bound(cdf_y_counts.begin(), cdf_y_counts.end(), rank_y) - cdf_y_counts.begin());
      s += Ys[idx_y];
    }
    results[b] = s / n2;
  }
  PutRNGstate();
  return results;
}