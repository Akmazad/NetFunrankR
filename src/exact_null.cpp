#include <Rcpp.h>
using namespace Rcpp;

// Convolve `acc` in place with kernel `k`, `times` times (direct, exact).
static void conv_times(std::vector<double>& acc, const std::vector<double>& k,
                       int times) {
  int kl = k.size();
  while (kl > 1 && k[kl - 1] == 0.0) kl--;
  int k0 = 0;
  while (k0 < kl && k[k0] == 0.0) k0++;
  for (int t = 0; t < times; t++) {
    int la = acc.size();
    std::vector<double> out(la + kl - 1, 0.0);
    for (int j = k0; j < kl; j++) {
      double w = k[j];
      if (w == 0.0) continue;
      double* o = out.data() + j;
      const double* a = acc.data();
      for (int i = 0; i < la; i++) o[i] += a[i] * w;
    }
    acc.swap(out);
  }
}

// Exact null distributions for one gene set.
//   f      : bins x (dmax+1) matrix, row b = distribution of d_Y over bin b
//   m      : number of query genes per bin
//   mnon   : number of non-member query genes per bin
// Returns list(O = pmf of overlap, S = pmf of total distance,
//              Sb = pmf of distance of non-members under conditional null)
// [[Rcpp::export]]
List exact_null_cpp(NumericMatrix f, IntegerVector m, IntegerVector mnon) {
  int B = f.nrow(), L = f.ncol();
  std::vector<double> pO(1, 1.0), pS(1, 1.0), pSb(1, 1.0);
  for (int b = 0; b < B; b++) {
    if (m[b] == 0) continue;
    std::vector<double> k(L);
    for (int j = 0; j < L; j++) k[j] = f(b, j);
    double pi = k[0];
    std::vector<double> bern(2);
    bern[0] = 1.0 - pi;
    bern[1] = pi;
    conv_times(pO, bern, m[b]);
    conv_times(pS, k, m[b]);
    if (mnon[b] > 0) {
      std::vector<double> g(k);
      g[0] = 0.0;
      double s = 0.0;
      for (int j = 0; j < L; j++) s += g[j];
      if (s > 0) {
        for (int j = 0; j < L; j++) g[j] /= s;
      } else {
        g.assign(std::max(L, 2), 0.0);
        g[1] = 1.0;
      }
      conv_times(pSb, g, mnon[b]);
    }
  }
  return List::create(_["O"] = wrap(pO), _["S"] = wrap(pS), _["Sb"] = wrap(pSb));
}
