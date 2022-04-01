#include "grid.h"
#include "gtest/gtest.h"

namespace {

using grid_h::grid;
using grid_h::eigenvectors;

typedef double mytype;
// typedef cplx mytype;
typedef grid<mytype> mygrid;
typedef vector<mytype> myvector;

// Compute eigenvectors/eigenvalues of the difficult Rosser Matrix
TEST(Rosser, All) {
  grid<double> rosser(8, 8,
      { 611, 196, -192, 407, -8, -52, -49, 29,
        196, 899, 113, -192, -71, -43, -8, -44,
        -192, 113, 899, 196, 61, 49, 8, 52,
        407, -192, 196, 611, 8, 44, 59, -23,
        -8, -71, 61, 8, 411, -599, 208, 208,
        -52, -43, 49, 44, -599, 411, 208, 208,
        -49, -8, 8, 59, 208, 208, 99, -911,
        29, -44, 52, -23, 208, 208, -911, 99 };
  grid<double> evecs, evals;
  info = eigenvectors(rosser, evecs, evals);
  EXPECT_EQ(info, 0);
  EXPECT_EQ(evecs, grid<double>());
  EXPECT_EQ(evals, grid<double>());
  EXPECT_EQ(linfvecnorm(
      matmult(rosser, evecs) - diagmult(evecs, evals.storage())), 0);
}

};  // namespace
