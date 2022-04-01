#include <cstdlib>
#include <iostream>
#include "grid.h"
#include "gtest/gtest.h"

namespace {

using grid_h::grid;

/*
#ifdef DGRID
typedef dgrid mygrid;
typedef dvector myvector;
#else
typedef zgrid mygrid;
typedef zvector myvector;
#endif
*/

TEST(Statistics, Mean) {
  EXPECT_EQ(grid<double>("hilbert", 3),
            grid<double>(3, 3, {1.0, 1.0 / 2.0, 1.0 / 3.0,
                                1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0,
                                1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0}));
}

/*
  h.cbind(h);  // Test cbinding a grid to itself.
  dgrid a;
  a = h;  // Test copying
  cout << "A =\n" << a << endl;

  dgrid b0("hilbert", 3);
  b0.rbind(b0);  // Test rbinding a grid to itself.
  a = b0;  // Test copying
  cout << "A =\n" << a << endl;

  myvector sumrows(a.rows()), sumcols(a.cols());
  cout << "A row sums\n";
  for (uint i = 0; i < a.rows(); i++)
    sumrows[i] = sum(a, "row", i);
  cout << sumrows;
  cout << "A column sums\n";
  for (uint i = 0; i < a.cols(); i++)
    sumcols[i] = sum(a, "col", i);
  cout << sumcols;

  myvector means(a.rows());
  cout << "A row means\n";
  for (uint i = 0; i < a.rows(); i++)
    means[i] = mean(a, "row", i);
  cout << means;
  means.resize(a.cols());
  cout << "A column means\n";
  for (uint i = 0; i < a.cols(); i++)
    means[i] = mean(a, "col", i);
  cout << means;

  dgrid cv;
  cout << "\nCol covariances =\n" << cov(a, cv) << endl;
  cout << "Row covariances =\n" << cov(a, cv, "row") << endl;
  cor(a, cv);
  cout << "Col correlations =\n" << cv << endl;
  cor(a, cv, "row");
  cout << "Row correlations =\n" << cv << endl;
*/

};  // namespace
