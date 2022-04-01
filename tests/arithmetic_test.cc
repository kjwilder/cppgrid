#include <iostream>
#include <vector>
#include "grid.h"
#include "gtest/gtest.h"

namespace {

using grid_h::grid;

TEST(Arithmetic, Basics) {
  grid<double> x(2, 3, {1, 2, 3, 4, 5, 6});
  grid<double> y(2, 3, {-1, 3, -5, 7, -9, 11});
  EXPECT_EQ(x + y,
            grid<double>(2, 3, {0, 5, -2, 11, -4, 17}));
  EXPECT_EQ(x - y,
            grid<double>(2, 3, {2, -1, 8, -3, 14, -5}));
  EXPECT_EQ(x * y,
            grid<double>(2, 3, {-1, 6, -15, 28, -45, 66}));
  EXPECT_EQ(x / y,
            grid<double>(2, 3, {-1, 2.0 / 3.0, -0.6, 4.0 / 7.0,
                                -5.0 / 9.0, 6.0 / 11.0}));
}

TEST(Arithmetic, Equations) {
  grid<double> x(2, 3, {1, 2, 3, 4, 5, 6});
  grid<double> y(2, 3, {-1, 3, -5, 7, -9, 11});
  EXPECT_EQ((x + 4) * y,
            grid<double>(2, 3, {-5, 18, -35, 56, -81, 110}));
  EXPECT_EQ(y / (x + 4),
            grid<double>(2, 3, {-0.2, 0.5, -5.0 / 7.0, 0.875, -1, 1.1}));
  EXPECT_EQ(x * x - y * y,
            grid<double>(2, 3, {0, -5, -16, -33, -56, -85}));
  EXPECT_EQ((x - y) * (x + y), x * x - y * y);
}

TEST(Arithmetic, Binding) {
  grid<double> x(2, 3, {1, 2, 3, 4, 5, 6});
  grid<double> y(2, 3, {-1, 3, -5, 7, -9, 11});
  auto z = x;
  z.cbind(y);
  EXPECT_EQ(z, grid<double>(2, 6, {1, 2, 3, 4, 5, 6, -1, 3, -5, 7, -9, 11}));
  z = x;
  z.rbind(y);
  EXPECT_EQ(z, grid<double>(4, 3, {1, 2, -1, 3, 3, 4, -5, 7, 5, 6, -9, 11}));
}

TEST(Arithmetic, Oddities) {
  grid<int> x(2, 3, {1, 2, 3, 4, 5, 6});
  grid<int> y(2, 3, {1, 3, 5, 7, 9, 11});
  grid<int> z(2, 3, {2, 0, -5, 3, -1, 4});
  grid<int> w = z;
  w += w += x;
  EXPECT_EQ(w, grid<int>(2, 3, {6, 4, -4, 14, 8, 20}));
  w = z;
  w -= w -= x;
  EXPECT_EQ(w, grid<int>(2, 3));
}

};  // namespace

