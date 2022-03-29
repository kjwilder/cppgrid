#include <chrono>
#include <iostream>
#include <thread>

#include "grid.h"
#include "gtest/gtest.h"

using grid_h::grid;

namespace {

TEST(Grid, ConstructorEmpty) {
  grid<int> g;
  EXPECT_EQ(g.rows(), 0);
  EXPECT_EQ(g.cols(), 0);
}

TEST(Grid, ConstructorRowsOnly) {
  grid<int> g(50);
  EXPECT_EQ(g.rows(), 50);
  EXPECT_EQ(g.cols(), 1);
}

TEST(Grid, ConstructorRowsAndColumns) {
  grid<int> g(10, 20);
  EXPECT_EQ(g.rows(), 10);
  EXPECT_EQ(g.cols(), 20);
}

TEST(Grid, ConstructorZeroColumns) {
  grid<int> g(10, 0);
  EXPECT_EQ(g.rows(), 10);
  EXPECT_EQ(g.cols(), 0);
  EXPECT_EQ(g.storage().size(), 0);
}

TEST(Grid, OrderingIsColumnMajor) {
  grid<int> g(2, 2, {1, 2, 3, 4});
  EXPECT_EQ(g(1, 0), 2);
  EXPECT_EQ(g(0, 1), 3);
}

TEST(Grid, Inverse) {
  grid<double> gi = grid<double>(2, 2, {2, 2, 1, 4}).inverse();
  EXPECT_NEAR(gi(0, 0), 2.0 / 3.0, 1e-6);
  EXPECT_NEAR(gi(1, 0), -1.0 / 3.0, 1e-6);
  EXPECT_NEAR(gi(0, 1), -1.0 / 6.0, 1e-6);
  EXPECT_NEAR(gi(1, 1), 1.0 / 3.0, 1e-6);
}

TEST(Grid, ScaleInt) {
  EXPECT_EQ(grid<int>({3, 5, 8}).scale(3),
            grid<int>({1, 1, 3}));
  EXPECT_EQ(grid<int>({3, 5, -8}).scale(3),
            grid<int>({1, 1, -3}));

  EXPECT_EQ(grid<int>().scale(5), grid<int>());
  EXPECT_EQ(grid<int>(1, 0).scale(5), grid<int>(1, 0));

  EXPECT_EQ(grid<int>(2, 2, {0, 0, 0, 0}).scale(3),
            grid<int>(2, 2, {0, 0, 0, 0}));
}

TEST(Grid, ScaleDouble) {
  EXPECT_EQ(grid<double>({3, 5, 8}).scale(3),
            grid<double>({9.0 / 8.0, 15.0 / 8.0, 3}));
}

TEST(Grid, TransformInt) {
  EXPECT_EQ(grid<int>({1, 2, 4, 8}).transform(5, 20),
            grid<int>({5, 7, 11, 20}));
}

TEST(Grid, TransformDouble) {
  auto g = grid<double>({1, 2, 4, 8}).transform(5, 20);
  EXPECT_NEAR(g(0), 5.0, 1e-6);
  EXPECT_NEAR(g(1), 50.0 / 7.0, 1e-6);
  EXPECT_NEAR(g(2), 80.0 / 7.0, 1e-6);
  EXPECT_NEAR(g(3), 20.0, 1e-6);
}

TEST(Grid, TransformConstants) {
  // Transform [1, 1, 1, 1] -> range [5, 20], transform to min value 5.
  EXPECT_EQ(grid<double>(4).fill(1).transform(5, 20), 5);
  // Transform [10, 10, 10, 10] -> range [5, 20], in-range value 10 unchanged.
  EXPECT_EQ(grid<double>(4).fill(10).transform(5, 20), 10);
  // Transform [25, 25, 25, 25] -> range [5, 20], transform to max value 20.
  EXPECT_EQ(grid<double>(4).fill(25).transform(5, 20), 20);
}

TEST(Grid, OperatorPlusEqual) {
  auto gd = grid<double>({1, 2, 3, 4});
  EXPECT_EQ(gd += grid<double>({2, 5, 7, 12}),
            grid<double>({3, 7, 10, 16}));
  EXPECT_EQ(gd, grid<double>({3, 7, 10, 16}));

  gd = grid<double>({1, 2, 3, 4});
  EXPECT_EQ(gd += 7, grid<double>({8, 9, 10, 11}));
  EXPECT_EQ(gd, grid<double>({8, 9, 10, 11}));
}

TEST(Grid, OperatorPlus) {
  auto gd = grid<double>({1, 2, 3, 4});
  EXPECT_EQ(gd + grid<double>({2, 5, 7, 12}),
            grid<double>({3, 7, 10, 16}));
  EXPECT_EQ(gd, grid<double>({1, 2, 3, 4}));

  gd = grid<double>({1, 2, 3, 4});
  EXPECT_EQ(gd + 7, grid<double>({8, 9, 10, 11}));
  EXPECT_EQ(gd, grid<double>({1, 2, 3, 4}));
}

TEST(Grid, OperatorMinusEqual) {
  EXPECT_EQ(grid<double>({1, 2, 3, 4}) -= grid<double>({2, 5, 7, 12}),
            grid<double>({-1, -3, -4, -8}));
}

}  // namespace
