#include <chrono>
#include <iostream>
#include <thread>

#include "grid.h"
#include "gtest/gtest.h"

namespace {

TEST(Grid, ConstructEmpty) {
  ucgrid g;
  EXPECT_EQ(g.rows(), 0);
  EXPECT_EQ(g.cols(), 0);
}

TEST(Grid, ConstructRowsOnly) {
  ucgrid g(50);
  EXPECT_EQ(g.rows(), 50);
  EXPECT_EQ(g.cols(), 1);
}

TEST(Grid, ConstructRowsAndColumns) {
  ucgrid g(10, 20);
  EXPECT_EQ(g.rows(), 10);
  EXPECT_EQ(g.cols(), 20);
}

TEST(Grid, Inverse) {
  dgrid g(2, 2);
  g(0, 0) = 2;
  g(0, 1) = 1;
  g(1, 0) = 2;
  g(1, 1) = 4;
  dgrid gi = g.inverse();
  EXPECT_NEAR(gi(0, 0), 2.0 / 3.0, 1e-6);
  EXPECT_NEAR(gi(0, 1), -1.0 / 6.0, 1e-6);
  EXPECT_NEAR(gi(1, 0), -1.0 / 3.0, 1e-6);
  EXPECT_NEAR(gi(1, 1), 1.0 / 3.0, 1e-6);
}

TEST(Grid, TransformInt) {
  igrid g(4);
  g(0) = 1;
  g(1) = 2;
  g(2) = 4;
  g(3) = 8;
  g.transform(5, 20);
  EXPECT_EQ(g(0), 5);
  EXPECT_EQ(g(1), 7);
  EXPECT_EQ(g(2), 11);
  EXPECT_EQ(g(3), 20);
}

TEST(Grid, TransformDouble) {
  dgrid g(4);
  g(0) = 1;
  g(1) = 2;
  g(2) = 4;
  g(3) = 8;
  g.transform(5, 20);
  EXPECT_NEAR(g(0), 5.0, 1e-6);
  EXPECT_NEAR(g(1), 7.1428571, 1e-6);
  EXPECT_NEAR(g(2), 11.4285714, 1e-6);
  EXPECT_NEAR(g(3), 20.0, 1e-6);
}

TEST(Grid, TransformLowConstant) {
  dgrid g(4);
  g.clear(1);
  g.transform(5, 20);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(g(i), 5);
  }
}

TEST(Grid, TransformBetweenConstant) {
  dgrid g(4);
  g.clear(10);
  g.transform(5, 20);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(g(i), 10);
  }
}

TEST(Grid, TransformHighConstant) {
  dgrid g(4);
  g.clear(25);
  g.transform(5, 20);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(g(i), 20);
  }
}

}  // namespace
