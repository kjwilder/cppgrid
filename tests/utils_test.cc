#include "utils.h"
#include "gtest/gtest.h"

namespace {

using utils_h::distance;
using utils_h::distance_squared;
using utils_h::linfdist;
using utils_h::angle;

TEST(Utils, Distance) {
  EXPECT_EQ(5.0, distance(0, 0, 3, 4));
}

TEST(Utils, Distance_squared_int) {
  EXPECT_EQ(25, distance_squared(0, 0, 3, 4));
}

TEST(Utils, Distance_squared_float) {
  EXPECT_EQ(25.0, distance_squared(0.0, 0.0, 3.0, 4.0));
}

TEST(Utils, Linfdist_int) {
  EXPECT_EQ(4, linfdist(0, 0, 3, 4));
}

TEST(Utils, Linfdist_float) {
  EXPECT_EQ(4.0, linfdist(0.0, 0.0, 3.0, 4.0));
}

TEST(Utils, Angle) {
  EXPECT_NEAR(-0.927, angle(0, 0, 3, 4), 0.001);
}

}  // namespace
