// Grid Test: Determinant functions.

#include <cstdlib>
#include <iostream>
#include "grid.h"
#include "gtest/gtest.h"

namespace {

#ifdef DGRID
typedef dgrid mygrid;
#else
typedef zgrid mygrid;
#endif

TEST(Determinant, All) {
  int size = 3;

  // Create a symmetric positive definite matrix
  mygrid a("hilbert", size);
  cout << "A = " << size << "x" << size << " Hilbert matrix\n" << a << endl;

  // Compute the determinant of A using various LU decompositions.
  cout << "gen determinant: " << det(a) << endl;
  cout << "sym determinant: " << det(a, symgrid) << endl;
  cout << "pos determinant: " << det(a, posgrid) << endl;
  cout << "tri determinant: (Should be wrong) " << det(a, trigrid) << endl;

  // Create a symmetric NEGATIVE!! definite matrix
  mygrid b(3, 3, 1);
  b(1, 0) = -1.0; b(0, 1) = -1.0;
  b(2, 0) = 0.494494; b(0, 2) = 0.494494;
  b(1, 2) = 0.474852; b(2, 1) = 0.474852;
  cout << "\n\nB =\n" << b << endl;

  // Compute the determinant of B using various LU decompositions.
  cout << "gen determinant: " << det(b) << endl;
  cout << "sym determinant: (Should NOT be wrong) " << det(b, symgrid) <<
endl;
  cout << "pos determinant: (Should be wrong) " << det(b, posgrid) <<
endl;
  cout << "tri determinant: (Should be wrong) " << det(b, trigrid) <<
endl;
}

};  // namespace
