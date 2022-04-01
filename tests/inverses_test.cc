// Grid Test: Inverse functions.

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

TEST(Inverses, All) {
  int size = 3;

  // Create a symmetric positive definite matrix
  mygrid a("hilbert", size), a_lu, a_inv;
  cout << "A = " << size << "x" << size << " Hilbert matrix\n" << a << endl;

  // Compute the inverse of A using "inv"
  cout << "gen inverse computed using inv():\n" << inv(a) << endl;
  cout << "sym inverse computed using inv():\n" << inv(a, symgrid) << endl;
  cout << "pos inverse computed using inv():\n" << inv(a, posgrid) << endl;

  // Compute the inverse of A using "LU/LUsolve"
  ivector pivots(size);
  int info = LU(a, a_lu, pivots);
  if (info) { cerr << "Error in LU: info = " << info << endl; }
  LUsolve(a_lu, mygrid("I", size), a_inv, pivots);
  if (info) { cerr << "Error in LUsolve: info = " << info << endl; }
  cout << "A inverse computed using LUsolve:\n" << a_inv << endl;

  info = LU(a, a_lu, pivots, symgrid);
  if (info) { cerr << "Error in LU(symgrid): info = " << info << endl; }
  LUsolve(a_lu, mygrid("I", size), a_inv, pivots, symgrid);
  if (info) { cerr << "Error in LUsolve(symgrid): info = " << info << endl; }
  cout << "A inverse computed using LUsolve(symgrid):\n" << a_inv << endl;

  info = LU(a, a_lu, pivots, posgrid);
  if (info) { cerr << "Error in LU(posgrid): info = " << info << endl; }
  LUsolve(a_lu, mygrid("I", size), a_inv, pivots, posgrid);
  if (info) { cerr << "Error in LUsolve(posgrid): info = " << info << endl; }
  cout << "A inverse computed using LUsolve(posgrid):\n" << a_inv << endl;

  // Compute the inverse of A using "solve"
  solve(a, mygrid("I", size), a_inv);
  cout << "gen inverse computed using solve():\n" << a_inv << endl;
  solve(a, mygrid("I", size), a_inv, symgrid);
  cout << "sym inverse computed using solve():\n" << a_inv << endl;
  solve(a, mygrid("I", size), a_inv, posgrid);
  cout << "pos inverse computed using solve():\n" << a_inv << endl;

  // Create a symmetric, but not positive definite, matrix
  a = mygrid("wilkinson", size);
  cout << "A = " << size << " by " << size << " Wilkinson matrix =\n"
       << a << endl;

  // Compute the inverse of A using "inv"
  cout << "gen inverse computed using inv():\n" << inv(a) << endl;
  cout << "sym inverse computed using inv():\n" << inv(a, symgrid) << endl;
  cout << "pos inverse computed using inv(): (Should be incorrect)\n"
       << inv(a, posgrid) << endl;

  // Do everything with a triangular grid.
  for (int i = 0; i != size; ++i)
    for (int j = 0; j != size; ++j)
      a(i, j) = (i >= j ? i - j + 1 : 0);
  cout << "A =\n" << a << endl;
  cout << "gen inverse computed using inv():\n" << inv(a) << endl;
  cout << "sym inverse computed using inv(): (Should be incorrect)\n"
       << inv(a, symgrid, 'L') << endl;
  cout << "pos inverse computed using inv(): (Should be incorrect)\n"
       << inv(a, posgrid, 'L') << endl;
  cout << "tri inverse computed using inv():\n" << inv(a, trigrid, 'L') << endl;

  return 0;
}

};  // namespace
