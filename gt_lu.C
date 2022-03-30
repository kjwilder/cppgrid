// Grid Test 9:
// Tests positive definite matrix functions, back/forward solve.

#include <cstdlib>
#include <iostream>
#include "grid.h"

using namespace std;

#ifdef DGRID
typedef double mytype;
#else
typedef cplx mytype;
#endif
typedef grid<mytype> mygrid;

int main(int argc, char** argv)
{
  int info, size = 3;
  if (argc > 1) size = atoi(argv[1]);

  // Create a symmetric positive definite matrix
  mygrid a("hilb", size);
  cout << "A = " << size << "x" << size << " Hilbert matrix\n" << a << endl;

  // Compute the upper (or lower) Cholesky factor of A
  mygrid a_chol;
  info = chol(a, a_chol);
  //info = chol(a, a_chol, 'L');
  if (info) { cerr << "Error computing Cholesky factor: info = "
	           << info << endl; exit(1); }
  cout << "a_chol = Upper Cholesky factor of A\n" << a_chol << endl;

  // Verify that t(B) %*% B and A are the same
  mygrid c = matmult(trans(a_chol), a_chol) - a;
#ifdef DGRID
  cout << "max(abs(t(a_chol) %*% a_chol - A)) = " << linfvecnorm(c) << "\n\n";
#endif

  // Compute solve(A) using the Cholesky factor of A
  mygrid a_inv;
  info = cholsolve(a_chol, mygrid("I", a_chol.rows()), a_inv);
  if (info) { cerr << "Error computing Cholesky inverse\n"; exit(1); }
  cout << "a_inv = solve(A) computed using a_chol\n" << a_inv << endl;

  // Compute solve(A) more directly.
  mygrid a_inv2 = inv(a, posgrid);
  cout << "a_inv2 = solve(A) using inv(posgrid)\n"
#ifdef DGRID
       << "max(abs(a_inv2 - a_inv)) = " << linfvecnorm(a_inv2 - a_inv) << endl
#endif
       << endl;

  // Compute solve(A) using back/forward subsitution with the
  // Cholesky factor of A.  Note that we avoid having to (explicitly)
  // transpose the upper Cholesky factor of A by supplying the
  // option 'U' (upper)  argument to forwardsolve.  Similarly,
  // backsolve has an optional 'L' argument.  It works fine and it
  // not particularly inefficient to use 'h' in both functions.
  backsolve(a_chol, forwardsolve(a_chol, mygrid("I", size), a_inv2, 'U'), a_inv2);
  cout << "a_inv2 = solve(A) using back/forward substitution with B\n"
#ifdef DGRID
       << "max(abs(a_inv2 - a_inv)) = " << linfvecnorm(a_inv2 - a_inv) << endl
#endif
       << endl;

  // Create an invertible square matrix
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      a(i, j) = (7 * i + 3 * j) / (10.0 * size);
  for (int i = 0; i < size; ++i)
    a(i, i) = size;
  cout << "A\n" << a << endl;

  // Create a right-hand-side matrix with three columns.
  mygrid res(size, 3);
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < 3; ++j)
      res(i, j) = (7 * i + 3 * j) / (10.0 * size);
  cout << "res\n" << res << endl;

  mygrid lu;
  ivector pivots;
  LU(a, lu, pivots);
  cout << "LU decomp of A with pivots\n" << lu << pivots << endl;

  // Tests:
  mygrid x;
  backsolve(lu, res, x);
  cout << "backsolve(lu, res, x)\n" << x << endl;
  backsolve(trans(lu), res, x, 'L');
  cout << "backsolve(trans(lu), res, x, 'L')\n" << x << endl;
  forwardsolve(lu, res, x);
  cout << "forwardsolve(lu, res, x)\n" << x << endl;
  forwardsolve(trans(lu), res, x, 'U');
  cout << "forwardsolve(trans(lu), res, x, 'U')\n" << x << endl;

  dgrid r_mat("hilbert", size);
  cout << "r_mat=\n" << r_mat << endl;
  dgrid cstar;
  info = chol(r_mat,cstar,'L');
  cout << "cstar=chol(r_mat)=\n" << cstar << endl;
  if (info) {
    std::cerr << "Error computing Cholesky factor: info = " << info << "\n";
    exit(1);
  }
  inv(cstar, cstar, trigrid, 'L');
  cout << "inv(cstar)=\n" << cstar << endl;
  chol(r_mat, r_mat, 'L');
  cout << "inv(cstar) %*% cstar=\n" << matmult(cstar, r_mat)  + 1e-10 << endl;

  return 0;
}

