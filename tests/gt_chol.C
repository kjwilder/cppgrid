// gt_chol : Test Cholesky functions.

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
typedef vector<mytype> myvector;

int main(int argc, char** argv)
{
  int info, size = 3;
  if (argc > 1) size = atoi(argv[1]);

  // Create a symmetric positive definite matrix
  mygrid a("hilb", size);
  cout << "A = " << size << "x" << size << " Hilbert matrix\n" << a << endl;

  mygrid a1;
  info = chol(a, a1);
  cout << "Upper Cholesky Factor of A =\n" << a1 << endl;

  myvector b1(size), b2;
  for (int i = 0; i < size; ++i) b1[i] = i;

  mygrid a2 = a1;
  b2 = b1;
  cout << "B =\n" << b2 << endl;
  cholupdate(a2, b2);
  cout << "Upper Cholesky Factor of A after updating with B =\n" << a2 << endl;

  a2 = a1;
  b2 = b1 * 2 + 0.5;
  cout << "B =\n" << b2 << endl;
  cholupdate(a2, b2);
  cout << "Upper Cholesky Factor of A after updating with B =\n" << a2 << endl;


  return 0;
}

