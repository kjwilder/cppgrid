// Grid Test: Compute eigenvectors/eigenvalues of the difficult Rosser Matrix

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

int main()
{
  int info;
  mytype rosser_values[] = 
    {611, 196, -192, 407, -8, -52, -49, 29, 196, 899, 113, -192, -71, -43,
     -8, -44, -192, 113, 899, 196, 61, 49, 8, 52, 407, -192, 196, 611, 8,
     44, 59, -23, -8, -71, 61, 8, 411, -599, 208, 208, -52, -43, 49, 44,
     -599, 411, 208, 208, -49, -8, 8, 59, 208, 208, 99, -911, 29, -44, 52,
     -23, 208, 208, -911, 99};

  mygrid evecs, mg;
  dgrid evals;

  mygrid rosser(rosser_values, 8, 8);
  cout << "(Difficult) 8x8 Rosser Matrix\n" << rosser << endl;

  info = eigenvectors(rosser, evecs, evals);
  if (info) { cerr << "Error computing eigenvectors\n"; exit(1); }
  cout << "Rosser Matrix Eigenvectors\n" << evecs << endl;
  cout << "Rosser Matrix Eigenvalues\n" << evals << endl;
  cout << "max(abs(rosser * evecs - evecs * evals)) is "
       << linfvecnorm(matmult(rosser, evecs) - diagmult(evecs, (const dvector&)evals))
       << endl;

  return 0;
}

