// Grid Test: Compute eigenvectors/eigenvalues of some difficult matrices.

#include <cstdlib>
#include <iostream>
#include "grid.h"

using namespace std;

#ifdef DGRID
typedef dgrid mygrid;
typedef dvector myvector;
#else
typedef zgrid mygrid;
typedef zvector myvector;
#endif

int main(int argc, char** argv)
{
  int size = 3;
  if (argc > 1) size = atoi(argv[1]);

  int info;
  mygrid evecs, mg;
  dgrid evals, evalssym;

  cout << size << "x" << size << 
	  " Hilbert Matrix Eigenvectors and Eigenvalues\n";
  mg = mygrid("hilb", size);
  info = eigenvalues(mg, evals, gengrid);
  if (info) { cerr << "Error computing eigenvalues\n"; exit(1); }
  cout << "eigenvalues(gengrid)\n" << evals << endl;
  info = eigenvalues(mg, evalssym, symgrid);
  if (info) { cerr << "Error computing eigenvalues\n"; exit(1); }
  cout << "eigenvalues(symgrid)\n" << evalssym << endl;
  cout << "Maxabs diff between the two set of eigenvalues: "
       << linfvecnorm(evals - evalssym) << endl;
  info = eigenvectors(mg, evecs, evals);
  if (info) { cerr << "Error computing eigenvectors\n"; exit(1); }
  cout << endl << evecs << endl << evals << endl;
  cout << "max(abs(Original * Eigenvectors - Eigenvectors * Eigenvalues)) = "
       << linfvecnorm(matmult(mg, evecs) - diagmult(evecs, (const dvector&)evals)) << endl;

  cout << "\n" << size << "x" << size 
	  << " Wilkinson Matrix Eigenvectors and Eigenvalues\n";
  mg = mygrid("wilk", size);
  info = eigenvalues(mg, evals, gengrid);
  if (info) { cerr << "Error computing eigenvalues\n"; exit(1); }
  cout << "eigenvalues(gengrid) - absolute values since the eigenvalues\n"
       << "are returned in a real grid but could be complex\n" << evals << endl;
  info = eigenvalues(mg, evalssym, symgrid);
  if (info) { cerr << "Error computing eigenvalues\n"; exit(1); }
  cout << "eigenvalues(symgrid)\n" << evalssym << endl;
  info = eigenvectors(mg, evecs, evals);
  if (info) { cerr << "Error computing eigenvectors\n"; exit(1); }
  cout << endl << evecs << endl << evals << endl;
  cout << "max(abs(Original * Eigenvectors - Eigenvectors * Eigenvalues)) = "
       << linfvecnorm(matmult(mg, evecs) - diagmult(evecs, (const dvector&)evals))
       << endl;

  return 0;
}

