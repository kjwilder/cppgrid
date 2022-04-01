// Grid Test: Compute eigenvectors/eigenvalues of some difficult matrices.

#include <cstdlib>
#include <iostream>
#include "grid.h"
#include "gtest/gtest.h"

namespace {

#ifdef DGRID
typedef dgrid mygrid;
typedef dvector myvector;
#else
typedef zgrid mygrid;
typedef zvector myvector;
#endif

TEST(Eigenvalues, All) {
  int size = 3;

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
       << linfvecnorm(
           matmult(mg, evecs) - diagmult(evecs, (const dvector&)evals)) << endl;

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
       << linfvecnorm(matmult(mg, evecs) - diagmult(
             evecs, (const dvector&)evals))
       << endl;
}

};  // namespace
