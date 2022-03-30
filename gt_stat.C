// Grid Test: Statistical Functions

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

  mygrid a0("hilbert", size);
  a0.cbind(a0);  // Test cbinding a grid to itself.
  mygrid a;
  a = a0;  // Test copying
  cout << "A =\n" << a << endl;

  mygrid b0("hilbert", size);
  b0.rbind(b0);  // Test rbinding a grid to itself.
  a = b0;  // Test copying
  cout << "A =\n" << a << endl;

  myvector sumrows(a.rows()), sumcols(a.cols());
  cout << "A row sums\n";
  for (uint i = 0; i < a.rows(); i++)
    sumrows[i] = sum(a, "row", i);
  cout << sumrows;
  cout << "A column sums\n";
  for (uint i = 0; i < a.cols(); i++)
    sumcols[i] = sum(a, "col", i);
  cout << sumcols;

  myvector means(a.rows());
  cout << "A row means\n";
  for (uint i = 0; i < a.rows(); i++)
    means[i] = mean(a, "row", i);
  cout << means;
  means.resize(a.cols());
  cout << "A column means\n";
  for (uint i = 0; i < a.cols(); i++)
    means[i] = mean(a, "col", i);
  cout << means;

  mygrid cv;
  cout << "\nCol covariances =\n" << cov(a, cv) << endl;
  cout << "Row covariances =\n" << cov(a, cv, "row") << endl;
  cor(a, cv);
  cout << "Col correlations =\n" << cv << endl;
  cor(a, cv, "row");
  cout << "Row correlations =\n" << cv << endl;

  return 0;
}

