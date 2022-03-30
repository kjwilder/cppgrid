// Grid Test: Compute various norms.

#include <iostream>
#include "grid.h"

using namespace std;

#ifdef DGRID
typedef dgrid mygrid;
#else
typedef zgrid mygrid;
#endif

int main(int argc, char** argv)
{
  uint size = 3;
  if (argc > 1) size = atoi(argv[1]);

  mygrid d("hilb", size), e = mygrid("seq", 1, size * size, 1), f;
  e.resize(size, size);
  e += d;
#ifndef DGRID
  e += cplx(0, 1);
#endif
  cout << "d =\n" << d << endl;

  cout << "L1 vector norm is " << l1vecnorm(d) << endl;
  cout << "L2 vector norm is " << l2vecnorm(d) << endl;
  cout << "Linf vector norm is " << linfvecnorm(d) << endl;
  cout << "L1 matrix norm is " << l1matnorm(d) << endl;
  cout << "L2 matrix norm is " << l2matnorm(d, symgrid) << endl;
  cout << "Linf matrix norm is " << linfmatnorm(d) << endl;
  cout << "Frobenius matrix norm is " << frobnorm(d) << endl;

  return 0;
}

