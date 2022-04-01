// Grid Test: Test standard math functions.
// One can assign, e.g., 'y' to be the sine of 'x', using either:
//   sin(x, y)
//   y = sin(x)
// The first version is preferred as it is more efficient.
// Both 'sin(x, x)' and 'x = sin(x)' will work as expected.

#include <cstdlib>
#include <iostream>
#include "grid.h"

namespace {

#ifdef DGRID
typedef dgrid mygrid;
typedef dvector myvector;
#else
typedef zgrid mygrid;
typedef zvector myvector;
#endif

int main(int argc, char** argv) {
  int size = 3;
  if (argc > 1) size = atoi(argv[1]);

  std::cout << std::setiosflags(std::ios::fixed);

  mygrid d(size), f, g;
  myvector vd, vf, vg;
  for (uint i = 0; i < d.size(); ++i)
    d[i] = i + 1;
  cout << "d = \n" << trans(d) << trans(exp(log(d)));

  exp(d, f); log(f, f); g = log(exp(d));
  cout << "log(exp(d)) = \n" << trans(f) << trans(g);
  for (uint i = 0; i != d.size(); ++i)
    cout << setw(10) << setprecision(4) << log(exp(d[i])) << " "; cout << endl;

  vd = d;
  sin(d, f); g = sin(d);
  sin(d, vf); vg = sin(d);
  cout << "sin(d) = \n" << trans(f) << trans(g);
  cout << vf << vg;
  sin(vd, f); g = sin(vd);
  sin(vd, vf); vg = sin(vd);
  cout << trans(f) << trans(g);
  cout << vf << vg;
  for (uint i = 0; i != d.size(); ++i)
    cout << setprecision(4) << setw(10) << sin(d[i]) << " "; cout << endl;

#ifdef DGRID
  atan(d, f); g = atan(d);
  cout << "atan(d) =\n" << trans(f) << trans(g);
  for (uint i = 0; i != d.size(); ++i)
    cout << setw(10) << setprecision(4) << atan(d[i]) << " "; cout << endl;

  tanh(d, f); g = tanh(d);
  cout << "tanh(d) =\n" << trans(f) << trans(g);
  for (uint i = 0; i != d.size(); ++i)
    cout << setw(10) << setprecision(4) << tanh(d[i]) << " "; cout << endl;
#endif

#ifndef DGRID
  pow(d, cplx(2.2), f); g = pow(d, cplx(2.2));
  cout << "pow(d,cplx(2.2)) =\n" << trans(f) << trans(g);
#endif
  return 0;
}

};  // namespace
