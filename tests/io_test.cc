// Grid Test: I/O Examples with grids and vectors of grids.

#include <iostream>
#include <string>
#include "grid.h"

namespace {

#ifdef DGRID
typedef dgrid mygrid;
typedef dgrids mygrids;
#else
typedef zgrid mygrid;
typedef zgrids mygrids;
#endif

int main(int argc, char** argv) {
  int size = 3;
  if (argc > 1) size = atoi(argv[1]);

  string file1 = "j1", file2 = "j2";

  mygrid a("hilbert", size);
  cout << "A =\n" << a << endl;

  if (!a.write(file1)) {
    cerr << "Unable to write a grid to file " << file1 << endl;
    exit(1);
  }

  mygrid b;
  if (!b.read(file1)) {
    cerr << "Unable to read a grid from file " << file1 << endl;
    exit(1);
  }
  cout << "B obtained by reading file " << file1 << " written to by A\n";
  cout << "B =\n" << b << endl;

  mygrids c(2), d(2);
  for (int i = 0; i < 2; ++i)
    c[i] = mygrid("hilbert", 4 * i + 3);

  cout << "Now, we create a vector of two grids, write them to a single\n"
       << "file, and then read them into another vector of two grids.\n";

  ofstream ofs(file2.c_str());
  if (!ofs) {
    cerr << "Unable to write a grids to file " << file2 << endl;
    exit(1);
  }
  for (int i = 0; i < 2; ++i)
    c[i].write(ofs, 1);
  ofs.close();

  ifstream ifs(file2.c_str());
  if (!ifs) {
    cerr << "Unable to read grids from file " << file2 << endl;
    exit(1);
  }
  for (int i = 0; i < 2; ++i)
    d[i].read(ifs, 1);
  for (int i = 0; i < 2; ++i)
    cout << "Grid d(" << i << "):\n" << d[i] << endl;

  return 0;
}

};  // namespace
