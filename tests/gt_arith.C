// Grid Test:  Basic arithmetic operators

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

  mygrid a("hilbert", size);
  mygrid b(size, size, 2);
  mygrid c;

  cout << "a\n" << a << endl;
  cout << "b\n" << b << endl;
  cout << "a + b\n" << a + b << endl;
  cout << "a - b\n" << a - b << endl;
  cout << "(a + 4) * b\n" << (a + 4.) * b << endl;
  cout << "a / (b + 3)\n" << a / (b + 3.) << endl;
  cout << "(a - b) * (a + b)\n" << (a - b) * (a + b) << endl;
  cout << "a^2 - b^2\n" << (a * a) - (b * b) << endl;

  c = a; c.cbind(b);
  try {
    a += c;
    a.rbind(c);
  } catch (mygrid::grid_error& g) {
    cout << "\nCaught Grid Error: [" << g.what() << "]\n";
    c = a; c.cbind(b);
  }
 
  cout << "c = a; c.cbind(b)\n" << c << endl;

  cout << "************************************************\n\n\n";

  myvector d(4);
  d[0] = 2; d[1] = 0; d[2] = -5; d[3] = 3;
  d += 2;

  cout << "a =\n" << a << endl;
  cout << "d =\n" << d << endl;
  b = a;
  a += a += d;
  cout << "a += a += d\n" << a << endl;
  a = b;
  a -= d;
  cout << "a -= d\n" << a << endl;
  a = b;
  a *= d;
  cout << "a *= d\n" << a << endl;

  return 0;
}

