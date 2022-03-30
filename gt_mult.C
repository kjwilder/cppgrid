// Grid Test 3:
// Matrix multiplications with vectors and grids.

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
  uint size = 3;
  if (argc > 1) size = atoi(argv[1]);

  mygrid d("hilbert", size);
  d.resize(size + 1, size - 1);

  cout << "D\n" << d << endl;

  myvector dvleft = d, dvright = d;
  dvleft.resize(size+1); dvright.resize(size-1);
  myvector vres;
  mygrid dres;

  cout << "Vector results\n";
  cout << "D %*% D[1:(sz-1)]\n" << matmult(d, dvright, vres) << endl;
  cout << "D[1:(sz+1)] %*% D\n" << matmult(dvleft, d, vres, 't') << endl;
  cout << "Grid results\n";
  cout << "D %*% D[1:(sz-1)]\n" << matmult(d, dvright, dres) << endl;
  cout << "D[1:(sz+1)] %*% D\n" << matmult(dvleft, d, dres, 't') << endl;
  cout << "D.diagmult(D[1:(sz-1)])\n" << diagmult(d, dvright, dres) << endl;
  cout << "D.diagmult(D[1:(sz+1)], left)\n" << diagmult(dvleft, d, dres) << endl;

  mygrid e(size - 2, size - 1, 0);
  myvector v(size);
  for (uint i = 0; i < size; ++i) v[i] = i + 1; 
  e += v;
  cout << "mygrid(sz-2, sz-1, 0) + myvector(c(1,...,sz))\n" << e << endl;

  // The most efficient way to compute C = alpha * A %*% B + beta * C
  mygrid m1("seq", 1, size * size + size, 1);
  m1.resize(size, size + 1);
  mygrid m2("seq", 11, 10 + size * size + size, 1);
  m2.resize(size + 1, size);
  mygrid m3("seq", 21, 20 + size * size, 1);
  m3.resize(size, size);
  cout << "m1\n" << m1 << "m2\n" << m2 << "m3\n" << m3 << endl;
  cout << "2 * m1 * m2 + 5 * m3\n" << matmult(m1, m2) * 2 + m3 * 5 << endl;
#ifdef USE_LAPACK
  gemm(m1, m2, m3, 2, 5);
  cout << "gemm: m3 = 2 * m1 * m2 + 5 * m3\n" << m3 << endl;
#endif

  m1.resize(2, 2); m1[0] = 0; m1[1] = 1; m1[2] = 2; m1[3] = 3;
  m2 = m1 + 1;
  cout << "Different ways to do a transpose\n";
  cout << "m1=\n" << m1 << endl;
  cout << "m2=\n" << m2 << endl;
  cout << "matmult(m1, trans(m2), m3)\n" <<  matmult(m1, trans(m2), m3) << endl;
  cout << "matmult(m1, m2, m3, 'N', 'T')\n" <<  matmult(m1, m2, m3, 'N', 'T') << endl;
  cout << "matmult(trans(m1), m2, m3)\n" <<  matmult(trans(m1), m2, m3) << endl;
  cout << "matmult(m1, m2, m3, 'T')\n" <<  matmult(m1, m2, m3, 'T') << endl;

#ifdef USE_LAPACK
  m1 = m2 = mygrid("hilbert", size);
  m3.resize(m1.rows(), m2.cols());
  cout << "Results of triangular matrix multiplications\n";
  cout << "m1=\n" << m1 << endl;
  cout << "m2=\n" << m2 << endl;
  cout << "m2 = lower_tri(m1) %*% m2\n" << trmm(m1, m2,1,'L','L') << endl;
  m2 = mygrid("hilbert", size);
  cout << "m2 = unit_diag_lower_tri(m1) %*% m2\n" << trmm(m1, m2,1,'L','L','N','U') << endl;
#endif

  mygrid ns1("seq",1,8,1), ns2("seq",11,18,1), ns3;
  ns1.resize(2,4);
  ns2.resize(2,4);
  cout << "ns1 = \n" << ns1 << endl;
  cout << "ns2 = \n" << ns2 << endl;
  cout << "ns1 * ns2 (tn) = \n" << matmult(ns1, ns2, ns3, 't', 'n') << endl;
  cout << "ns1 * ns2 (nt) = \n" << matmult(ns1, ns2, ns3, 'n', 't') << endl;
  cout << "ns1 * ns2 (cn) = \n" << matmult(ns1, ns2, ns3, 'c', 'n') << endl;
  cout << "ns1 * ns2 (nc) = \n" << matmult(ns1, ns2, ns3, 'n', 'c') << endl;
  return 0;
}

