#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <cstdio>
#include <cassert>
#include <fstream>
#include <vector>

// __________________________________________________________________________
// grid of elements of type T

template <class T>
class grid
{

 private:
  int nr, nc; // number of rows, columns
  T *sto;

  // Initialization Member Functions
  void initnew(int r, int c) { 
    assert(r >= 0 && c >= 0); 
    nr = r; nc = c; sto = (r == 0 || c == 0) ? 0 : new T[r * c];
  }
  void freegrid() { assert(invariant()); nr = nc = 0; delete [] sto; sto = 0; }

  // Consistency Checking Functions
  int inrange(int r, int c) const { 
    return (r >= 0 && c >= 0 && r < nr && c < nc && sto != 0); 
  }
  int invariant() const { 
    return (((nr == 0 || nc == 0) && sto == 0 && nr >= 0 && nc >= 0) || 
	    ((nr > 0 && nc > 0) && sto != 0));
  }

 public:
  // Constructors, Destructor, Operator=
  grid() { initnew(0, 0); }
  grid(int r) { initnew(r, 1); }
  grid(int r, int c) { initnew(r, c); }
  grid(const grid &m) : nr(0), nc(0), sto(0) { *this = m; }
  grid<T>& operator=(const grid &m);
  void subgrid(grid& m, int r, int c, int numrows, int numcols) const;

  ~grid() { assert (invariant()); delete [] sto; }

  // Basic Member Access Functions
  int rows() const { return nr; }
  int cols() const { return nc; }
  T& operator()(int r) const { 
    assert(inrange(r, 0)); return sto[r]; 
  }
  T& operator()(int r, int c) const { 
    assert(inrange(r, c)); return sto[c * nr + r]; 
  }
  T& operator[](int r) const { return sto[r]; }
  void set(int r, int c, T val) { if (inrange(r, c)) (*this)(r, c) = val; }
  void set(double r, double c, T val) { set(int(r), int(c), val); }
  void set(int r, double c, T val) { set(r, int(c), val); }
  void set(double r, int c, T val) { set(int(r), c, val); }
  T get(int r, int c) const { return (inrange(r, c) ? (*this)(r, c) : 0); }
  T get(double r, double c) const { return get(int(r), int(c)); }
  T get(int r, double c) const { return get(r, int(c)); }
  T get(double r, int c) const { return get(int(r), c); }

  // Initialization Member Functions
  void init() { init(0, 0); }
  void init(int r) { init(r, 1); }
  void init(int r, int c) { 
    assert(invariant());
    if (r > 0 && c > 0 && r * c <= nr * nc && r * c * 2 >= nr * nc) 
    { nr = r; nc = c; }
    else 
    { freegrid() ; initnew(r, c); } 
    assert(invariant()); 
  }
  void clear(const T& val = 0)
  { assert(invariant()); for (int i = 0; i < nr * nc; ++i) sto[i] = val; }

  // I/O Functions
  void write(const char *file);
  int write(std::ofstream& os);
  int read(const char* file);
  int read(std::ifstream& is);
  int loadpgm(const char* pgmname);
  int savepgm(const char* pgmname);
  int size() const;
  void dump(int max = -1) const;
  void dump2(int max = -1) const;

  // Useful Utility Functions
  grid<T>& operator<<(grid &m);
  grid<T>& operator+=(const grid &m);
  grid<T>& operator-=(const grid &m);
  grid<T>& operator+=(T i);
  grid<T>& operator-=(T i) { return *this += -i;  }
  grid<T> operator*(const grid &m) const;
  void normalize(T val);
  void normalize(T minval, T maxval);
  const grid<T> transpose() const;
  grid<T> LU() const;
  grid<T> inverse() const;
  int offpixels() { 
    int c = 0; for (int i = 0; i < nr * nc; ++i) c += (sto[i] == 0); return c;
  }
  int onpixels() { return(nr * nc - offpixels()); }
  void sort(int col = 0, int left = -1, int right = -1);

}; // class grid

// __________________________________________________________________________
// Dump a grid; Will only work for classes with ostream<<(const T&).

template <class T>
void grid<T>::dump(int max) const
{
  if (max == -1)
    max = nc;
  for (int j = 0; j < max; ++j) {
    for (int i = 0; i < nr; ++i)
      std::cout << ((*this)(i, j)) << " ";
    std::cout << "\n";
  }

} // grid::dump

// __________________________________________________________________________
// Dump a grid; Will only work for classes with ostream<<(const T&).

template <class T>
void grid<T>::dump2(int max) const
{
  if (max == -1)
    max = nr;
  for (int i = 0; i < max; ++i) {
    for (int j = 0; j < nc; ++j)
      std::cout << ((*this)(i, j)) << " ";
    std::cout << "\n";
  }

} // grid::dump2

// __________________________________________________________________________
// Write out a grid to a file

template<class T>
void grid<T>::write(const char *file)
{
  assert(invariant());
  Ofstream(ofs, file);
  if (!ofs) {
    std::cerr << "Unable to write a grid to file [" << file << "]." << std::endl;
    return;
  }

  ofs.write("GR11", 4);
  varwrite(ofs, nr);
  varwrite(ofs, nc);
  arraywrite(ofs, sto, nr * nc);

} // grid::write

// __________________________________________________________________________
// Write out a grid to an open output stream

template<class T>
int grid<T>::write(std::ofstream& ofs)
{
  assert(invariant());

  //if (!ofs.is_open())
  if (!ofs)
    return 0;

  varwrite(ofs, nr);
  varwrite(ofs, nc);
  arraywrite(ofs, sto, nr * nc);

  return 1;

} // grid::write

// __________________________________________________________________________
// Calculate the size of a grid which is to be written out

template<class T>
int grid<T>::size() const
{
  assert(invariant());
  return (sizeof(nr) + sizeof(nc) + nr * nc * sizeof(T));

} // grid::size

// __________________________________________________________________________
// Read in a grid which was previously written to a file

template<class T>
int grid<T>::read(const char *file)
{
  assert(invariant());
  Ifstream(ifs, file);
  if (!ifs)
    return 0;

  char version[4];
  ifs.read(version, 4);
  if (memcmp(version, "GR11", 4) != 0 && memcmp(version, "GR12", 4) != 0)
  {
    if (loadpgm(file))
      return 1;
    else {
      std::cerr << "The file [" << file << "] is not a grid or pgm file" << std::endl;
      return 0;
    }
  }

  freegrid();
  if (!memcmp(version, "GR11", 4))
  {
    varread(ifs, nr);
    varread(ifs, nc);
    assert(nr >= 0 && nc >= 0);
    if (nr > 0 && nc > 0) {
      sto = new T[nr * nc];
      arrayread(ifs, sto, nr * nc);
    }
  }
  else if (!memcmp(version, "GR12", 4))
  {
    ifs >> nr;
    ifs >> nc;
    assert(nr >= 0 && nc >= 0);
    if (nr > 0 && nc > 0) {
      sto = new T[nr * nc];
      for (int j = 0; j < nc; ++j)
	for (int i = 0; i < nr; ++i)
          ifs >> (*this)(i, j);
    }
  }
  assert(invariant());

  return 1;

} // grid::read

// __________________________________________________________________________
// Read in a grid from an open istream

template<class T>
int grid<T>::read(std::ifstream& is)
{
  assert(invariant());

  //if (!is.is_open() || is.eof())
  if (!is || is.eof())
    return 0;

  freegrid();
  varread(is, nr);
  varread(is, nc);
  assert(nr >= 0 && nc >= 0);
  if (nr > 0 && nc > 0) {
    sto = new T[nr * nc];
    arrayread(is, sto, nr * nc);
  }

  assert(invariant());

  return 1;

} // grid::read

// __________________________________________________________________________
// Set a grid equal to a grid 'm'; leave 'm' unchanged.

template<class T>
grid<T>& grid<T>::operator=(const grid& m)
{
  assert(invariant() && m.invariant());
  if (this != &m) {
    init(m.nr, m.nc);
    if (nr > 0 && nc > 0) 
      memcpy(sto, m.sto, nr * nc * sizeof(T));
  }
  assert(invariant() && m.invariant());

  return *this;

} // grid::operator=

// __________________________________________________________________________
// Set 'm' equal to a subgrid of *this.

template<class T>
void grid<T>::subgrid(grid& m, int r, int c, int numrows, int numcols) const
{
  assert(invariant() && m.invariant());
  assert(r >= 0 && c >= 0 && numrows >= 0 && numcols >= 0);
  assert(r + numrows <= rows() && c + numcols <= cols());

  if (this != &m) {
    m.init(numrows, numcols);
    if (nr > 0 && nc > 0) 
      for (int i = 0; i < m.nr; ++i)
	for (int j = 0; j < m.nc; ++j)
	  m(i, j) = (*this)(r + i, c + j);
  }
  else {
    grid tmp;
    subgrid(tmp, r, c, numrows, numcols);
    m << tmp;
  }

  assert(invariant() && m.invariant());

} // grid::subgrid

// __________________________________________________________________________
// Set a grid equal to 'm'; obliterate 'm'.

template<class T>
grid<T>& grid<T>::operator<<(grid &m)
{
  assert(invariant() && m.invariant());
  if (this != &m) {
    freegrid(); 
    nr = m.nr; nc = m.nc; sto = m.sto;
    m.initnew(0, 0);
  }
  assert(invariant() && m.invariant());

  return *this;

} // grid::operator<<

// __________________________________________________________________________
// Add another grid to the current grid.

template<class T>
grid<T>& grid<T>::operator+=(const grid& m)
{
  assert(invariant() && m.invariant());
  assert(nr == m.nr && nc == m.nc);

  for (int i = 0; i < nr * nc; ++i)
    sto[i] += m.sto[i];

  assert(invariant() && m.invariant());

  return *this;

} // grid::operator+=

// __________________________________________________________________________
// Subtract a grid from the current grid.

template<class T>
grid<T>& grid<T>::operator-=(const grid& m)
{
  assert(invariant() && m.invariant());
  assert(nr == m.nr && nc == m.nc);

  for (int i = 0; i < nr * nc; ++i)
    sto[i] -= m.sto[i];

  assert(invariant() && m.invariant());

  return *this;

} // grid::operator-=

// __________________________________________________________________________
// Add a fixed value to each element of a grid

template<class T>
grid<T>& grid<T>::operator+=(T ii)
{
  for (int i = 0; i < nr * nc; ++i)
    sto[i] += ii;
  return *this;

} // grid::operator+=

// __________________________________________________________________________

template<class T>
void grid<T>::normalize(T val)
{
  int i, j;
  T gridmax = 0;

  assert(invariant());
  for (i = 0; i < nr; ++i)
    for (j = 0; j < nc; ++j)
      if ((*this)(i, j) > gridmax)
        gridmax = (*this)(i, j);
  if (gridmax > 0) {
    for (i = 0; i < nr; ++i)
      for (j = 0; j < nc; ++j)
        (*this)(i, j) = (*this)(i, j) * val / gridmax;
  }
  assert(invariant());

} // grid::normalize

// __________________________________________________________________________
// Perform a linear transformation on the values of a grid so that the
// values range from 'val1' to 'maxval'

template<class T>
void grid<T>::normalize(T val1, T val2)
{
  if (rows() == 0 || cols() == 0)
    return;

  T gridmin = this->sto[0];
  T gridmax = this->sto[0];

  assert(invariant());
  int i, j;
  for (i = 0; i < nr; ++i)
    for (j = 0; j < nc; ++j) {
      if ((*this)(i, j) > gridmax)
        gridmax = (*this)(i, j);
      if ((*this)(i, j) < gridmin)
        gridmin = (*this)(i, j);
    }

  T range = gridmax - gridmin;
  if (range > 0) {
    for (i = 0; i < nr; ++i)
      for (j = 0; j < nc; ++j)
        (*this)(i, j) = ((*this)(i, j) * (val2 - val1) +
			 (gridmax * val1 - gridmin * val2)) / range;
  }
  else if (gridmin < val1)
    clear(val1);
  else if (gridmax > val2)
    clear(val2);

  assert(invariant());

} // grid::normalize

// __________________________________________________________________________
// Multiply two grids together using matrix multiplication.

template<class T>
grid<T> grid<T>::operator*(const grid &m) const
{
  assert(invariant() && m.invariant());
  assert(nc == m.nr);

  grid<T> tmp(nr, m.nc);
  tmp.clear();
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < m.nc; ++j)
      for (int k = 0; k < nc; ++k)
	tmp(i, j) += (*this)(i, k) * m(k, j);

  assert(invariant() && m.invariant() && tmp.invariant());
  return tmp;

} // grid<T>::operator*

// __________________________________________________________________________
// Determine the LU decomposition of a square grid.

template<class T>
grid<T> grid<T>::LU() const
{
  assert(invariant());
  assert(nr == nc);

  grid<T> tmp = *this;
  int i, j;
  for (i = 0; i < nr - 1; ++i)
  {
    for (j = i + 1; j < nr; ++j)
      tmp(j, i) /= tmp(i, i);
    for (j = i + 1; j < nr; ++j)
      for (int k = i + 1; k < nr; ++k)
	tmp(j, k) -= tmp(j, i) * tmp(i, k);
  }

  assert(invariant() && tmp.invariant());
  return tmp;

} // grid::LU

// __________________________________________________________________________
// Calculate the matrix inverse of a grid.

template<class T>
grid<T> grid<T>::inverse() const
{
  assert(invariant());
  assert(nr == nc);

  grid<T> tmp = *this;

  grid<int> p(nr);
  int i, j, k;
  for (j = 0; j < nr; ++j)
    p(j) = j;
  grid<double> hv(nr);
  hv.clear();

  for (j = 0; j < nr; ++j)
  {
    T max = fabs(tmp(j, j));
    int r = j;
    for (i = j + 1; i < nr; ++i) {
      if (fabs(tmp(i, j)) > max) {
	max = fabs(tmp(i, j));
	r = i;
      }
    }
    if (max == 0.0) {
      std::cerr << "Unable to invert a matrix" << std::endl;
      exit(1);
    }
    if (r > j) {
      for (k = 0; k < nr; ++k)
	std::swap(tmp(j, k), tmp(r, k));
      std::swap(p(j), p(r));
    }
    T hr = 1 / tmp(j, j);
    for (i = 0; i < nr; ++i)
      tmp(i, j) *= hr;
    tmp(j, j) = hr;
    for (k = 0; k < nr; ++k) {
      if (k == j)
	continue;
      for (i = 0; i < nr; ++i) {
	if (i == j)
	  continue;
	tmp(i, k) -= tmp(i, j) * tmp(j, k);
      }
      tmp(j, k) *= (-hr);
    }
  }
  for (i = 0; i < nr; ++i) {
    for (k = 0; k < nr; ++k)
      hv(p(k)) = tmp(i, k);
    for (k = 0; k < nr; ++k)
      tmp(i, k) = hv(k);
  }

  assert(invariant() && tmp.invariant());
  return tmp;

} // grid::inverse

// __________________________________________________________________________

template<class T>
const grid<T> grid<T>::transpose() const
{
  assert(invariant());
  grid<T> tmp(nc, nr);
  for (int i = 0; i < nc; ++i)
    for (int j = 0; j < nr; ++j)
      tmp(i, j) = (*this)(j, i);
  assert(invariant() && tmp.invariant());

  return tmp;

} // grid::transpose

// __________________________________________________________________________

template<class T>
void grid<T>::sort(int col, int left, int right)
{
  if (left == -1 && right == -1) {
    left = 0;
    right = nr - 1;
  }

  int i = left;
  int j = right;

  T midval = (*this)((left + right) / 2, col);

  do {
    while ((*this)(i, col) < midval && i < right)
      ++i;
    while (midval < (*this)(j, col) && j > left)
      --j;

    if (i <= j) 
    {
      for (int k = 0; k < nc; ++k) {
	T tmpval = (*this)(i, k);
	(*this)(i, k) = (*this)(j, k);
	(*this)(j, k) = tmpval;
      }
      ++i;
      --j;
    }

  } 
  while (i <= j);

  if (left < j)
    sort(col, left, j);
  if (i < right)
    sort(col, i, right);

} // grid::sort

// __________________________________________________________________________

template <class T>
int grid<T>::loadpgm(const char* pgmname)
{
  // Open the pgm file
  Ifstream(ifs, pgmname);
  if (!ifs)
    return 0;

  // Find the first line that doesn't begin with white space.
  char buf[255];
  char pchar;
  int mode, r, c, maxval, matches;
  do {
    ifs.Getline(buf, 255);
  } while (!ifs.eof() && (buf[0] == 0 || buf[0] == ' '));
  if (ifs.eof())
    return 0;

  // Make sure the file is a pgm file.  Determine the pgm mode, the
  // dimensions, and the range of pixel values
  matches = sscanf(buf, "%c%1d%d%d%d", &pchar, &mode, &r, &c, &matches);
  if (matches < 2 || mode < 2 || mode > 6)
    return 0;

  if (matches < 5)
  {
    ifs.Getline(buf, 255);
    while (!ifs.eof() && (buf[0] == '#' || buf[0] == 0))
    ifs.Getline(buf, 255);
    sscanf(buf, "%d %d", &r, &c);
    ifs.Getline(buf, 255);
    sscanf(buf, "%d", &maxval);
    if (ifs.eof() || r <= 0 || c <= 0 || maxval <= 0)
      return 0;
  }

  init(r, c);
  int i, j;
  switch(mode)
  {
   case 2:
    for(j = 0; j < c; ++j) {
      int tmpi;
      for(i = 0; i < r; ++i) {
	ifs >> tmpi;
	(*this)(i, j) = tmpi;
      }
    };
    break;

   case 3:
    for(j = 0; j < c; j++) {
      int tmpi[3];
      for(i = 0; i < r; ++i) {
	ifs >> tmpi[0] >> tmpi[1] >> tmpi[2];
	(*this)(i, j) = uchar(0.212671 * tmpi[0] + 0.715160 * tmpi[1] + 
                              0.072169 * tmpi[2]);
      }
    };
    break;

   case 5:
    for(j = 0; j < c; ++j) {
      for(i = 0; i < r; ++i)
	ifs.read((char *)&((*this)(i, j)), 1);
    };
    break;

   case 6:
    for(j=0; j < c; j++) {
      uchar tmpc[3];
      for(i = 0; i < r; ++i) {
        ifs.read((char *)tmpc, 3);
	(*this)(i, j) = uchar(0.212671 * tmpc[0] + 0.715160 * tmpc[1] + 
                              0.072169 * tmpc[2]);
      }
    };
    break;

   default:
    return 0;
    break;

  };

  return 1;

} // grid::loadpgm

// __________________________________________________________________________

template <class T>
int grid<T>::savepgm(const char* pgmname)
{
  // Open the pgm file
  Ofstream(ofs, pgmname);
  if (!ofs)
    return 0;

  // Print the pgm header
  ofs << "P5\n" << nr << " " << nc << "\n255\n";

  for(int j = 0; j < nc; ++j) {
    for(int i = 0; i < nr; ++i)
      ofs.write((char *)&((*this)(i, j)), 1);
  };

  return 1;

} // grid::savepgm

// __________________________________________________________________________
//
template<class T>
const grid<T> operator+(const grid<T>& m, const grid<T>& n) 
{ grid<T> p = m; return (p += n); }

template<class T>
const grid<T> operator-(const grid<T>& m, const grid<T>& n) 
{ grid<T> p = m; return (p -= n); }

typedef grid<char> cgrid;
typedef grid<unsigned char> ucgrid;
typedef grid<int> igrid;
typedef grid<unsigned int> uigrid;
typedef grid<long> lgrid;
typedef grid<float> fgrid;
typedef grid<double> dgrid;
typedef std::vector<cgrid> cgrids;
typedef std::vector<ucgrid> ucgrids;
typedef std::vector<ucgrids> ucgridss;
typedef std::vector<igrid> igrids;
typedef std::vector<uigrid> uigrids;
typedef std::vector<lgrid> lgrids;
typedef std::vector<fgrid> fgrids;
typedef std::vector<dgrid> dgrids;

// grid.h

#endif // GRID_H

