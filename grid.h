#ifndef GRID_H_
#define GRID_H_

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "utils.h"

namespace grid_h {

typedef std::complex<double> cplx;

enum { gengrid = 0, posgrid = 1, symgrid = 2, trigrid = 3};

template <class T>
class grid {
 private:
  size_t nr, nc;  // number of rows, columns
  std::vector<T> sto;  // column major grid storage

  // Consistency Checking Functions
  int consistent() const { return sto.size() == nr * nc; }
  int inrange(size_t r, size_t c) const { return r < nr && c < nc; }

 public:
  // Constructors, Operator=
  grid() : nr(0), nc(0) { }
  explicit grid(size_t r) : nr(r), nc(1), sto(r) { }
  grid(size_t r, size_t c) : nr(r), nc(c), sto(r * c) { }
  grid(size_t r, size_t c, const std::vector<T>& v) : nr(r), nc(c), sto(v) {
    assert(r * c == v.size()); }
  grid(const grid &m) : nr(m.rows()), nc(m.cols()), sto(m.storage()) { }
  explicit grid(const std::vector<T>& v) : nr(v.size()), nc(1), sto(v) { }
  grid(const std::string& s, size_t N);  // Special Matrices
  grid(const std::string& s, T start, T end, T step);  // Sequences

  grid& operator=(const grid &g);
  template <class S> void operator=(const grid<S>& g);
  bool operator==(const grid& g) const {
    if (this == &g) {
      return true;
    }
    return nr == g.nr && nc == g.nc && sto == g.sto;
  }
  bool operator==(const T& val) const {  // TODO(wilder): Write tests!
    return !std::any_of(
        begin(), end(), [val](const T& el){ return el != val; });
  }
  void subgrid(
      grid* m, size_t r, size_t c, size_t numrows, size_t numcols) const;
  grid& rbind(const grid& g);
  grid& cbind(const grid& g);
  typename std::vector<T>::iterator begin() { return sto.begin(); }
  typename std::vector<T>::iterator end() { return sto.end(); }
  typename std::vector<T>::const_iterator begin() const { return sto.begin(); }
  typename std::vector<T>::const_iterator end() const { return sto.end(); }

  // Basic Member Access Functions
  size_t rows() const { return nr; }
  size_t cols() const { return nc; }
  const std::vector<T>& storage() const { return sto; }

  const T& operator()(size_t r) const {
    assert(inrange(r, 0));
    return sto[r];
  }
  T& operator()(size_t r) {
    assert(inrange(r, 0));
    return sto[r];
  }
  const T& operator()(size_t r, size_t c) const {
    assert(inrange(r, c));
    return sto[c * nr + r];
  }
  T& operator()(size_t r, size_t c) {
    assert(inrange(r, c));
    return sto[c * nr + r];
  }
  T& operator[](size_t r) { return sto[r]; }
  const T& operator[](size_t r) const { return sto[r]; }
  void set(size_t r, size_t c, T val) {
    if (inrange(r, c)) {
      (*this)(r, c) = val;
    }
  }
  T get(size_t r, size_t c) const { return inrange(r, c) ? (*this)(r, c) : 0; }

  void resize(size_t r, size_t c) { nr = r; nc = c; sto.resize(r * c); }
  void resize(size_t r) { resize(r, 1); }
  void resize() { resize(0, 0); }
  grid& fill(const T& val) { std::fill(begin(), end(), val); return *this; }
  grid& clear() { return fill(0); }

  // I/O Functions
  void write(const char *file);
  int write(std::ofstream& os);
  int read(const char* file);
  int read(std::ifstream& is);
  int loadpgm(const std::string& pgmname);
  int savepgm(const std::string& pgmname);
  size_t memory_size() const;  // TODO(wilder): Remove this?
  friend std::ostream& operator<<(std::ostream& os, const grid& g) {
    os << "Dim (" << g.nr << ", " << g.nc << ") Storage {" << g.storage()[0];
    for (auto el = g.storage().begin() + 1; el != g.storage().end(); ++el) {
      os << ", " << *el;
    }
    os << "}\n";
    return os;
  }
  void dump(size_t max = 0, bool invert = false) const;

  // Useful Utility Functions
  grid& operator<<(grid &m);
  grid& operator+=(const grid &m);
  grid& operator-=(const grid &m);
  grid& operator*=(const grid &m);
  grid& operator/=(const grid &m);
  grid operator+(const grid &m);
  grid operator-(const grid &m);
  grid operator*(const grid &m);
  grid operator/(const grid &m);
  grid& operator+=(T val);
  grid& operator-=(T val) { return *this += -val;  }
  grid operator+ (T val) const;
  grid operator*(const grid &m) const;
  grid& scale(T val);
  grid& transform(T minval, T maxval);
  grid transpose() const;
  grid LU() const;
  grid inverse() const;
  size_t offpixels() {
    size_t c = 0;
    for (auto& el : sto) {
      c += el == 0;
    }
    return c;
  }
  int onpixels() { return(sto.size() - offpixels()); }
  void sort_rows(size_t col);
};

// __________________________________________________________________________
// Dump a grid - for classes with ostream<<(const T&).

template <class T>
void grid<T>::dump(size_t max, bool invert) const {
  if (invert) {
    if (max == 0) {
      max = nr;
    }
    max = std::min(nr, max);
    for (size_t i = 0; i < max; ++i) {
      for (size_t j = 0; j < nc; ++j) {
        std::cout << ((*this)(i, j)) << " ";
      }
      std::cout << "\n";
    }
  } else {
    if (max == 0) {
      max = nc;
    }
    max = std::min(nc, max);
    for (size_t j = 0; j < max; ++j) {
      for (size_t i = 0; i < nr; ++i) {
        std::cout << ((*this)(i, j)) << " ";
      }
      std::cout << "\n";
    }
  }
}

// __________________________________________________________________________
// Write out a grid to a file

template<class T>
void grid<T>::write(const char *file) {
  assert(consistent());
  Ofstream(ofs, file);
  if (!ofs) {
    std::cerr << "Unable to write a grid to file [" << file << "]."
      << std::endl;
    return;
  }

  ofs.write("GR11", 4);
  varwrite(ofs, nr);
  varwrite(ofs, nc);
  arraywrite(ofs, sto, sto.size());
}

// __________________________________________________________________________
// Write out a grid to an open output stream

template<class T>
int grid<T>::write(std::ofstream& ofs) {
  if (!ofs) {
    return 0;
  }
  varwrite(ofs, nr);
  varwrite(ofs, nc);
  arraywrite(ofs, sto, sto.size());
  return 1;
}

// __________________________________________________________________________
// Calculate the storage size of a grid to be written out

template<class T>
size_t grid<T>::memory_size() const {
  return sizeof(nr) + sizeof(nc) + sto.size() * sizeof(T);
}

// __________________________________________________________________________
// Read in a grid which was previously written to a file

template<class T>
int grid<T>::read(const char *file) {
  assert(consistent());
  Ifstream(ifs, file);
  if (!ifs)
    return 0;

  char version[4];
  ifs.read(version, 4);
  if (memcmp(version, "GR11", 4) != 0 && memcmp(version, "GR12", 4) != 0) {
    if (loadpgm(file)) {
      return 1;
    } else {
      std::cerr << "The file [" << file << "] is not a grid or pgm file"
        << std::endl;
      return 0;
    }
  }

  resize(0, 0);
  if (!memcmp(version, "GR11", 4)) {
    varread(ifs, nr);
    varread(ifs, nc);
    if (nr > 0 && nc > 0) {
      resize(nr * nc);
      arrayread(ifs, sto, nr * nc);
    }
  } else if (!memcmp(version, "GR12", 4)) {
    ifs >> nr;
    ifs >> nc;
    assert(nr >= 0 && nc >= 0);
    if (nr > 0 && nc > 0) {
      resize(nr * nc);
      for (size_t j = 0; j < nc; ++j) {
        for (size_t i = 0; i < nr; ++i) {
          ifs >> (*this)(i, j);
        }
      }
    }
  }
  assert(consistent());

  return 1;
}

// __________________________________________________________________________
// Read in a grid from an open istream

template<class T>
int grid<T>::read(std::ifstream& is) {
  assert(consistent());

  // if (!is.is_open() || is.eof())
  if (!is || is.eof())
    return 0;

  resize(0, 0);
  varread(is, nr);
  varread(is, nc);
  assert(nr >= 0 && nc >= 0);
  if (nr > 0 && nc > 0) {
    resize(nr * nc);
    arrayread(is, sto, nr * nc);
  }
  assert(consistent());
  return 1;
}

// ==========================================================================
// Special matrix constructors

template<class T>
grid<T>::grid(const std::string& s, size_t n) : nr(n), nc(n), sto(n * n, 0) {
  if (s == "I" || s.substr(0, 3) == "ide") {  // Identity Matrix
    for (size_t i = 0; i != n; ++i)
      (*this)(i, i) = 1;
  } else if (s.substr(0, 3) == "hil") {  // Badly conditioned Hilbert matrix
    for (size_t i = 0; i != n; ++i)
      for (size_t j = 0; j != n; ++j)
        (*this)(i, j) = T(1.0 / (i + j + 1.0));
  } else if (s.substr(0, 3) == "wil") {  // Wilkinson eigenvalue matrix
    for (size_t i = 0; i != n - 1; ++i) {
      (*this)(i, i + 1) = 1;
      (*this)(i + 1, i) = 1;
    }
    for (size_t i = 0; i != n; ++i)
      (*this)(i, i) = T(std::abs(i - (n - 1) / 2.0));
  } else {
    // throw grid_error("Unknown special matrix initialization");
  }
}

template<class T>
grid<T>::grid(const std::string& s, T start, T end, T step) {
  if (s.substr(0, 3) == "seq") {
    size_t count = std::abs((end - start) / step) + 1;
    resize(count);
    for (size_t i = 0; i != count; ++i) {
      sto[i] = start;
      start += step;
    }
  } else {
    // throw grid_error("grid(string&, T, T, T): Unknown initialization");
  }
}

// __________________________________________________________________________
// Set a grid equal to a grid 'g'; leave 'g' unchanged.

template<class T>
grid<T>& grid<T>::operator=(const grid& g) {
  if (this != &g) {
    nr = g.nr;
    nc = g.nc;
    sto = g.storage();
  }
  return *this;
}

template<class T>
template<class S>
void grid<T>::operator=(const grid<S>& g) {
  assert(this != &g);
  nr = g.nr;
  nc = g.nc;
  std::copy(g.begin(), g.end(), begin());
}

// __________________________________________________________________________
// Set '*m' equal to a subgrid of *this.

template<class T>
void grid<T>::subgrid(
    grid* m, size_t r, size_t c, size_t numrows, size_t numcols) const {
  assert(r >= 0 && c >= 0 && numrows >= 0 && numcols >= 0);
  assert(r + numrows <= nr && c + numcols <= nc);

  if (this != m) {
    m->resize(numrows, numcols);
    if (nr > 0 && nc > 0) {
      for (size_t i = 0; i < m->nr; ++i) {
        for (size_t j = 0; j < m->nc; ++j) {
          (*m)(i, j) = (*this)(r + i, c + j);
        }
      }
    }
  } else {
    grid tmp;
    subgrid(&tmp, r, c, numrows, numcols);
    *m << tmp;
  }
}

// ==========================================================================
// Append the rows or columns of 'g' to '*this'. Works when &g == this.

template <class T>
grid<T>& grid<T>::cbind(const grid<T>& g) {
  assert(nr == g.rows());
  size_t old_size = sto.size();
  size_t old_g_size = g.storage().size();
  resize(nr, nc + g.cols());
  std::copy(g.begin(), g.begin() + old_g_size, begin() + old_size);
  return *this;
}

template<class T>
grid<T>& grid<T>::rbind(const grid<T>& g) {
  assert(nc == g.cols());
  grid<T> combined(nr + g.rows(), nc);
  auto this_pos = begin();
  auto g_pos = g.begin();
  auto combined_pos = combined.begin();
  for (size_t i = 0; i != nc; ++i) {
    std::copy(this_pos, this_pos + nr, combined_pos);
    std::copy(g_pos, g_pos + g.rows(), combined_pos + nr);
    this_pos += nr;
    g_pos += g.rows();
    combined_pos += nr + g.rows();
  }
  return *this << combined;
}

// __________________________________________________________________________
// Set a grid equal to 'm'; obliterate 'm'.

template<class T>
grid<T>& grid<T>::operator<<(grid &m) {
  if (this != &m) {
    resize(0, 0);
    std::swap(*this, m);
  }
  return *this;
}

// __________________________________________________________________________
// Add another grid to the current grid.

template<class T>
grid<T>& grid<T>::operator+=(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  std::transform(begin(), end(), m.sto.begin(), begin(), std::plus<T>());
  return *this;
}

template<class T>
grid<T> grid<T>::operator+(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  return grid<T>(*this) += m;
}

// __________________________________________________________________________
// Subtract a grid from the current grid.

template<class T>
grid<T>& grid<T>::operator-=(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  std::transform(begin(), end(), m.sto.begin(), begin(), std::minus<T>());
  return *this;
}

template<class T>
grid<T> grid<T>::operator-(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  return grid<T>(*this) -= m;
}

// __________________________________________________________________________
// Multiply a grid to the current grid.

template<class T>
grid<T>& grid<T>::operator*=(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  std::transform(begin(), end(), m.sto.begin(), begin(), std::multiplies<T>());
  return *this;
}

template<class T>
grid<T> grid<T>::operator*(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  return grid<T>(*this) *= m;
}

// __________________________________________________________________________
// Divide the current grid by another grid.

template<class T>
grid<T>& grid<T>::operator/=(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  std::transform(begin(), end(), m.sto.begin(), begin(), std::divides<T>());
  return *this;
}

template<class T>
grid<T> grid<T>::operator/(const grid& m) {
  assert(nr == m.nr && nc == m.nc);
  return grid<T>(*this) /= m;
}

// __________________________________________________________________________
// Add a fixed value to each element of a grid

template<class T>
grid<T>& grid<T>::operator+=(T val) {
  for (auto& el : sto) { el += val; }
  return *this;
}

// __________________________________________________________________________
// Add a fixed value to each element of a grid

template<class T>
grid<T> grid<T>::operator+(T val) const {
  return grid<T>(*this) += val;
}

// __________________________________________________________________________
// Scale the values of a grid so that the largest magnitude is 'val'.

template<class T>
grid<T>& grid<T>::scale(T val) {
  if (nr > 0 && nc > 0) {
    const T gmin = *std::min_element(begin(), end());
    const T gmax = *std::max_element(begin(), end());
    const T absmax = std::max(std::abs(gmin), std::abs(gmax));
    if (absmax != 0) {
      for (auto& el : *this) { el = el * val / absmax; }
    }
  }
  return *this;
}

// __________________________________________________________________________
// Linearly transform the values of a grid so that the values range
// from 'val1' to 'val2'.

template<class T>
grid<T>& grid<T>::transform(T val1, T val2) {
  if (nr == 0 || nc == 0)
    return *this;
  const T gmin = *std::min_element(begin(), end());
  const T gmax = *std::max_element(begin(), end());
  const T oldrange = gmax - gmin;
  const T newrange = val2 - val1;
  if (oldrange != 0) {
    for (auto& el : *this) {
      el = (el - gmin) * newrange / oldrange + val1;
    }
  } else if (gmin < val1) {
    fill(val1);
  } else if (gmax > val2) {
    fill(val2);
  }
  return *this;
}

// __________________________________________________________________________
// Multiply two grids together using matrix multiplication.

/*
template<class T>
grid<T> grid<T>::operator*(const grid &m) const {
  grid<T> tmp(nr, m.nc);
  assert(nc == m.nr);
  tmp.clear();
  for (size_t i = 0; i < nr; ++i)
    for (size_t j = 0; j < m.nc; ++j)
      for (size_t k = 0; k < nc; ++k)
        tmp(i, j) += (*this)(i, k) * m(k, j);
  return tmp;
}
*/

// __________________________________________________________________________
// Determine the LU decomposition of a square grid.

template<class T>
grid<T> grid<T>::LU() const {
  grid<T> tmp = *this;
  assert(nr == nc);
  if (nr > 0) {
    for (size_t i = 0; i < nr - 1; ++i) {
      for (size_t j = i + 1; j < nr; ++j)
        tmp(j, i) /= tmp(i, i);
      for (int j = i + 1; j < nr; ++j)
        for (int k = i + 1; k < nr; ++k)
          tmp(j, k) -= tmp(j, i) * tmp(i, k);
    }
  }
  return tmp;
}

// __________________________________________________________________________
// Calculate the matrix inverse of a grid.

template<class T>
grid<T> grid<T>::inverse() const {
  grid<T> tmp = *this;
  assert(nr == nc);
  grid<int> p(nr);
  int i, j, k;
  for (j = 0; j < nr; ++j)
    p(j) = j;
  grid<double> hv(nr);
  hv.clear();

  for (j = 0; j < nr; ++j) {
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
  return tmp;
}

// __________________________________________________________________________

template<class T>
grid<T> grid<T>::transpose() const {
  grid<T> transpose(nc, nr);
  auto this_pos = begin();
  for (size_t i = 0; i < nc; ++i)
    for (size_t j = 0; j < nr; ++j)
      transpose(i, j) = *(this_pos++);
  return transpose;
}

// __________________________________________________________________________
// Sort the rows of a grid according to the values in a specified column.

template<class T>
void grid<T>::sort_rows(size_t col) {
  if (sto.size() == 0) {
    return;
  }
  const auto col_start = begin() + col * rows();
  std::vector<size_t> permutation(rows());
  std::iota(permutation.begin(), permutation.end(), 0);
  std::sort(permutation.begin(), permutation.end(),
      [&col_start](size_t i, size_t j) {
        return *(col_start + i) < *(col_start + j); });
  std::vector<T> tmp(rows());
  auto pos = begin();
  for (size_t j = 0; j < cols(); ++j) {
    std::copy(pos, pos + rows(), tmp.begin());
    for (auto ip = permutation.begin(); ip != permutation.end(); ++ip) {
      *(pos++) = tmp[*ip];
    }
  }
}

// __________________________________________________________________________

template <class T>
int grid<T>::loadpgm(const std::string& pgmname) {
  Ifstream(ifs, pgmname);
  if (!ifs) {
    return 0;
  }

  // Find the first line that doesn't begin with white space.
  char buf[255];
  char pchar;
  int mode, r, c, maxval, matches;
  do {
    ifs.Getline(buf, 255);
  } while (!ifs.eof() && (buf[0] == 0 || buf[0] == ' '));
  if (ifs.eof()) {
    return 0;
  }

  // Make sure the file is a pgm file.  Determine the pgm mode, the
  // dimensions, and the range of pixel values
  matches = sscanf(buf, "%c%1d%d%d%d", &pchar, &mode, &r, &c, &matches);
  if (matches < 2 || mode < 2 || mode > 6) {
    return 0;
  }

  if (matches < 5) {
    ifs.Getline(buf, 255);
    while (!ifs.eof() && (buf[0] == '#' || buf[0] == 0)) {
      ifs.Getline(buf, 255);
    }
    sscanf(buf, "%d %d", &r, &c);
    ifs.Getline(buf, 255);
    sscanf(buf, "%d", &maxval);
    if (ifs.eof() || r <= 0 || c <= 0 || maxval <= 0) {
      return 0;
    }
  }

  resize(r, c);
  int i, j;
  switch (mode) {
    case 2:
      for (j = 0; j < c; ++j) {
        int tmpi;
        for (i = 0; i < r; ++i) {
          ifs >> tmpi;
          (*this)(i, j) = tmpi;
        }
      }
      break;

    case 3:
      for (j = 0; j < c; j++) {
        int tmpi[3];
        for (i = 0; i < r; ++i) {
          ifs >> tmpi[0] >> tmpi[1] >> tmpi[2];
          (*this)(i, j) = static_cast<unsigned char>(
              0.212671 * tmpi[0] + 0.715160 * tmpi[1] + 0.072169 * tmpi[2]);
        }
      }
    break;

    case 5:
      char tmp;
      for (auto& el : sto) {
        ifs.read(&tmp, 1);
        el = c;
      }
      break;

    case 6:
      for (j=0; j < c; j++) {
        std::vector<char> tmpc(3);
        for (i = 0; i < r; ++i) {
          ifs.read(&tmpc[0], 3);
          (*this)(i, j) = static_cast<unsigned char>(
              0.212671 * static_cast<unsigned char>(tmpc[0]) +
              0.715160 * static_cast<unsigned char>(tmpc[1]) +
              0.072169 * static_cast<unsigned char>(tmpc[2]));
        }
      }
      break;

    default:
      return 0;
      break;
  }
  return 1;
}

// __________________________________________________________________________

template <class T>
int grid<T>::savepgm(const std::string& pgmname) {
  Ofstream(ofs, pgmname);
  if (!ofs) {
    return 0;
  }
  ofs << "P5\n" << nr << " " << nc << "\n255\n";
  for (const auto& el : sto) {
    char c = static_cast<char>(el);
    ofs.write(&c, 1);
  }
  return 1;
}

// __________________________________________________________________________

template<class T>
const grid<T> operator+(const grid<T>& m, const grid<T>& n) {
  auto p = m; return p += n;
}

template<class T>
const grid<T> operator-(const grid<T>& m, const grid<T>& n) {
  auto p = m; return p -= n;
}

};  // namespace grid_h

#endif  // GRID_H_

/* 


// ==========================================================================
// Grid (matrix) of elements of type T.
// A grid is essentially a standard library vector with added functionality.

template<class T>
class grid {
 protected:
  vector<T> sto;
  uint nr, nc; // number of rows and columns

 public:
  typedef typename vector<T>::iterator iterator;
  typedef typename vector<T>::const_iterator const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;

  // Constructors, Destructor, operator= and operator==
  grid() : nr(0), nc(0) { }
  explicit grid(uint r, uint c=1) : sto(r * c), nr(r), nc(c) { }
  grid(uint r, uint c, T v) : sto(r * c, v), nr(r), nc(c) { }
  explicit grid(const string& s) : nr(0), nc(0) { read(s); } // File input
  grid(const T* tp, uint r, uint c) : sto(r * c), nr(r), nc(c)
    { std::copy(tp, tp + r * c, begin()); }
  grid(const vector<T>& v) : sto(v), nr(v.size()), nc(1) { }
  grid(const grid<T>& g) : sto(g.sto), nr(g.nr), nc(g.nc) { }
  ~grid() { }
  grid<T>& operator=(const grid<T> &g);
  template<class S> void operator=(const grid<S> &g);
  bool operator==(const grid<T> &g) const;
  operator const vector<T>&() const { return sto; }
  // The following non-const conversion can be dangerous.
  // operator vector<T>&() { return sto; }
  vector<T>& vec() { return sto; }

  // Basic Member Functions
  void clear() { nr = nc = 0; sto.clear(); }
  uint size() const { return (uint) sto.size(); }
  T& front() { return sto.front(); }
  const T& front() const { return sto.front(); }
  iterator begin() { return sto.begin(); }
  iterator end() { return sto.end(); }
  const_iterator begin() const { return sto.begin(); }
  const_iterator end() const { return sto.end(); }
  void resize(uint r, uint c = 1) { resize(r * c); nr = r; nc = c; }
  void resize(uint r, uint c, const T& v)
    { resize(r * c, v); nr = r; nc = c; }
  void fill(const T& val = 0) { std::fill(begin(), end(), val); }
  void push_back(const T& t)
    { if (nc == 1) ++nr; else if (nr == 1) ++nc; else return;
      sto.push_back(t); }

  // I/O Functions
  bool write(const string& file, bool header = 1, bool binary = 1);
  bool write(std::ofstream& os, bool binary = 1);
  bool read(const string& file, bool header = 1);
  bool read(std::ifstream& is, bool binary = 1);
  int memory_size() const;

  // Manipulate the contents of grids.
  void swap(vector<T>& v)
    { nr = v.size(); nc = 1; sto.swap(v); }
  void swap(grid<T>& g) {
    if (this != &g)
      { sto.swap(g.sto); std::swap(nr, g.nr); std::swap(nc, g.nc); }
  }
  grid<T>& operator<<(vector<T> &v)
    { swap(v); v.clear(); return *this; }
  grid<T>& operator<<(grid<T> &g) {
    if (this != &g) { swap(g); g.sto.clear(); g.nr = g.nc = 0; }
    return *this;
  }

  void subgrid(grid<T>& g, uint r, uint c, uint numrows, uint numcols);

  // Useful Utility Functions
  void normalize(T val1, T val2);
  grid<T>& sort() { std::sort(begin(), end()); return *this; }

  // Exceptions
  class grid_error : public std::exception {
   private:
    string msg;
   public:
    virtual const char* what() const throw() { return(msg.c_str()); }
    grid_error(const string& s = "") : msg(s) { }
    ~grid_error() throw() { }
  };

}; // class grid<T>

// ==========================================================================
// Useful typedefs

typedef vector<bool> bvector;
typedef vector<char> cvector;
typedef vector<unsigned char> ucvector;
typedef vector<int> ivector;
typedef vector<uint> uivector;
typedef vector<long> lvector;
typedef vector<float> fvector;
typedef vector<double> dvector;
typedef vector<cplx> zvector;
typedef vector<string> strvector;

typedef bvector::iterator biter;
typedef cvector::iterator citer;
typedef ucvector::iterator uciter;
typedef ivector::iterator iiter;
typedef uivector::iterator uiiter;
typedef lvector::iterator liter;
typedef fvector::iterator fiter;
typedef dvector::iterator diter;
typedef zvector::iterator ziter;
typedef strvector::iterator striter;

typedef bvector::const_iterator bconstiter;
typedef cvector::const_iterator cconstiter;
typedef ucvector::const_iterator ucconstiter;
typedef ivector::const_iterator iconstiter;
typedef uivector::const_iterator uiconstiter;
typedef lvector::const_iterator lconstiter;
typedef fvector::const_iterator fconstiter;
typedef dvector::const_iterator dconstiter;
typedef zvector::const_iterator zconstiter;
typedef strvector::const_iterator strconstiter;

typedef grid<bool> bgrid;
typedef grid<char> cgrid;
typedef grid<unsigned char> ucgrid;
typedef grid<int> igrid;
typedef grid<uint> uigrid;
typedef grid<long> lgrid;
typedef grid<float> fgrid;
typedef grid<double> dgrid;
typedef grid<cplx> zgrid;
typedef grid<string> strgrid;

typedef vector<bgrid> bgrids;
typedef vector<cgrid> cgrids;
typedef vector<ucgrid> ucgrids;
typedef vector<igrid> igrids;
typedef vector<uigrid> uigrids;
typedef vector<lgrid> lgrids;
typedef vector<fgrid> fgrids;
typedef vector<dgrid> dgrids;
typedef vector<zgrid> zgrids;

typedef vector<bgrids> bgridss;
typedef vector<cgrids> cgridss;
typedef vector<ucgrids> ucgridss;
typedef vector<igrids> igridss;
typedef vector<uigrids> uigridss;
typedef vector<lgrids> lgridss;
typedef vector<fgrids> fgridss;
typedef vector<dgrids> dgridss;
typedef vector<zgrids> zgridss;

//============================================================================

template<> inline grid<cplx>& grid<cplx>::sort() { return *this; }

// ==========================================================================
// Print a grid; Will only work for classes with ostream<<(const T&).
// If 'tr' is non-zero, the output will be transposed, and if 'max'
// is non-zero, at most 'max' rows will be output.  This function
// provides somewhat finer control over the output than operator<<.
// By default, fields are separated by spaces.  An alternative separator
// can be specified, but this is designed for compatibility with other
// programs as then the file can not be easily read back using this code.

template<class T>
void dump(const grid<T>& g, std::ostream& os=std::cout,
          bool tr=false, size_t max=0, uint prec=4, char sep=' ')
{
  if (g.size() == 0) { os << "<Empty grid>\n"; return; }
  int wid = prec + 6;
  if (tr) {
    if (max == 0 || max > g.cols()) max = g.cols();
    for (size_t j = 0; j != max; ++j) {
      for (size_t i = 0; i != g.rows(); ++i)
        os << std::setw(wid) << std::setprecision(prec)
                  << (g(i, j)) << sep;
      os << "\n";
    }
  } else {
    if (max == 0 || max > g.rows()) max = g.rows();
    for (size_t i = 0; i != max; ++i) {
      for (size_t j = 0; j != g.cols(); ++j)
        os << std::setw(wid) << std::setprecision(prec)
                  << (g(i, j)) << sep;
      os << "\n";
    }
  }
}

// Call the above grid<T>::dump with a fixed field output.

template<class T>
void dumpfixed(const grid<T>& g, std::ostream& os=std::cout,
               uint tr=0, uint max=0, int prec=4, char sep=' ')
{

  std::ios::fmtflags f = os.setf(std::ios::fixed);
  dump(g, os, tr, max, prec, sep);
  os.flags(f);
}

// Print a vector/grid to an ostream using a default format similar to
// the one used in matlab.  A newline is printed before the vector/grid.

template<class T>
std::ostream& operator<<(std::ostream& os, const vector<T>& v)
{
  if (v.size() == 0) { os << "<Empty vector>\n"; return os; }
  const int wid = 10, prec = 4;
  std::ios::fmtflags f = os.setf(std::ios::fixed);
  for (uint i = 0; i != v.size(); ++i)
    os << std::setw(wid) << std::setprecision(prec) << v[i] << " ";
  os << "\n";
  os.flags(f);
  return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const grid<T>& g)
  { dumpfixed(g, os, 0, 0, 4); return os; }

// ==========================================================================
// operator<< specialized to unsigned char, defined in grid.C

template<> std::ostream& operator<<(std::ostream& os, const ucgrid& g);

// ==========================================================================
// (Binary) write out a grid to an open output stream.

template<class T>
bool grid<T>::write(std::ofstream& ofs, bool binary)
{
  if (!ofs)
    return false;

  if (binary) {
    ofs.write((char*) &nr, sizeof(nr));
    ofs.write((char*) &nc, sizeof(nc));
    for (uint i = 0; i != sto.size(); ++i)
      ofs.write((char*) &sto[i], sizeof(T));
  } else {
    ofs << nr << " " << nc << "\n";
    for (uint i = 0; i != rows(); ++i) {
      for (uint j = 0; j != cols(); ++j)
        ofs << grid<T>::operator()(i, j);
      ofs << "\n";
    }
  }
  return true;

} // grid::write

// ==========================================================================
// (Binary) write out a grid to a file.

template<class T>
bool grid<T>::write(const string& file, bool header, bool binary)
{
  std::ofstream ofs(file.c_str());
  if (!ofs) {
    return false;
  }

  if (header) {
    ofs.write("GR11", 4);
  }
  return write(ofs, binary);

} // grid::write

// ==========================================================================
// Read in a grid which was previously written to a file

template<class T> bool grid<T>::read(const string& file, bool header)
{
  std::ifstream ifs(file.c_str());
  if (!ifs)
    return false;

  char version[4];
  if (header == 0)
    memcpy(version, "GR12", 4);
  else
    ifs.read(version, 4);
  if (memcmp(version, "GR11", 4) != 0 && memcmp(version, "GR12", 4) != 0) {
    if (loadpgm(*this, file))
      return true;
    else {
      std::cerr << "[" << file << "] is not a grid or pgm file" << std::endl;
      return false;
    }
  }

  // GR11 indicates a platform dependent binary file
  if (!memcmp(version, "GR11", 4))
    read(ifs, 1);
  // GR12 should be a platform independent text file
  else if (!memcmp(version, "GR12", 4))
    read(ifs, 0);
  return true;

} // grid::read

// ==========================================================================
// Read in a grid from an open istream

template<class T>
bool grid<T>::read(std::ifstream& ifs, bool binary)
{
  if (!ifs || ifs.eof()) return false;

  if (binary) {
    ifs.read((char*) &nr, sizeof(nr));
    ifs.read((char*) &nc, sizeof(nc));
    resize(nr, nc);
    for (uint i = 0; i != nr * nc; ++i)
      ifs.read((char*) &sto[i], sizeof(T));
  } else {
    ifs >> nr;
    ifs >> nc;
    resize(nr, nc);
    for (uint i = 0; i != rows(); ++i)
      for (uint j = 0; j != cols(); ++j)
        ifs >> grid<T>::operator()(i, j);
  }
  return true;

} // grid::read

//============================================================================

template<class T>
bool loadpgm(grid<T>& a, const string& pgmname)
{
  // Open the pgm file
  std::ifstream ifs(pgmname.c_str());
  if (!ifs)
    return false;

  // Find the first line that doesn't begin with white space.
  char buf[255];
  char pchar;
  int mode, r, c, maxval, matches;
  do {
    ifs.getline(buf, 255);
  } while (!ifs.eof() && (buf[0] == 0 || buf[0] == ' '));
  if (ifs.eof())
    return false;

  // Make sure the file is a pgm file.  Determine the pgm mode, the
  // dimensions, and the range of pixel values
  matches = sscanf(buf, "%c%1d%d%d%d", &pchar, &mode, &r, &c, &matches);
  if (matches < 2 || mode < 2 || mode > 6)
    return false;

  if (matches < 5) {
    ifs.getline(buf, 255);
    while (!ifs.eof() && (buf[0] == '#' || buf[0] == 0))
    ifs.getline(buf, 255);
    sscanf(buf, "%d %d", &r, &c);
    ifs.getline(buf, 255);
    sscanf(buf, "%d", &maxval);
    if (ifs.eof() || r <= 0 || c <= 0 || maxval <= 0)
      return false;
  }

  a.resize(r, c);
  int i, j;
  switch(mode)
  {
   case 2:
    for(j = 0; j < c; ++j) {
      int tmpi;
      for(i = 0; i < r; ++i) {
        ifs >> tmpi;
        a(i, j) = tmpi;
      }
    };
    break;

   case 3:
    for(j = 0; j < c; j++) {
      int tmpi[3];
      for(i = 0; i < r; ++i) {
        ifs >> tmpi[0] >> tmpi[1] >> tmpi[2];
        a(i, j) = (unsigned char)
          (0.212671 * tmpi[0] + 0.715160 * tmpi[1] + 0.072169 * tmpi[2]);
      }
    };
    break;

   case 5:
    for(j = 0; j < c; ++j) {
      for(i = 0; i < r; ++i)
        ifs.read((char *)&a(i, j), 1);
    };
    break;

   case 6:
    for(j=0; j < c; j++) {
      unsigned char tmpc[3];
      for(i = 0; i < r; ++i) {
        ifs.read((char *)tmpc, 3);
        a(i, j) = (unsigned char)
          (0.212671 * tmpc[0] + 0.715160 * tmpc[1] + 0.072169 * tmpc[2]);
      }
    };
    break;

   default:
    return false;
    break;
  };
  return true;

} // loadpgm

// ==========================================================================
// Set 'g' equal to a subgrid of *this.  ('r', 'c') is the starting point
// in 'g', and the dimension of the subgrid is ('rsize', 'csize').

template<class T>
void grid<T>::subgrid(grid<T>& g, uint r, uint c, uint rsize, uint csize)
{
  if (rsize == 0 || csize == 0 || r >= rows() || c >= cols())
    { g.clear(); return; }

  if (rsize > rows() - r) rsize = rows() - r;
  if (csize > cols() - c) csize = cols() - c;
  if (this == &g)
    { grid tmp; subgrid(tmp, r, c, rsize, csize); swap(tmp); }
  else {
    g.resize(rsize, csize);
    T *gp = &g[0], *thisp = &(operator()(r, c));
    for (uint j = 0; j != csize; ++j) {
      for (uint i = 0; i != rsize; ++i)
        *gp++ = *thisp++;
      thisp += rows() - rsize;
    }
  }

} // grid::subgrid

// ==========================================================================
// Add/Subtract/Multiply/Divide a vector 'v' to/by a vector or grid 'g'.
// Column order, recycle vector elements, so arguments do not need to
// have the same size.

#define VGOP1(VG,OP1,OP2) \
template<class T> VG <T>& \
operator OP1(VG <T>& g, const vector<T>& v) { \
  if (g.size() == 0 || v.size() == 0) return g; \
  T* gp = &g[0]; const T* vp = &v[0]; \
  for (uint i = 0; i != g.size(); ++i) { \
    if (i % v.size() == 0) vp = &v[0]; \
    *gp++ OP1 *vp++; \
  } \
  return g; \
} \
template<class T> VG <T> \
operator OP2(const VG <T>& g, const vector<T>& v) \
  { VG <T> tmp = g; return operator OP1(tmp, v); }

#define VGOP(OP1, OP2) VGOP1(grid, OP1, OP2) VGOP1(vector, OP1, OP2)

VGOP(+=, +) VGOP(-=, -) VGOP(*=, *) VGOP(/=, /)

#undef VGOP1
#undef VGOP

// Add/Subtract/Multiply/Divide a grid 'g2' to/by a grid 'g1'
// For these operations the grids must be conformable.

#define GGOP(OP1,OP2) \
template<class T> grid<T>& \
operator OP1(grid<T>& g1, const grid<T>& g2) { \
  if (g1.rows() != g2.rows() || g1.cols() != g2.cols()) \
    throw typename grid<T>::grid_error("operator" # OP1 ": Grids are not conformable"); \
  T *g1p = &g1[0]; const T* g2p = &g2[0]; \
  for (uint i = 0; i != g1.size(); ++i) *g1p++ OP1 *g2p++; \
  return g1; \
} \
template<class T> grid<T> \
operator OP2(const grid<T>& g1, const grid<T>& g2) \
  { grid<T> tmp = g1; return operator OP1(tmp, g2); }

GGOP(+=, +) GGOP(-=, -) GGOP(*=, *) GGOP(/=, /)

#undef GGOP

// Add/Subtract/Multiply/Divide each element of a vector or grid 'g' to/by
// a fixed scalar value 't'

#define VGSOP1(VG, OP1, OP2) \
template<class T, class S> VG <T>& \
operator OP1(VG <T>& g, const S& t) { \
  T *gp = &g[0]; \
  for (uint i = 0; i != g.size(); ++i) *gp++ OP1 t; \
  return g; \
} \
template<class T, class S> VG <T> \
operator OP2(const VG <T>& g, const S& t) \
  { VG <T> tmp = g; return operator OP1(tmp, t); }

#define VGSOP(OP1, OP2) VGSOP1(vector, OP1, OP2) VGSOP1(grid, OP1, OP2)

VGSOP(+=, +) VGSOP(-=, -) VGSOP(*=, *) VGSOP(/=, /)

#undef VGSOP1
#undef VGSOP

// ==========================================================================
// Mathematical Functions.  Why doesn't 'transform' work for these functions?
// More efficient: function(a, b) --- Less efficient: b = function(a)

// log(grid, grid), log(grid, vector), grid = log(grid), vector = log(vector)
#define VGFN1(FN) \
template <class T> grid <T>& FN(const grid <T>& m1, grid <T>& m2) { \
  m2.resize(m1.rows(), m1.cols()); \
  typename vector<T>::const_iterator m1i = m1.begin(); \
  typename vector<T>::iterator m2i = m2.begin(); \
  while (m1i != m1.end()) *m2i++ = std::FN(*m1i++); \
  return m2; \
} \
template <class T> vector <T>& FN(const grid <T>& m1, vector <T>& m2) { \
  m2.resize(m1.size()); \
  typename vector<T>::const_iterator m1i = m1.begin(); \
  typename vector<T>::iterator m2i = m2.begin(); \
  while (m1i != m1.end()) *m2i++ = std::FN(*m1i++); \
  return m2; \
}

// log(vector, grid), log(vector, vector)
#define VGFN2(VG, FN) \
template <class T> VG <T>& FN(const vector <T>& m1, VG <T>& m2) { \
  m2.resize(m1.size()); \
  typename vector<T>::const_iterator m1i = m1.begin(); \
  typename vector<T>::iterator m2i = m2.begin(); \
  while (m1i != m1.end()) *m2i++ = std::FN(*m1i++); \
  return m2; \
} \
template<class T> VG <T> FN(const VG <T>& m) \
  { VG <T> tmp; return FN(m, tmp); }

#define VGFN(FN) VGFN1(FN) VGFN2(grid, FN) VGFN2(vector, FN)

VGFN(abs) VGFN(exp) VGFN(log) VGFN(log10) VGFN(sin) VGFN(cos) VGFN(tan)
VGFN(asin) VGFN(acos) VGFN(atan) VGFN(sinh) VGFN(cosh) VGFN(tanh)

#undef VGFN1
#undef VGFN2
#undef VGFN

template <class T>
grid<T>& pow(const grid<T>& m1, T p, grid<T>& m2) {
  m2.resize(m1.rows(), m1.cols());
  typename vector<T>::const_iterator m1i = m1.begin();
  typename vector<T>::iterator m2i = m2.begin();
  while (m1i != m1.end()) *m2i++ = std::pow(*m1i++, p);
  return m2;
}

template<class T> grid<T> pow(const grid<T>& m, T p)
  { grid<T> tmp; return pow(m, p, tmp); }

// ==========================================================================
// LU Decomposition with partial pivoting and associated functions.
// For use when LAPACK is not available.  The pivots indexing is
// consistent with LAPACK.

// Can be safely called as LU(a, a), which is efficient but replaces
// 'a' with its LU decomposition.

template<class T>
int LU(const grid<T>& a, grid<T>& g, ivector& pivots,
       int gridtype=gengrid, char ul='U')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("LU: Grid is not square");

  uint n = a.rows();
  pivots.resize(n);
  if (gridtype == posgrid)
    return chol(a, g, ul);

  g = a;
  if (gridtype == trigrid)
    return 0;

  int info = 0;
  for (uint i = 0; i < n - 1; ++i)
  {
    uint maxrow = i;
    double maxval = std::abs(g(i, i));
    for (uint j = i + 1; j < n; ++j)
      if (std::abs(g(j, i)) > maxval) {
        maxval = std::abs(g(j, i));
        maxrow = j;
      }
    if (maxval == 0) {
      info = 1;
      pivots[i] = 1;
    } else {
      pivots[i] = maxrow + 1;
      if (i != maxrow)
        for (uint j = 0; j < n; ++j)
          std::swap(g(i, j), g(maxrow, j));
      for (uint j = i + 1; j < n; ++j) {
        g(j, i) /= g(i, i);
        for (uint k = i + 1; k < n; ++k)
          g(j, k) -= g(j, i) * g(i, k);
      }
    }
  }
  if (g(n - 1, n - 1) == T(0)) {
    info = 1;
    pivots[n - 1] = 1;
  } else
    pivots[n - 1] = n;
  return info;

} // LU

template<class T>
int LUsolve(const grid<T>& a, const grid<T>& b, grid<T>& sol,
    const ivector& pivots, int gridtype=gengrid, char ul='U',
            char tr='N', char dg='N')
{
  if (gridtype == trigrid)
    return trsolve(a, b, sol, ul, tr, dg);

  if (gridtype == posgrid)
    return cholsolve(a, b, sol, ul);

  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("LUsolve: Grid is not square");
  if (a.rows() != pivots.size())
    throw typename grid<T>::grid_error("LUsolve: Pivots vector is the wrong size");
  if (a.rows() != b.rows())
    throw typename grid<T>::grid_error("LUsolve: Grids are not conformable");

  if (&sol == &a || &sol == &b) {
    grid<T> tmp;
    LUsolve(a, b, tmp, pivots, gridtype, ul, tr, dg);
    sol << tmp;
    return 0;
  }

  if (&a == &b) {
    sol = grid<T>("I", a.rows());
    return 0;
  }

  sol = b;
  for (uint i = 0; i < pivots.size(); ++i)
    if ((uint)pivots[i] != i + 1)
     for (uint j = 0; j != sol.cols(); ++j)
        std::swap(sol(i, j), sol(pivots[i] - 1, j));

  grid<T> y(sol.rows(), sol.cols());
  for (uint j = 0; j != sol.cols(); ++j) {
    for (uint i = 0; i != sol.rows(); ++i) {
      y(i, j) = sol(i, j);
      for (uint k = 0; k < i; ++k)
        y(i, j) -= a(i, k) * y(k, j);
    }
  }

  for (uint j = 0; j != sol.cols(); ++j) {
    for (uint i = sol.rows(); i > 0; --i) {
      sol(i - 1, j) = y(i - 1, j) / a(i - 1, i - 1);
      for (uint k = i; k != sol.rows(); ++k)
        sol(i - 1, j) -= a(i - 1, k) * sol(k, j) / a(i - 1, i - 1);
    }
  }
  return 0;
}

template<class T>
int solve(const grid<T>& a, const grid<T>& b, grid<T>& sol,
          int gridtype=gengrid, char ul='U', char tr='N', char dg='N')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("solve: Grid is not square");
  if (a.rows() != b.rows())
    throw typename grid<T>::grid_error("solve: Grids are not conformable");

  ivector pivots; int info = 0;
  if (gridtype == gengrid || gridtype == symgrid || gridtype == posgrid) {
    grid<T> LUthis;
    info = LU(a, LUthis, pivots, gridtype, ul);
    if (info == 0)
      info = LUsolve(LUthis, b, sol, pivots, gridtype, ul, tr, dg);
  } else if (gridtype == trigrid) {
    info = LUsolve(a, b, sol, pivots, gridtype, ul, tr, dg);
  }
  return info;
}

template<class T>
T det(const grid<T>& a, int gridtype=gengrid)
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("det: Grid is not square");
  if (gridtype == trigrid)
    return detLU(a, gridtype);
  if (gridtype == gengrid || gridtype == symgrid || gridtype == posgrid) {
    grid<T> lu; ivector pivots;
    uint info = LU(a, lu, pivots, gridtype);
    if (info != 0) {
      std::cerr << "Warning: LU routine returned info=" << info
           << " with gridtype=" << gridtype << " in det()\n";
      return 0;
    }
    return detLU(lu, gridtype, pivots);
  } else
    throw typename grid<T>::grid_error("det: Unknown grid type");
  return 0;
}

template<class T>
T detLU(const grid<T>& lu, int gridtype=gengrid, const ivector pivots=ivector())
{
  if (lu.rows() != lu.cols())
    throw typename grid<T>::grid_error("detLU: Grid is not square");
  T res = 1;
  if (gridtype == gengrid || gridtype == symgrid) {
#ifdef USE_LAPACK
    if (gridtype == gengrid) {
#else
    if (gridtype == gengrid || gridtype == symgrid) {
#endif
      int mult = 1;
      for (uint i = 0; i < lu.rows(); ++i) {
        res *= lu(i, i);
        if ((uint)pivots[i] != 1+i) mult = -mult;
      }
      res *= mult;
    } else // THIS DOES NOT WORK!!!
      for (uint i = 0; i < lu.rows(); ++i) res *= lu(i, i);
  } else if (gridtype == posgrid) {
    for (uint i = 0; i < lu.rows(); ++i) res *= lu(i, i);
    res *= res;
  } else if (gridtype == trigrid) {
    for (uint i = 0; i < lu.rows(); ++i) res *= lu(i, i);
  } else
    throw typename grid<T>::grid_error("detLU: Unknown grid type");
  return res;
}

template<class T>
int inv(const grid<T>& a, grid<T>& x, int gridtype=gengrid,
        char ul='U', bool fill=true, char dg='N')
{
  grid<T> lu;
  ivector pivots;
  if (LU(a, lu, pivots, gridtype, ul))
    throw typename grid<T>::grid_error("inv: Grid is singular");
  return LUsolve(lu, grid<T>("I", a.rows()), x, pivots, gridtype, ul);
}

template<class T>
grid<T> inv(const grid<T>& a, int gridtype=gengrid, char ul='U',
                bool fill=true, char dg='N')
  { grid<T> tmp; inv(a, tmp, gridtype, ul, dg); return tmp; }

template<class T> grid<T>&
backsolve(const grid<T>& a, const grid<T>& b, grid<T>& x, char ul='U')
{
  x = b;
  if (ul == 'U') {
    for (uint i = 0; i < x.rows(); ++i) {
      const uint curr = x.rows() - 1 - i;
      if (a(curr, curr) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(curr, j) /= a(curr, curr);
        for (uint k = 0; k != curr; ++k)
          x(k, j) -= a(k, curr) * x(curr, j);
      }
    }
  }
  else {
    for (uint i = 0; i < x.rows(); ++i) {
      const uint curr = x.rows() - 1 - i;
      if (a(curr, curr) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(curr, j) /= a(curr, curr);
        for (uint k = 0; k != curr; ++k)
          x(k, j) -= a(curr, k) * x(curr, j);
      }
    }
  }
  return x;
}

template<class T> grid<T>&
forwardsolve(const grid<T>& a, const grid<T>& b, grid<T>& x, char ul='L')
{
  x = b;
  if (ul == 'L') {
    for (uint i = 0; i < x.rows(); ++i) {
      if (a(i, i) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(i, j) /= a(i, i);
        for (uint k = i + 1; k < x.rows(); ++k)
          x(k, j) -= a(k, i) * x(i, j);
      }
    }
  }
  else {
    for (uint i = 0; i < x.rows(); ++i) {
      const uint curr = i;
      if (a(curr, curr) == T(0))
        return x;
      for (uint j = 0; j < x.cols(); ++j) {
        x(curr, j) /= a(curr, curr);
        for (uint k = i + 1; k < x.rows(); ++k)
          x(k, j) -= a(i, k) * x(i, j);
      }
    }
  }
  return x;
}

// ==========================================================================
// Cholesky Decomposition Functions for use when LAPACK is not available.
// Can be called as chol(a, a), which is efficient but replaces 'a' with
// its Cholesky decomposition.  'ul' can be either 'U' or 'L' to
// indicate whether the upper or lower half is computed.

template<class T>
int chol(const grid<T>& a, grid<T>& chol, char ul='U')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("chol: Grid is not square");
  uint n = a.rows();
  chol = a;
  if (ul == 'U') {
    for (uint i = 0; i != n - 1; ++i)
      memset(&chol[0] + i * n + i + 1, 0, sizeof(T) * (n - i - 1));
  } else {
    for (uint i = 1; i != n; ++i)
      memset(&chol[0] + i * n, 0, sizeof(T) * i);
  }
  if (ul == 'U') {
    for (uint i = 0; i != n; ++i)
    {
      for (uint j = 0; j != i; ++j)
        chol(i, i) -= chol(j, i) * chol(j, i);
      if (chol(i, i) == T(0))
        return i;
      chol(i, i) = std::sqrt((T)chol(i, i));
      for (uint j = i + 1; j != n; ++j) {
        for (uint k = 0; k != i; ++k)
          chol(i, j) -= chol(k, i) * chol(k, j);
        chol(i, j) /= chol(i, i);
      }
    }
  } else {
    for (uint i = 0; i != n; ++i)
    {
      for (uint j = 0; j != i; ++j)
        chol(i, i) -= chol(i, j) * chol(i, j);
      if (chol(i, i) == T(0))
        return i;
      chol(i, i) = std::sqrt((T)chol(i, i));
      for (uint j = i + 1; j != n; ++j) {
        for (uint k = 0; k != i; ++k)
          chol(j, i) -= chol(i, k) * chol(j, k);
        chol(j, i) /= chol(i, i);
      }
    }
  }
  return 0;

} // chol

template<class T>
int cholsolve(const grid<T>& a, const grid<T>& b, grid<T>& x, char ul='U')
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("cholsolve: Grid is not square");
  grid<T> y = b;
  forwardsolve(a, b, y, ul);
  backsolve(a, y, x, ul);
  return 0;

} // cholsolve

template<class T>
int cholupdate(grid<T>& a, vector<T>& x)
  { throw typename grid<T>::grid_error("cholupdate: Calling non-Lapack stub"); return 0; }

template<class T>
int trsolve(const grid<T>& a, const grid<T>& b, grid<T>& sol,
            char ul, char tr, char dg)
{
  if (a.rows() != a.cols())
    throw typename grid<T>::grid_error("trsolve: Grid is not square");
  if (ul == 'U')
    forwardsolve(a, b, sol, ul);
  else
    backsolve(a, b, sol, ul);
  return 0;

} // trsolve

//============================================================================
// Various vector and matrix norms

template<class T>
double l1vecnorm(const grid<T>& a) {
  double sum = 0; const T* tp = &a.front();
  for (uint i = 0; i < a.size(); ++i)
    sum += std::abs(*tp++);
  return sum;
}

template<class T>
double l2vecnorm(const grid<T>& a) {
  double sum = 0, absterm; const T* tp = &a.front();
  for (uint i = 0; i < a.size(); ++i)
    { absterm = std::abs(*tp++); sum += absterm * absterm; }
  return std::sqrt(sum);
}

template<class T>
double linfvecnorm(const grid<T>& a) {
  double max = 0, absterm; const T* tp = &a.front();
  for (uint i = 0; i < a.size(); ++i)
    { if ((absterm = std::abs(*tp++)) > max) max = absterm; }
  return max;
}

template<class T>
double l1matnorm(const grid<T>& a) {
  double max = 0;
  for (uint j = 0; j < a.cols(); ++j) {
    double sum = 0; const T* tp = &a(0, j);
    for (uint i = 0; i < a.rows(); ++i)
      sum += std::abs(*tp++);
    if (sum > max) max = sum;
  }
  return max;
}

template<class T>
double l2matnorm(const grid<T>& a, int gridtype) {
  grid<T> evecs; dgrid evals;
  eigenvectors(a, evecs, evals, gridtype);
  return linfvecnorm(evals);
}

template<class T>
double linfmatnorm(const grid<T>& a) {
  double max = 0;
  for (uint i = 0; i < a.rows(); ++i) {
    double sum = 0;
    for (uint j = 0; j < a.cols(); ++j)
      sum += std::abs(a(i, j));
    if (sum > max) max = sum;
  }
  return max;
}

template<class T>
double frobnorm(const grid<T>& a) { return l2vecnorm(a); }

// ==========================================================================
// Matrix Multiplication Functions.
//
// Multiply two grids and/or vectors together using matrix multiplication
// as is done with '%*%' in R.  Ideally, code that uses this class should
// be compiled using -DUSE_LAPACK and linked with BLAS/LAPACK, but these
// implementations will allow the code to work regardless.

template<class T> grid<T>&
matmult(const grid<T>& a, const grid<T>& b, grid<T>& c,
        char ta='N', char tb='N')
{
  if (a.cols() != b.rows())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if (&c == &a || &c == &b)
    { grid<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(a.rows(), b.cols()); c.fill(0);
    for (uint i = 0; i < a.rows(); ++i)
      for (uint j = 0; j < b.cols(); ++j) {
        const T *ap = &a(i, 0), *bp = &b(0, j);
        T *cp = &c(i, j);
        for (uint k = 0; k < a.cols(); ++k)
          { *cp += *ap * *bp++; ap += a.rows(); }
      }
  }
  return c;
}

template<class T> grid<T>&
matmult(const grid<T>& a, const vector<T>& b, grid<T>& c, char ta='N', char tb='N')
{
  if (a.cols() != b.size())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if (&c == &a || (const vector<T>*)&c == &b)
    { grid<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(a.rows(), 1); c.fill(0);
    const T *ap = &a[0], *bp = &b[0];
    for (uint i = 0; i < a.cols(); ++i) {
      T* cp = &c[0];
      for (uint j = 0; j < a.rows(); ++j) *cp++ += *ap++ * *bp;
      ++bp;
    }
  }
  return c;
}

template<class T> grid<T>&
matmult(const vector<T>& a, const grid<T>& b, grid<T>& c, char ta='N', char tb='N')
{
  if (a.size() != b.rows())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if ((const vector<T>*)&c == &a || &c == &b)
    { grid<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(1, b.cols()); c.fill(0);
    const T *bp = &b[0]; T* cp = &c[0];
    for (uint j = 0; j < b.cols(); ++j) {
      const T *ap = &a[0];
      for (uint i = 0; i < b.rows(); ++i) *cp += *ap++ * *bp++;
      ++cp;
    }
  }
  return c;
}

template<class T> vector<T>&
matmult(const grid<T>& a, const vector<T>& b, vector<T>& c, char ta='N', char tb='N')
{
  if (a.cols() != b.size())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (ta != 'N' || tb != 'N')
    throw typename grid<T>::grid_error("matmult: Transposes not supported in non-BLAS matmult");
  if (&c == (const vector<T>*)&a || &c == &b)
    { vector<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(a.rows()); std::fill(c.begin(), c.end(), 0);
    const T *ap = &a[0], *bp = &b[0];
    for (uint i = 0; i < a.cols(); ++i) {
      T* cp = &c[0];
      for (uint j = 0; j < a.rows(); ++j) *cp++ += *ap++ * *bp;
      ++bp;
    }
  }
  return c;
}

template<class T> vector<T>&
matmult(const vector<T>& a, const grid<T>& b, vector<T>& c, char ta='N', char tb='N')
{
  if (a.size() != b.rows())
    throw typename grid<T>::grid_error("matmult: Grids are not conformable");
  if (&c == &a || &c == (const vector<T>*)&b)
    { vector<T> tmp; c.swap(matmult(a, b, tmp, ta, tb)); }
  else {
    c.resize(b.cols()); std::fill(c.begin(), c.end(), 0);
    T* cp = &c[0];
    const T *bp = &b[0];
    for (uint j = 0; j < b.cols(); ++j) {
      const T* ap = &a[0];
      for (uint i = 0; i < b.rows(); ++i)
        *cp += *ap++ * *bp++;
      ++cp;
    }
  }
  return c;
}

template<class T, class S> grid<T>&
diagmult(const grid<T>& a, const vector<S>& b, grid<T>& c)
{
  if (a.cols() != b.size())
    throw typename grid<T>::grid_error("diagmult: Grids are not conformable");
  if (&c == &a || (const vector<S>*)&c == &b)
    { grid<T> tmp; c.swap(diagmult(a, b, tmp)); }
  else {
    c.resize(a.rows(), a.cols());
    const T* ap = &a[0]; const S* bp = &b[0];
    T* cp = &c[0];
    for (uint i = 0; i < a.size(); ) {
      *cp++ = *ap++ * *bp;
      if (++i % a.rows() == 0) ++bp;
    }
  }
  return c;
}

template<class T, class S> grid<T>&
diagmult(const vector<S>& a, const grid<T>& b, grid<T>& c)
{
  if (a.size() != b.rows())
    throw typename grid<T>::grid_error("diagmult: Grids are not conformable");
  if ((const vector<S>*)&c == &a || &c == &b)
    { grid<T> tmp; c.swap(diagmult(a, b, tmp)); }
  else {
    c.resize(b.rows(), b.cols());
    const S* ap = &a[0]; const T* bp = &b[0];
    T* cp = &c[0];
    for (uint i = 0; i < b.size(); ) {
      *cp++ = *ap++ * *bp++;
      if (++i % b.rows() == 0) ap = &a[0];
    }
  }
  return c;
}

template<class T> grid<T>
matmult(const grid<T>& a, const grid<T>& b)
  { grid<T> tmp; return matmult(a, b, tmp); }

template<class T> grid<T>
matmult(const grid<T>& a, const vector<T>& b)
  { grid<T> tmp; return matmult(a, b, tmp); }

template<class T> grid<T>
matmult(const vector<T>& a, const grid<T>& b)
  { grid<T> tmp; return matmult(a, b, tmp); }

template<class T, class S> grid<T>
diagmult(const grid<T>& a, const vector<S>& b)
  { grid<T> tmp; return diagmult(a, b, tmp); }

template<class T, class S> grid<T>
diagmult(const vector<S>& a, const grid<T>& b)
  { grid<T> tmp; return diagmult(a, b, tmp); }

// ==========================================================================
// Statistical Functions.  Most of these functions accept an argument 'rc'
// that is either "row", "col" or "".  If "row", the next argument(s)
// is the row(s) of the grid that the particular statistical function
// should be applied to.  Similarly, for "col".  If "rc" is empty,
// the function will be applied to all the values in the grid.  If
// an additional argument 'l' is supplied, it is a list (subset)
// of the elements within the row or column that are to be used in
// the computation.

template <class T>
T sum(const grid<T>& g, const string& rc="", uint ind=0)
{
  T total = 0.0;
  const T* curr;
  if (rc == "row") {
    curr = &g(ind, 0);
    for (uint i = 0; i != g.cols(); ++i) {
      total += *curr; curr += g.rows();
    }
  } else if (rc == "col") {
    curr = &g(0, ind);
    total = std::accumulate(curr, curr + g.rows(), T(0));
  } else {
    total = std::accumulate(g.begin(), g.end(), T(0));
  }
  return total;
}

template <class T>
T sum(const grid<T>& g, const string& rc, uint ind, const uivector& l) {
  T total = 0.0;
  if (rc == "row")
    for (uint i = 0; i < l.size(); ++i)
      total += g(ind, l[i]);
  else if (rc == "col")
    for (uint i = 0; i < l.size(); ++i)
      total += g(l[i], ind);
  else
    for (uint i = 0; i < l.size(); ++i)
      total += g(l[i]);
  return total;

}

template <class T>
T mean(const grid<T>& g, const string& rc="", uint ind=0) {
  T total = 0.0;
  if (rc == "row" && g.cols() != 0)
    total = sum(g, rc, ind) / T(g.cols());
  else if (rc == "col" && g.rows() != 0)
    total = sum(g, rc, ind) / T(g.rows());
  else if (g.size() != 0)
    total = sum(g) / T(g.size());
  return total;
}

template <class T>
T mean(const grid<T>& g, const string& rc, uint ind, const uivector& l)
{
  T total = 0.0;
  if (l.size() != 0) {
    total = sum(g, rc, ind, l);
    total /= l.size();
  }
  return total;
}

template <class T>
grid<T>& mean(const grid<T>& g, grid<T>& m, const string& rc="col")
{
  if (rc == "row") m.resize(g.rows());
  else if (rc == "col") m.resize(g.cols());
  else throw typename grid<T>::grid_error("mean: Need to specify rc (row or col)");

  for (uint i = 0; i != m.size(); ++i) m[i] = mean(g, rc, i);
  return m;
}

template <class T>
T var(const grid<T>& g, const string& rc="", uint ind=0)
{
  T total = 0, m = mean(g, rc, ind);
  const T* curr;
  if (rc == "row" && g.cols() > 1) {
    curr = &g(ind, 0);
    for (uint i = 0; i != g.cols(); ++i)
      { total += (*curr - m) * (*curr - m); curr += g.rows(); }
    total /= (g.cols() - 1.0);
  } else if (rc == "col" && g.rows() > 1) {
    curr = &g(0, ind);
    for (uint i = 0; i != g.rows(); ++i)
      { total += (*curr - m) * (*curr - m); curr++; }
    total /= (g.rows() - 1.0);
  } else if (g.size() > 1) {
    curr = &g.front();
    for (uint i = 0; i != g.size(); ++i)
      { total += (*curr - m) * (*curr - m); curr++; }
    total /= (g.size() - 1.0);
  }
  return total;
}

template <class T>
T var(const grid<T>& g, const string& rc, uint ind, const uivector& l)
{
  if (l.size() < 2) return 0.0;
  T total = 0, m = mean(g, rc, ind);
  if (rc == "row")
    for (uint i = 0; i < l.size(); ++i)
      total += (g(ind, l[i]) - m) * (g(ind, l[i]) - m);
  else if (rc == "col")
    for (uint i = 0; i < l.size(); ++i)
      total += (g(l[i], ind) - m) * (g(l[i], ind) - m);
  else
    for (uint i = 0; i < l.size(); ++i)
      total += (g[l[i]] - m) * (g[l[i]] - m);
  total /= (l.size() - 1.0);
  return total;
}

template <class T>
grid<T>& var(const grid<T>& g, grid<T>& v, const string& rc="col")
{
  if (rc == "row") v.resize(g.rows());
  else if (rc == "col") v.resize(g.cols());
  else throw typename grid<T>::grid_error("var: Need to specify rc (row or col)");

  for (uint i = 0; i != v.size(); ++i) v[i] = var(g, rc, i);
  return v;
}

template <class T>
T sd(const grid<T>& g, const string& rc="", uint ind=0)
  { return std::sqrt(var(g, rc, ind)); }

template <class T>
T sd(const grid<T>& g, const string& rc, uint ind, const uivector& l)
  { return std::sqrt(var(g, rc, ind, l)); }

template <class T>
grid<T>& sd(const grid<T>& g, grid<T>& s, const string& rc="col")
{
  if (rc == "row") s.resize(g.rows());
  else if (rc == "col") s.resize(g.cols());
  else throw typename grid<T>::grid_error("sd: Need to specify rc (row or col)");

  for (uint i = 0; i != s.size(); ++i) s[i] = sd(g, rc, i);
  return s;
}

template <class T>
T cov(const grid<T>& g, const string& rc, uint i1, uint i2)
{
  T total = 0.0, m1 = mean(g, rc, i1), m2 = mean(g, rc, i2);
  const T *curr1, *curr2;
  if (rc == "row") {
    curr1 = &g(i1, 0); curr2 = &g(i2, 0);
    for (uint i = 0; i < g.cols(); ++i) {
      total += (*curr1 - m1) * (*curr2 - m2);
      curr1 += g.rows(); curr2 += g.rows();
    }
    total /= (g.cols() - 1.0);
  } else if (rc == "col") {
    curr1 = &g(0, i1); curr2 = &g(0, i2);
    for (uint i = 0; i < g.rows(); ++i)
      total += (*curr1++ - m1) * (*curr2++ - m2);
    total /= (g.rows() - 1.0);
  } else
    throw typename grid<T>::grid_error("cov: Need to specify rc (row or col)");
  return total;
}

template <class T>
T cov(const grid<T>& g, const string& rc, uint i1, uint i2, const uivector& l)
{
  T total = 0.0, m1 = mean(g, rc, i1), m2 = mean(g, rc, i2);
  if (rc == "row") {
    for (uint i = 0; i < l.size(); ++i)
      total += (g(i1, l[i]) - m1) * (g(i2, l[i]) - m2);
    total /= (l.size() - 1.0);
  } else if (rc == "col") {
    for (uint i = 0; i < l.size(); ++i)
      total += (g(l[i], i1) - m1) * (g(l[i], i2) - m2);
    total /= (l.size() - 1.0);
  } else
    throw typename grid<T>::grid_error("cov: Need to specify rc (row or col)");
  return total;
}

template <class T>
grid<T>& cov(const grid<T>& g, grid<T>& d, const string& rc="col")
{
  int nvars = 0, nobs = 0;
  if (rc == "row")
    { nvars = g.rows(); nobs = g.cols(); matmult(g, trans(g), d); }
  else if (rc == "col")
    { nvars = g.cols(); nobs = g.rows(); matmult(trans(g), g, d); }
  else
    throw typename grid<T>::grid_error("cov: Need to specify rc (row or col)");

  grid<T> m(nvars, 1); mean(g, m, rc);
  ((d /= nobs) -= matmult(m, trans(m))) *= (nobs / (nobs - 1.0));
  return d;
}

template <class T>
T cor(const grid<T>& g, const string& rc, uint i1, uint i2)
  { return cov(g, rc, i1, i2) / (sd(g, rc, i1) * sd(g, rc, i2)); }

template <class T>
T cor(const grid<T>& g, const string& rc,
      uint i1, uint i2, const uivector& l)
  { return cor(g, rc, i1, i2, l) / (sd(g, rc, i1, l) * sd(g, rc, i2, l)); }

template <class T>
void cor(const grid<T>& g, grid<T>& d, const string& rc="col")
{
  uint sz = 0;
  if (rc == "row") sz = g.rows();
  else if (rc == "col") sz = g.cols();
  else throw typename grid<T>::grid_error("cor: Need to specify rc (row or col)");
  d.resize(sz, sz);
  for (uint i = 0; i != sz; ++i)
    for (uint j = 0; j != sz; ++j)
      d(i, j) = cor(g, rc, i, j);
}

// ==========================================================================
// Stubs for functions that require BLAS/LAPACK for their implementation.

template<class T>
int eigenvalues(const grid<T>& a, grid<double>& m,
                int gridtype=symgrid, char ul='U')
  { grid<T> e; return eigenvectors(a, e, m, gridtype, ul); }

template<class T>
int eigenvalues(const grid<T>& a, grid<cplx>& m, int gridtype=gengrid,
                char ul='U')
  { grid<T> e; dgrid d; return eigenvectors(a, e, d, gridtype, ul); }

template<class T> int eigenvectors(const grid<T>& a, grid<T>& E, grid<double>& D,
                                   int gridtype=symgrid, char ul='U') {
  throw typename grid<T>::grid_error("eigenvectors: Calling non-Lapack stub");
  return 0;
}

template<class T>
int singularvalues(const grid<T>& a, grid<double>& evals)
  { grid<T> evecs; return singularvectors(a, evecs, evals); }

template<class T>
int singularvectors(const grid<T>& a, grid<T>& E, grid<double>& D) {
  throw typename grid<T>::grid_error("singularvectors: Calling non-Lapack stub");
  return 0;
}

// ==========================================================================
// Specializations to dgrid/zgrid and additional functions for use
// when BLAS/LAPACK libraries are available

#ifdef USE_LAPACK

template<> dgrid& matmult(const dgrid& a, const dgrid& b, dgrid& c, char ta, char tb);
template<> dgrid& matmult(const dgrid& a, const dvector& b, dgrid& c, char ta, char tb);
template<> dgrid& matmult(const dvector& a, const dgrid& b, dgrid& c, char ta, char tb);
template<> dvector& matmult(const dgrid&, const dvector&, dvector& c, char ta, char tb);
template<> dvector& matmult(const dvector&, const dgrid&, dvector& c, char ta, char tb);
template<> zgrid& matmult(const zgrid& a, const zgrid& b, zgrid& c, char ta, char tb);
template<> zgrid& matmult(const zgrid& a, const zvector& b, zgrid& c, char ta, char tb);
template<> zgrid& matmult(const zvector& a, const zgrid& b, zgrid& c, char ta, char tb);
template<> zvector& matmult(const zgrid&, const zvector&, zvector& c, char ta, char tb);
template<> zvector& matmult(const zvector&, const zgrid&, zvector& c, char ta, char tb);

template<> int LU(const dgrid&, dgrid&, ivector&, int, char);
template<> int LUsolve(const dgrid&, const dgrid&, dgrid&, const ivector&,
                       int, char, char, char);
template<> int inv(const dgrid&, dgrid&, int, char, bool, char);
template<> dgrid inv(const dgrid&, int, char, bool, char);
template<> dgrid& backsolve(const dgrid&, const dgrid&, dgrid&, char);
template<> dgrid& forwardsolve(const dgrid&, const dgrid&, dgrid&, char);
template<> int chol(const dgrid&, dgrid&, char);
template<> int cholsolve(const dgrid&, const dgrid&, dgrid&, char);
template<> int cholupdate(dgrid& a, dvector& x);
template<> int cholupdate(zgrid& a, zvector& x);
template<> int trsolve(const dgrid&, const dgrid&, dgrid&, char, char, char);
template<> int eigenvalues(const dgrid&, dgrid&, int, char ul);
template<> int eigenvalues(const dgrid&, zgrid&, int, char ul);
template<> int eigenvectors(const dgrid&, dgrid&, dgrid&,
                            int gridtype, char ul);
template<> int singularvalues(const dgrid&, dgrid&);
template<> int singularvectors(const dgrid&, dgrid&, dgrid&);

template<> int LU(const zgrid&, zgrid&, ivector&, int, char);
template<> int LUsolve(const zgrid&, const zgrid&, zgrid&, const ivector&, int,
  char, char, char);
template<> int inv(const zgrid&, zgrid&, int, char, bool, char);
template<> zgrid inv(const zgrid&, int, char, bool, char);
template<> int chol(const zgrid&, zgrid&, char);
template<> int cholsolve(const zgrid&, const zgrid&, zgrid&, char);
template<> int trsolve(const zgrid&, const zgrid&, zgrid&, char, char, char);

template<> zgrid& backsolve(const zgrid&, const zgrid&, zgrid&, char);
template<> zgrid& forwardsolve(const zgrid&, const zgrid&, zgrid&, char);
template<> int eigenvalues(const zgrid&, dgrid&, int, char ul);
template<> int eigenvalues(const zgrid&, zgrid&, int, char ul);
template<> int eigenvectors(const zgrid&, zgrid&, dgrid&, int, char ul);
template<> int singularvalues(const zgrid&, dgrid&);
template<> int singularvectors(const zgrid&, zgrid&, dgrid&);

// ==========================================================================
// BLAS routines

dgrid& axpy(const dgrid& x, double alpha, dgrid& y);
dgrid& gemm(const dgrid& a, const dgrid& b, dgrid& c, double alpha=1,
            double beta=0, char tra='N', char trb='N');
dgrid& trmm(const dgrid& a, dgrid& b, double alpha=1, const char side='L',
            const char uplo='U', const char tra='N', const char diag='N');
dgrid& trsm(const dgrid& a, const dgrid& b, dgrid& x, double alpha,
            char sd='L', char ul='U', char tr='N', char dg='N');

zgrid& axpy(const zgrid& x, const cplx alpha, zgrid& y);
zgrid& gemm(const zgrid& a, const zgrid& b, zgrid& c, const cplx& alpha=1,
            const cplx& beta=0, char tra='N', char trb='N');
zgrid& trmm(const zgrid& a, zgrid& b, const cplx& alpha=1, const char side='L',
            const char uplo='U', const char tra='N', const char diag='N');
zgrid& trsm(const zgrid& a, const zgrid& b, zgrid& x, const cplx& alpha,
            char sd='L', char ul='U', char tr='N', char dg='N');

// ==========================================================================
// LAPACK routines

// Routines to solve AX=B.  A is replaced by some LU type of
// decomposition, while B is replaced by the solution X.

int gels(dgrid& a, dgrid& b, char trans='N');
int gesv(dgrid& a, dgrid& b, ivector& pivots);
int sysv(dgrid& a, dgrid& b, ivector& pivots, char ul);
int posv(dgrid& a, dgrid& b, char ul);

int gels(zgrid& a, zgrid& b, char trans='N');
int gesv(zgrid& a, zgrid& b, ivector& pivots);
int sysv(zgrid& a, zgrid& b, ivector& pivots, char ul);
int hesv(zgrid& a, zgrid& b, ivector& pivots, char ul);
int posv(zgrid& a, zgrid& b, char ul);

#endif

template class grid<unsigned char>;
template class grid<int>;
template class grid<double>;
template class grid<std::complex>;

*/
