%module grid
%{
#include "grid.h"
%}

%include "std_string.i"
%include "std_vector.i"

%rename(opequal) operator=(const grid<T>& g);
%ignore grid::operator=;
%ignore operator<<;
%ignore operator[](uint r);

%extend grid { T __getitem__(int x) { return self->operator[](x); } };
%extend grid { void __setitem__(int x, int z=1, T y=0) { self->operator()(x,z) = y; } };
%extend grid { string __repr__()
{ std::stringstream s;
  for (uint i = 0; i < self->rows(); ++i) {
    for (uint j = 0; j < self->cols(); ++j)
      s << self->operator()(i, j) << " ";
    s << "\n";
  }
  return s.str();}
};

%include "grid.h"

%template(ucvector) std::vector<unsigned char>;
%template(ucgrid) grid<unsigned char>;

%template(dvector) std::vector<double>;
%template(dgrid) grid<double>;

