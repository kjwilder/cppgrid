CC=g++
#CC=g++-3.4
#CC=g++-4.0 -pedantic -ansi
#CC=g++-2.95 -pedantic -ansi

#CCOPTS = -Wall -g -DUSE_ASSERTS
CCOPTS = -Wall -O
#CCOPTS = -Wall -O3 -march=pentium4
#CCOPTS += -DGRIDTEST

# Blas/Lapack options
CCOPTS += -DUSE_LAPACK
LDLIBS = -L/usr/lib/ATLAS -llapack -lf77blas -lcblas -latlas -lg2c
LDLIBS = -L/usr/lib/atlas/sse2 -llapack -lblas

HEADERS = grid.h 
SOURCES = grid.C
OBJECTS = $(SOURCES:.C=.o)

# Uncomment the next line to compile the examples using real grids
# instead of complex grids.
CCOPTS += -DDGRID
CCSWIG = -I/usr/include/python2.3 -DUSE_ASSERTS -DUSE_LAPACK

TESTSRCS = gt_mult.C gt_norm.C gt_eigen.C gt_rosser.C gt_arith.C gt_io.C \
           gt_lu.C gt_chol.C gt_math.C gt_inv.C gt_det.C gt_stat.C
TESTOBJS = $(TESTSRCS:.C=.o)
TESTEXECS = $(TESTSRCS:.C=)

TARLIST = $(TESTSRCS) $(SOURCES) $(HEADERS) grid.i Makefile
TARFILES = $(TARLIST:%=grid/%)

#all : grid.o _grid.so
all : $(TESTEXECS)

%.o : %.C $(HEADERS)
	$(CC) $(CCOPTS) -c $<

gt_mult : gt_mult.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_norm : gt_norm.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_rosser : gt_rosser.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_eigen : gt_eigen.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_arith : gt_arith.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_io : gt_io.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_stat : gt_stat.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_covtest : gt_covtest.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_chol : gt_chol.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_lu : gt_lu.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_math : gt_math.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_inv : gt_inv.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

gt_det : gt_det.o grid.o $(HEADERS)
	$(CC) $(CCOPTS) $< $(OBJECTS) -o $@ $(LDLIBS)

#$(TESTSRCS) : $<.o

clean ::
	/bin/rm -f $(OBJECTS) $(TESTEXECS) $(TESTOBJS) joe joegrids

checkin ::
	ci -l grid.i $(TESTSRCS) $(SOURCES) $(HEADERS) Makefile

archive : checkin
	ln -s . grid
	tar chf grid.tar $(TARFILES)
	rm grid

grid : $(OBJECTS)
	$(CC) $(CCOPTS) $(OBJECTS) -o grid $(LDLIBS)

swig : _grid.so

_grid.so : grid_wrap.o $(OBJECTS)
	$(CC) -Wall -shared grid_wrap.o $(OBJECTS) -o _grid.so $(LDLIBS)

grid_wrap.o : grid_wrap.cxx
	$(CC) $(CCSWIG) -c grid_wrap.cxx

grid_wrap.cxx : grid.i $(HEADERS)
	swig $(CCSWIG) -c++ -python grid.i

swigclean ::
	/bin/rm -f grid_wrap* grid.py* _grid.so

