all : libsparsekit.a

.SUFFIXES: .f90 

# COMPILATEUR
FC = gfortran

# Options d'optimisation
OPTFLAGS = 

# Option de deboggage
DEBUGFLAGS =
#

FLAGS = $(DEBUGFLAGS)
FLAGS = $(OPTFLAGS)

# Repertoire contenant les librairies LAPACK  (64 bits)
#LPATH = -L/Developer/opt/intel/composerxe/mkl/lib

# Liste des librairies utilisees (64 bits)
#LIB = -lmkl_lapack -lmkl_em64t -lmkl_core -lpthread #-liomp5  
LIB =


# Repertoire contenant les libariries LAPACK  (32 bits)
#BLAS_PATH = /opt/intel/mkl/10.0.011/lib/32
#BLAS_LIB = -L$(BLAS_PATH) -lmkl_intel -lmkl_core -lguide -lpthread -Wl,-rpath=$(BLAS_PATH)
#BLAS_SEQ_LIB = -L$(BLAS_PATH) -lmkl_intel -lmkl_sequential -lmkl_core -Wl,-rpath=$(BLAS_PATH)

#BLAS_PATH = -L/usr/lib                         # ifc7
#BLAS_LIB = -latlas -lblas                      # ifc7


# MKL PATH and INCLUDE
#MKL_PATH = /Developer/opt/intel/composerxe/mkl
#MKL_INC = $(MKL_PATH)/include

# MKL SOLVERS (PARDISO, DSS, ISS)
MKL_SOLVER_LIB = 

# LAPACK 
#LAPACK_PATH = $(MKL_PATH)
#LAPACK_LIB = -L$(LAPACK_PATH) -lmkl_lapack -Wl,-rpath=$(LAPACK_PATH)
LAPACK_PATH = 
LAPACK_LIB = 

#LAPACK_PATH = $(BLAS_PATH)                             # ifc7
#LAPACK_LIB = -L$(LAPACK_PATH) -llapack                 # ifc7




########################################################################

LIBS = 

libobj = f90_kind.o tri.o unary_mod.o formats.o blas1.o blassm.o ilut.o \
	 iters.o gestionerr.o initsolver.o resol.o matcsr2.o genmat.o functns.o

libsparsekit.a : f90_kind.o $(libobj)
#	libtool -static  -o libsparsekit.a $(libobj)
	ar -r  libsparsekit.a $(libobj)
#	ranlib libsparsekit.a

tri.o : tri.f90
	$(FC) $(FLAGS)  -c tri.f90

unary_mod.o : unary_mod.f
	$(FC) $(FLAGS) -FI -c ${@:.o=.f}

genmat.o : genmat.f
	$(FC) $(FLAGS) -FI -c ${@:.o=.f}

iters.o: iters.f
	$(FC) $(FLAGS) -c ${@:.o=.f}

blas1.o : blas1.f
	$(FC) $(FLAGS) -c ${@:.o=.f}

blassm.o : blassm.f
	$(FC) $(FLAGS) -c ${@:.o=.f}

ilut.o: ilut.f
	$(FC) $(FLAGS) -c ${@:.o=.f}

formats.o : formats.f
	$(FC) $(FLAGS) -FI -c ${@:.o=.f}

functns.o : functns.f
	$(FC) $(FLAGS) -FI -c ${@:.o=.f}

resol.o : resol.f90
	$(FC) $(FLAGS) -c ${@:.o=.f90}

f90_kind.o : f90_kind.f90
	$(FC) $(FLAGS) -c ${@:.o=.f90}

matcsr2.o : matcsr2.f90
	$(FC) $(FLAGS) -c ${@:.o=.f90}

initsolver.o : initsolver.f90
	$(FC) $(FLAGS) -c ${@:.o=.f90}

gestionerr.o : gestionerr.f90
	$(FC) $(FLAGS) -c ${@:.o=.f90}


clean :
	rm -f *.o *.mod *.a
cleanall :
	make clean; rm -f *~ 
echo :
	@echo "MKL_SOLVER_LIB :" $(MKL_SOLVER_LIB)
	@echo "LAPACK_LIB     :" $(LAPACK_LIB)
