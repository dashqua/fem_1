FC=gfortran  # Compilation
LD=gfortran  # Linker
LIBS=  -L./libraries -lsparsekit
INCLUDE= -I./libraries 
FFLAGS = -fcheck=all -O3 -unroll=4
EXE=pw2mts.exe

OBJS = EFLib.o mod_lec_fic.o pw2mts.o

all:   $(OBJS) 
	@echo "Links edition"
	$(LD) $(INCLUDE) $(LDFLAGS) $(OBJS) $(LIBS)  -o $(EXE) 

pw2mts.o : pw2mts.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c pw2mts.f90

EFLib.o : EFLib.f90 mod_lec_fic.o
	$(FC) $(INCLUDE) $(FFLAGS) -c EFLib.f90

mod_lec_fic.o : mod_lec_fic.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c mod_lec_fic.f90

clean: 
	rm *~ -f *.o *.mod *.dat $(EXE) 
 
%.o : %.mod 

