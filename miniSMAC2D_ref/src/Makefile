EXE = smac2d

FC=mpif90
FFLAGS=-O3 -c -Dmpimode -DD_PRECISION -g
CPP=/lib/cpp
CPPFLAGS=-P -traditional -Dmpimode -DD_PRECISION
LD=mpif90
LDFLAGS=-O3
INCDIR=../include
MPIDIR=-I/home/projects/openmpi/1.7.3/gnu/4.7.2/include

# delete all previous suffixes
.SUFFIXES:
# assign new suffixes for defining transformation rule '.f.o'
.SUFFIXES: .f .o

OBJ =   bc.o \
	btlu3.o \
        flux.o \
        geom.o \
        initia.o \
        lreltur.o \
        metric.o \
        rbc.o \
        smac2d.o \
        tursa.o \
        blas.o \
        conup.o \
        fm.o \
	ilusol.o \
        lrelax.o \
        lusgs.o \
        preproc.o \
        rgrid.o \
        step.o \
        visc.o

#---------------- compilation rules ----------------------
.f.o:
	$(CPP) $(CPPFLAGS) $*.f -I$(INCDIR) $(MPIDIR) > $*F.f
	$(FC) $(FFLAGS) $(*F)F.f -o $*.o
	rm $*F.f

$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ)

# clean up
clean:
	rm -f *F.f *.o *.f~ *F.f.*
