LIBFITSRC = ../lib/dierckx

#------------------------------------
# default Unix fortran
FC = f77
FTNLIB =
FFLAGS = -O

#------------------------------------
# Intel Fortran Compiler
#FC = ifort
#FTNLIB = -Vaxlib
#FTNLIB = -Vaxlib /usr/lib/C-ctype.o /usr/lib/C_name.o /usr/lib/ctype-info.o

#------------------------------------
#FC = g77
#FFLAGS = -O -dbl

#------------------------------------
# Dependencies
FITOBJECTS = curfit.o  fpbspl.o  fpcurf.o  fpgivs.o \
	fprati.o  splev.o splder.o \
	fpback.o  fpchec.o  fpdisc.o  \
	fpknot.o  fprota.o  splder.o

#------------------------------------

fitlib: $(FITOBJECTS)
	ar rvs -o fitlib.a $(FITOBJECTS)

mncurf: mncurf.o $(FITOBJECTS)
	$(FC) -o $@ mncurf.o $(FITOBJECTS)

#-----------
%.o: $(LIBFITSRC)/%.f
	$(FC) -c $(FFLAGS) $<

.PHONY: clean

clean:
	rm -fr *.o mncurf
