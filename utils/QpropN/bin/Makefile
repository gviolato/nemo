SRC = ../src
VLSRC = ../vlsrc
FFLAGS = -O 

#------------------------------------
# default Unix fortran
FC = f77
FTNLIB =

#------------------------------------
# Intel Fortran Compiler
#FC = ifort
#FTNLIB = -Vaxlib
#FTNLIB = -Vaxlib /usr/lib/C-ctype.o /usr/lib/C_name.o /usr/lib/ctype-info.o

#------------------------------------
#FC = g77
#FFLAGS = -O -dbl

#------------------------------------

qprop: polardatamodule.mod qprop.o polar.o gvcalc.o cdfun.o tqcalc.o motor.o spline.o io.o qcget.o
	$(FC) -o qprop qprop.o gvcalc.o cdfun.o tqcalc.o motor.o \
spline.o io.o qcget.o polar.o polardata.o $(FTNLIB)

qmil: qmil.o tpdes.o gvcalc.o cdfun.o tqcalc.o bnsolv.o \
spline.o io.o qcget.o
	$(FC) -o qmil qmil.o tpdes.o gvcalc.o cdfun.o tqcalc.o bnsolv.o \
spline.o io.o qcget.o $(FTNLIB)

vlprop: vlprop.o cdfun.o vlset.o motor.o aic.o \
spline.o io.o qcget.o aread.o lefinda.o
	$(FC) -o vlprop vlprop.o cdfun.o vlset.o motor.o aic.o \
spline.o io.o qcget.o aread.o lefinda.o $(FTNLIB)



qprop.o: $(SRC)/qprop.f $(SRC)/QDEF.INC
	$(FC) -c $(FFLAGS) $(SRC)/qprop.f

qmil.o: $(SRC)/qmil.f $(SRC)/QDEF.INC
	$(FC) -c $(FFLAGS) $(SRC)/qmil.f

tpdes.o: $(SRC)/tpdes.f
	$(FC) -c $(FFLAGS) $(SRC)/tpdes.f

gvcalc.o: $(SRC)/gvcalc.f
	$(FC) -c $(FFLAGS) $(SRC)/gvcalc.f

cdfun.o: $(SRC)/cdfun.f
	$(FC) -c $(FFLAGS) $(SRC)/cdfun.f

tqcalc.o: $(SRC)/tqcalc.f
	$(FC) -c $(FFLAGS) $(SRC)/tqcalc.f

polardatamodule.mod: $(SRC)/polardata.f
	$(FC) -c $(SRC)/polardata.f

polar.o: $(SRC)/polar.f
	$(FC) -c $(FFLAGS) $(SRC)/polar.f

motor.o: $(SRC)/motor.f
	$(FC) -c $(FFLAGS) $(SRC)/motor.f

bnsolv.o: $(SRC)/bnsolv.f
	$(FC) -c $(FFLAGS) $(SRC)/bnsolv.f

spline.o: $(SRC)/spline.f
	$(FC) -c $(FFLAGS) $(SRC)/spline.f

io.o: $(SRC)/io.f
	$(FC) -c $(FFLAGS) $(SRC)/io.f

qcget.o: $(SRC)/qcget.f
	$(FC) -c $(FFLAGS) $(SRC)/qcget.f


vlprop.o: $(VLSRC)/vlprop.f
	$(FC) -c $(FFLAGS) $(VLSRC)/vlprop.f

vlset.o: $(VLSRC)/vlset.f
	$(FC) -c $(FFLAGS) $(VLSRC)/vlset.f

helxyz.o: $(VLSRC)/helxyz.f
	$(FC) -c $(FFLAGS) $(VLSRC)/helxyz.f

helcam.o: $(VLSRC)/helcam.f
	$(FC) -c $(FFLAGS) $(VLSRC)/helcam.f

aread.o: $(VLSRC)/aread.f
	$(FC) -c $(FFLAGS) $(VLSRC)/aread.f

lefinda.o: $(VLSRC)/lefinda.f
	$(FC) -c $(FFLAGS) $(VLSRC)/lefinda.f

aic.o: $(VLSRC)/aic.f
	$(FC) -c $(FFLAGS) $(VLSRC)/aic.f

.PHONY: clean

clean:
	rm -fr *.o *.mod qprop
