#-----------------------------------------------------------
#  Makes test* linearization-checking programs
#-----------------------------------------------------------
#
# default Unix fortran
FC = f77
FTNLIB =

#------------------------------------------
# Intel fortran compiler
FC = ifort
FTNLIB = -Vaxlib 

#------------------------------------------

test1: test1.o gvcalc.o
	$(FC) -o test1 test1.o gvcalc.o $(FTNLIB)

test2: test2.o gvcalc.o cdfun.o tqcalc.o
	$(FC) -o test2 test2.o gvcalc.o cdfun.o tqcalc.o $(FTNLIB)

testm: testm.o motor.o
	$(FC) -o testm testm.o motor.o $(FTNLIB)

testc: testc.o cdfun.o
	$(FC) -o testc testc.o cdfun.o $(FTNLIB)

test1.o: test1.f
	$(FC) -c -O1 -r8 test1.f

test2.o: test2.f
	$(FC) -c -O1 -r8 test2.f

testc.o: testc.f
	$(FC) -c -O1 -r8 testc.f

testm.o: testm.f
	$(FC) -c -O1 -r8 testm.f

gvcalc.o: gvcalc.f
	$(FC) -c -O1 -r8 gvcalc.f

cdfun.o: cdfun.f
	$(FC) -c -O1 -r8 cdfun.f

tqcalc.o: tqcalc.f
	$(FC) -c -O1 -r8 tqcalc.f

motor.o: motor.f
	$(FC) -c -O1 -r8 motor.f
