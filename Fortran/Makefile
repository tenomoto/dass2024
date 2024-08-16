.SUFFIXES:
.SUFFIXES: .o .f90

FC = gfortran
LDLIBS =
LDFLAGS =
FFLAGS = -O2

WINDSRCS = run_wind_ens.f90 random_module.f90 enkf_module.f90 eakf_module.f90
WINDOBJS = $(WINDSRCS:.f90=.o)
L63ENSSRCS = run_l63_ens.f90 l63_module.f90 step_module.o rk4_module.f90 random_module.f90 enkf_module.f90 eakf_module.f90
L63ENSOBJS = $(L63ENSSRCS:.f90=.o)
L63VARSRCS = run_l63_var.f90 l63_module.f90 step_module.o rk4_module.f90 random_module.f90
L63VAROBJS = $(L63VARSRCS:.f90=.o)
L96PREPSRCS = prep_l96.f90 l96_module.f90 step_module.o rk4_module.f90 random_module.f90
L96PREPOBJS = $(L96PREPSRCS:.f90=.o)
L96VARSRCS = run_l96_var.f90 l96_module.f90 step_module.o rk4_module.f90 random_module.f90
L96VAROBJS = $(L96VARSRCS:.f90=.o)
L96ENSSRCS = run_l96_ens.f90 l96_module.f90 step_module.o rk4_module.f90 random_module.f90 enkf_module.f90 eakf_module.f90
L96ENSOBJS = $(L96ENSSRCS:.f90=.o)
TESTSRCS = test.f90  random_module.f90 l63_module.f90 l96_module.f90 step_module.o rk4_module.o io_module.f90
TESTOBJS = $(TESTSRCS:.f90=.o)
EXECS = run_wind_ens run_l63_ens run_l63_var prep_l96 run_l96_var run_l96_ens test

all : $(EXECS)

run_wind_ens: $(WINDOBJS)
	$(FC) $(LDFLAGS) -o $@ $(WINDOBJS) $(LDLIBS)

run_l63_ens: $(L63ENSOBJS)
	$(FC) $(LDFLAGS) -o $@ $(L63ENSOBJS) $(LDLIBS)

run_l63_var: $(L63VAROBJS)
	$(FC) $(LDFLAGS) -o $@ $(L63VAROBJS) $(LDLIBS)

prep_l96: $(L96PREPOBJS)
	$(FC) $(LDFLAGS) -o $@ $(L96PREPOBJS) $(LDLIBS)

run_l96_var: $(L96VAROBJS)
	$(FC) $(LDFLAGS) -o $@ $(L96VAROBJS) $(LDLIBS)
	
run_l96_ens: $(L96ENSOBJS)
	$(FC) $(LDFLAGS) -o $@ $(L96ENSOBJS) $(LDLIBS)
	
test : $(TESTOBJS)
	$(FC) $(LDFLAGS) -o $@ $(TESTOBJS) $(LDLIBS)

clean:
	rm -f *.o *.mod $(EXECS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

l63_module.o : step_module.o random_module.o
l96_module.o : step_module.o random_module.o
step_module.o : rk4_module.o
run_l63_ens.o : l63_module.o enkf_module.o eakf_module.o random_module.o step_module.o
run_l63_var.o : l63_module.o random_module.o step_module.o
prep_l96.o : l96_module.o rk4_module.o  random_module.o
run_l96_var.o : l96_module.o step_module.o random_module.o step_module.o
run_l96_ens.o : l96_module.o enkf_module.o eakf_module.o
run_wind_ens.o :  random_module.o enkf_module.o eakf_module.o
test.o :  random_module.o l63_module.o l96_module.o io_module.o
