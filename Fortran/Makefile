.SUFFIXES:
.SUFFIXES: .o .f90

FC = gfortran
LDLIBS =
LDFLAGS =
FFLAGS = -O2

WINDSRCS = run_wind_enkf.f90 io_module.f90 random_module.f90 enkf_module.f90 eakf_module.f90
WINDOBJS = $(WINDSRCS:.f90=.o)
L96VARSRCS = run_l96_var.f90 l96_module.f90 rk4_module.f90 io_module.f90 random_module.f90
L96VAROBJS = $(L96VARSRCS:.f90=.o)
TESTSRCS = test.f90 io_module.f90 random_module.f90 rk4_module.f90 l96_module.f90
TESTOBJS = $(TESTSRCS:.f90=.o)
EXECS = run_wind_enkf run_l96_var test

all : $(EXECS)

run_wind_enkf: $(WINDOBJS)
	$(FC) $(LDFLAGS) -o $@ $(WINDOBJS) $(LDLIBS)

run_l96_var: $(L96VAROBJS)
	$(FC) $(LDFLAGS) -o $@ $(L96VAROBJS) $(LDLIBS)
	
test : $(TESTOBJS)
	$(FC) $(LDFLAGS) -o $@ $(TESTOBJS) $(LDLIBS)

clean:
	rm -f $(WINDOBJS) $(L96VAROBJS) $(EXECS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

l96_module.o : rk4_module.o random_module.o
run_l96_var.o : io_module.o random_module.o l96_module.o
run_wind_enkf.o : io_module.o random_module.o enkf_module.o eakf_module.o
test.o : io_module.o random_module.o l96_module.o
