# CPMOD makefile

# Path to object files, source files and executable
OBJ=.
SRC=../src
DIR=..

# Name of the executable
EXE=cpmod

# Compilation with GNU Fortran
# Settings for release (optimised):
FC = gfortran
FLAGS = -O3 -ffree-line-length-none -cpp
# Settings for debug:
# FLAGS = -O0 -fcheck=all -g -fbacktrace -ffree-line-length-none -cpp

# Code
ICE = module_thermo.o module_ice_microphysics.o
MODPBE = PBE_general.o PBE_discrete.o
MAIN = CPMOD_pbe.o
PBE = PBE_solver.o PBE_growth.o PBE_agg.o PBE_agg_kernels.o \
PBE_test.o PBE_frag.o
PSR = PSR_pbe.o

# Compilation order
MODULES= $(MODPBE) $(ICE)  
OBJECTS= $(MODULES) $(MAIN) $(PBE) $(PSR)

# Make options
all: $(EXE)

# Build executable & place it in directory
$(EXE): $(OBJECTS)
	$(FC) $(FLAGS) $(addprefix $(OBJ)/,$(OBJECTS)) -o $@
	mv $(EXE) $(DIR)
	@echo Creating $(EXE) ..

$(OBJ)/%.o: $(SRC)/%.f
	$(FC) $(FLAGS) -c $<;

$(OBJ)/%.o: $(SRC)/%.f90
	$(FC) $(FLAGS) -c $<;

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod ../cpmod
