#USE_OpenMP = TRUE
USE_OpenMP = FALSE

FC = gfortran

ifeq ($(USE_OpenMP),TRUE)
	FLAGS = -Wall -fopenmp -cpp
	SUFFIX = _OpenMP
else
	FLAGS = -Wall -cpp
	SUFFIX = 
endif

MAIN = ComputeCharacteristicStrainForBHBmergers

SRC = UtilitiesModule.f90 $(MAIN).f90
OBJ = UtilitiesModule.o $(MAIN).o

all: $(MAIN)

$(MAIN):
	$(FC) $(FLAGS) -c $(SRC)
	$(FC) $(FLAGS) -o $(MAIN)$(SUFFIX) $(OBJ)

clean:
	rm -f $(MAIN)$(SUFFIX) *.o *.mod
