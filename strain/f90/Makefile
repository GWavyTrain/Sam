USE_OpenMP = TRUE
#USE_OpenMP = FALSE

FC = gfortran

ifeq ($(USE_OpenMP),TRUE)
	FLAGS = -Wall -fopenmp -cpp
else
	FLAGS = -Wall -cpp
endif

MAIN = ComputeCharacteristicStrainForBHBmergers

all: $(MAIN)

ifeq ($(USE_OpenMP),TRUE)
$(MAIN):
	$(FC) $(FLAGS) -o $(MAIN)_OpenMP $(MAIN).f90
else
$(MAIN):
	$(FC) $(FLAGS) -o $(MAIN) $(MAIN).f90
endif

clean:
	rm -f $(MAIN)
