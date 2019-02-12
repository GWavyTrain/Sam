USE_OpenMP = TRUE
#USE_OpenMP = FALSE

FC = gfortran

ifeq ($(USE_OpenMP),TRUE)
	FLAGS = -Wall -fopenmp -cpp
else
	FLAGS = -Wall -cpp
endif

MAIN = ComputeCharacteristicStrainForBHBmergers

SRC = UtilitiesModule.f90 $(MAIN).f90
OBJ = UtilitiesModule.o $(MAIN).o

all: $(MAIN)

$(MAIN):
	$(FC) $(FLAGS) -c $(SRC)
	$(FC) $(FLAGS) -o $(MAIN) $(OBJ)

clean:
	rm -f $(MAIN) *.o *.mod
