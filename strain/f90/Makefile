FLAGS  = -cpp

ifeq ($(DEBUG),TRUE)
    FLAGS += -Wall
endif

ifeq ($(USE_OMP),TRUE)
    FLAGS += -fopenmp
endif

MAIN = ComputeCharacteristicStrainForBHBmergers
TEST = ModuleTester

SRC  = $(DIR)/UtilitiesModule.f90 $(DIR)/$(MAIN).f90
OBJ  = UtilitiesModule.o $(MAIN).o
SRC2 = $(DIR)/UtilitiesModule.f90 $(DIR)/$(TEST).f90
OBJ2 = UtilitiesModule.o $(TEST).o

all: $(MAIN) $(TEST)

$(MAIN):
	$(FC) $(FLAGS) -c $(SRC)
	$(FC) $(FLAGS) -o $(MAIN).exe $(OBJ)

$(TEST):
	$(FC) $(FLAGS) -c $(SRC2)
	$(FC) $(FLAGS) -o $(TEST).exe $(OBJ2)

clean:
	rm -rf *.exe *.dSYM *.o *.mod
