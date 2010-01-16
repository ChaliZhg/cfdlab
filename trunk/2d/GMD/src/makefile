FC = ifort
FFLAGS = -diag-disable vec -r8 -free -Tf

TARGETS = vort

ALL: $(TARGETS)

SRC = $(wildcard *.f95)
OBJ = $(patsubst %.f95,%.o,$(SRC))

vort: $(OBJ)
	$(FC) -o vort $(OBJ)

%.o: %.f95
	$(FC) -c $(FFLAGS) $*.f95

clean:
	rm -f $(OBJ) $(TARGETS) *.mod
