FC = gfortran

all: grid

grid: naca.pts naca.ctr
	$(TAXIS_HOME)/extern/delaundo/std/delaundo

naca.pts: naca.f
	$(FC) -o naca naca.f
	./naca

clean:
	rm -f naca naca.pts naca.dpl naca.log
