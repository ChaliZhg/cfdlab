include ../makefile.in

PROGRAM=fvm

SRC=$(wildcard *.f)
OBJ=$(patsubst %.f,%.o,$(SRC))

OBJS=areas.o \
	  common.o \
	  conres.o \
     convert.o \
     check_positivity.o \
	  datastruct.o\
	  farfield_flux.o \
     finalize.o \
     flength.o \
     hllc_flux.o \
     kfvs_flux.o \
     roe_flux.o \
     geometric.o \
	  gradient.o \
     green_gauss.o \
     isocont.o \
     main.o \
     mayavi.o \
	  misc.o \
     prep_gnuplot.o \
     recon.o \
     residu.o \
     result.o \
     sa_model.o \
	  solid_flux.o \
     tri_coord.o \
     update.o \
     vigie.o \
     viscosity.o \
     wall_dist.o \
     check_string.o \
     read_input.o \
     read_grid.o \
     initialize.o \
     clcd.o \
     time_step.o \

.f.o:
	$(F77) $(CFLAGS) -c $*.f

$(PROGRAM) :$(OBJ)
	              $(F77) -o $(PROGRAM) $(OBJ)

clean:;    rm -f $(OBJ) core
