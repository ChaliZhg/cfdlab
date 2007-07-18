#The variables F77 and CFLAGS are assumed to be set in your csh or bash
#preferences. Otherwise set them here
#F77=
#CFLAGS=
PROGRAM=fvm

OBJS=areas.o \
     bound_flux.o \
     clcd.o \
     create_edge.o \
     convert.o \
     check_positivity.o \
     divergence.o \
     divergence_hllc.o \
     divergence_kfvs.o \
     divergence_roe.o \
     dtlength.o \
     el_surr_point.o \
     el_surr_edge.o \
     finalize.o \
     flength.o \
     fvm.o \
     fvm_residual.o \
     geometric.o \
     green_gauss.o \
     initialize.o \
     isocont.o \
     math.o \
     mayavi.o \
	  misc.o \
     prep_gnuplot.o \
     pt_surr_pt.o \
     read_input.o \
     read_grid.o \
     recon.o \
     residue.o \
     result.o \
     sa_model.o \
     save_old.o \
     time_step.o \
     tri_orient.o \
     update.o \
     vigie.o \
     viscosity.o \
     wall_dist.o \
     write_grid.o

.f.o:
	$(F77) $(CFLAGS) -c $*.f

$(PROGRAM) :$(OBJS)
	              $(F77) -o $(PROGRAM) $(OBJS)

clean:;    rm -f $(OBJS) core
