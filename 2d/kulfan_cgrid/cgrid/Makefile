###############################################################################
#
#  Makefile for CGRID
#  ==================
#
#  (c) J. Blazek, Feb. 14, 1995
#  Modified: May 10, 2004
#
###############################################################################
#
#  EXEC = path to executable
#  O    = extension of object files (Unix, Linux, g77: .o / DOS: .obj)
#  RM   = command to delete files
#  LIB  = library
##
include ../Makefile.in

EXEC = cgrid
O = .o
RM = rm -f
LIB = ../srccom/libgrid.a

#### PC - GNU Fortran 77
#
LD = $(F77)
FC = $(F77)
#FFLAGS = -O3 -ff77 -fpedantic -Wunused -malign-double -c
FFLAGS = -c
LDFLAGS = -o

#### Mac OS X - f2c
#
# LD = f77
# FC = f77
# FFLAGS = -O -u -trapuv -c
# LDFLAGS = -O -trapuv -o

#### IBM - RS6000
#
# LD = xlf
# FC = xlf
# FFLAGS = -O3 -c
# LDFLAGS = -o

#### SGI - Octane
#
# LD = f77
# FC = f77
# FFLAGS = -64 -mips4 -O2 -c
# LDFLAGS = -64 -mips4 -O2 -o

#### SUN Ultra
#
# LD = f77
# FC = f77
# FFLAGS = -fast -xarch=v8plusa -c
# LDFLAGS = -fast -xarch=v8plusa -o

SOURCES = cgrid.f decas.f

OBJECTS = ${SOURCES:.f=$(O)}
.SUFFIXES: .f $(O)

.f$(O):
	$(FC) $(FFLAGS) $*.f

$(EXEC): $(OBJECTS) $(LIB)
	$(LD) $(LDFLAGS) $(EXEC) $(OBJECTS) $(LIB)

clean:
	$(RM) $(OBJECTS) $(EXEC)
