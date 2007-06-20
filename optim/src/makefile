#
# fortran flags
#
CC       = gcc
CC       = mpicc
FC       = gfortran
FFLAGS   = ${CFLAGS}

INDENT   = indent -br -brs -npcs -npsl -nsaf -nsai -nsaw -cli3 -l78 -i3 -ts0

TARGETS = optim optim_mpi

#
# program dependencies
#

SRCC  = $(wildcard *.c)
OBJC  = $(patsubst %.c,%.o,$(SRCC))
optim = $(OBJC) mpi/nompi.o
optim_mpi = $(OBJC) mpi/mpi.o

#
# programs
#

ALL:  $(TARGETS)

optim:  ${optim}
		${CC} -o optim ${optim} -lm

optim_mpi:  ${optim_mpi}
		${CC} -o optim_mpi ${optim_mpi} -lm

##############################################################################
# transformation rules
##############################################################################

.c.o:
		${CC} ${CFLAGS} -Wall \
	      -c $*.c

#.f.o:     ; ${FC} ${FFLAGS} -c $*.f
#.c.o:     ; ${CC} ${CFLAGS} -Wall -c $*.c

mpi/%.o:  mpi/%.c
		cd mpi; ${CC} -c $*.c

##############################################################################
indent:
	${INDENT} *.c
	${RM} *.c~
##############################################################################
# clean things up
##############################################################################

clean:
	${RM} *.o $(TARGETS)
