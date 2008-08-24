#
# fortran flags
#
FC       = gfortran
FFLAGS   = ${CFLAGS} -fdefault-real-8
FFLAGS   = ${CFLAGS}

TARGETS = grid solve solveg

#Some colorful printing
RED=\033[30;43m
GREEN=\033[30;42m
BLACK=\033[0m
PRINT := $(shell \
	[ -z `echo -e` ] && PRINT='echo -e' || PRINT='echo';\
	echo $$PRINT;\
	)

#
# program dependencies
#

grid  =  grid.o

solve =  solve.o

solveg =  solveg.o

#
# programs
#

ALL:  $(TARGETS)

grid: ${grid}
		${FC} -o grid ${grid}

solve:  ${solve}
		${FC} -o solve ${solve} -lgrafic -lX11

solveg:  ${solveg}
		${FC} -o solveg ${solveg} -lgrafic -lX11

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.F.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

solve.o: solve.F
		${FC} ${FFLAGS} -c -o solve.o solve.F

solveg.o: solve.F
		${FC} ${FFLAGS} -D_GRAFIC -c -o solveg.o solve.F

##############################################################################
# clean things up
##############################################################################

clean:	
	${RM} -f *.o $(TARGETS)
