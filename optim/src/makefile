#
# fortran flags
#
CC       = gcc
FC       = gfortran
FFLAGS   = ${CFLAGS}

INDENT   = indent -br -brs -npcs -npsl -nsaf -nsai -nsaw -cli3 -l78 -i3 -ts0

TARGETS = optim

#
# program dependencies
#

SRCC  = $(wildcard *.c)
OBJC  = $(patsubst %.c,%.o,$(SRCC))
optim = $(OBJC)

#
# programs
#

ALL:  $(TARGETS)

optim:  ${optim}
		${CC} -o optim ${optim} -lm

##############################################################################
# transformation rules
##############################################################################

.f.o:     ; ${FC} ${FFLAGS} -c $*.f
.c.o:     ; ${CC} ${CFLAGS} -c $*.c

##############################################################################
indent:
	${INDENT} *.c
	${RM} *.c~
##############################################################################
# clean things up
##############################################################################

clean:
	${RM} *.o $(TARGETS)
