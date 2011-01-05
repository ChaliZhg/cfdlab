FC = ifort
CC = gcc

FFLAGS = -r8

TPN =  ${TAPENADE_HOME}/bin/tapenade -r8

TARGETS=grid flo adj ccre

FLO = flo.o bc.o residu.o source.o costfun.o read_input.o

ADJ = adj.o bc.o sourceis_ba.o costfun_bq.o residu_bq.o adStack.o adBuffer.o \
	  	residu.o read_input.o

all: $(TARGETS)

grid: grid.f
	$(FC) $(FFLAGS) -o grid grid.f

flo: ${FLO}
	$(FC) -o flo ${FLO}

adj: ${ADJ}
	$(FC) -o adj ${ADJ}

ccre: ccre.f
	$(FC) $(FFLAGS) -o ccre ccre.f

residu_bq.f: residu.f
	$(TPN) -b -head residu \
		    -vars "q res" \
	       -outvars "q res" \
			 -difffuncname "_bq" \
			 -o residu $<

residu_bq.o: residu_bq.f
	$(FC) -c $(FFLAGS) $<

costfun_bq.f: costfun.f
	$(TPN) -b -head costfun \
		    -vars "q cost" \
	       -outvars "q cost" \
			 -difffuncname "_bq" \
			 -o costfun $<

costfun_bq.o: costfun_bq.f
	$(FC) -c $(FFLAGS) $<

sourceis_ba.f: sourceis.f
	$(TPN) -b -head sourceis \
		    -vars "a res" \
	       -outvars "a res" \
			 -difffuncname "_ba" \
			 -o sourceis $<

sourceis_ba.o: sourceis_ba.f
	$(FC) -c $(FFLAGS) $<

adBuffer.o: $(TAPENADE_HOME)/ADFirstAidKit/adBuffer.f
	$(FC) -c $<

adStack.o: $(TAPENADE_HOME)/ADFirstAidKit/adStack.c
	$(CC) -c $<

.f.o:     ; ${FC} ${FFLAGS} -c $*.f

clean:
	rm -f *.o $(TARGETS)

allclean:
	rm -f *.dat *.o *_ba.f *_bq.f *_bq.h *.msg *~ $(TARGETS)
