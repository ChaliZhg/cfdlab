CXX    = g++
CFLAGS = -O2 -Wall

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,%.o,$(SRC))

TARGET = reservoir

all: $(TARGET)

reservoir: $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(CXX) $(CFLAGS) -c $*.cc

clean:
	rm -f *.o $(TARGET)

allclean:
	make clean
	rm -f *.vtk
