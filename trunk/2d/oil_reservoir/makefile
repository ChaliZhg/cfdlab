CXX    = c++
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
	rm -f *.o $(TARGET) *.vtk
