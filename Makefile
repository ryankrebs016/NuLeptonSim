SRC= ./include
CXXFLAGS = -fPIC -w -fopenmp -Wall -std=c++0x -I$(SRC)

OBJECTS= Earth.o Table.o

all: Earth.o Table.o ARA_sec

Earth.o: $(SRC)/Earth.cc
	$(CXX) -c $(SRC)/Earth.cc -o Earth.o $(CXXFLAGS)
Table.o: $(SRC)/Table.cc
	$(CXX) -c $(SRC)/Table.cc -o Table.o $(CXXFLAGS)
ARA_sec: ARA_sec.cxx $(OBJECTS)
	$(CXX) ARA_sec.cxx -o ARA_sec $(CXXFLAGS) $(OBJECTS)
clean:
	rm *o ARA_sec

