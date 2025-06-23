SRC= ./include
CXXFLAGS = -fPIC -w -fopenmp -Wall -std=c++0x -I$(SRC) -O2 #-DDBG

OBJECTS= Earth.o Table.o Utilities.o

all: Earth.o Table.o Utilities.o Propagator Simu_elost

Earth.o: $(SRC)/Earth.cc
	$(CXX) -c $(SRC)/Earth.cc -o Earth.o $(CXXFLAGS)
Table.o: $(SRC)/Table.cc
	$(CXX) -c $(SRC)/Table.cc -o Table.o $(CXXFLAGS)
Utilities.o: $(SRC)/Utilities.cxx
	$(CXX) -c $(SRC)/Utilities.cxx -o Utilities.o $(CXXFLAGS)
Propagator: propagator.cxx $(OBJECTS)
	$(CXX) propagator.cxx -o Propagator $(CXXFLAGS) $(OBJECTS)
Simu_elost: Simu_elost.cxx $(OBJECTS)
	$(CXX) Simu_elost.cxx -o Simu_elost $(CXXFLAGS) $(OBJECTS)
clean:
	rm *o Simu_elost

