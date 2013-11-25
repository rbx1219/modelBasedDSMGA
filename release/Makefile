
CXX = g++
#CXXFLAGS = -O0 -g -Wall
#CXXFLAGS = -O0 -pg -Wall -march=nocona 
CXXFLAGS = -O3 -Wall
INCLUDE = 
TLIB = 

#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .cc and .o files

.SUFFIXES: .cpp

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

.cpp :
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ -lm $(TLIB) -lg++

#-----File Dependencies----------------------

SRC = $(SRC1) $(SRC2)

SRC1 = fastcounting.cpp bbi.cpp dsmclusteringchromosome.cpp global.cpp bitwisedistance.cpp dsmga.cpp dsmgaMain.cpp mt19937ar.cpp chromosome.cpp myrand.cpp tablelookup.cpp

SRC2 = fastcounting.cpp bisection.cpp bbi.cpp dsmclusteringchromosome.cpp global.cpp bitwisedistance.cpp dsmga.cpp mt19937ar.cpp chromosome.cpp myrand.cpp tablelookup.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

OBJ1 = $(addsuffix .o, $(basename $(SRC1)))
OBJ2 = $(addsuffix .o, $(basename $(SRC2)))

all: DSMGA bisection


DSMGA: $(OBJ1)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ1)

bisection: $(OBJ2) DSMGA
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ2)

#-----Other stuff----------------------------
depend:
	makedepend -Y. $(SRC)

clean:
	rm -f $(OBJ)

# DO NOT DELETE

fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h tablelookup.h
fastcounting.o: fastcounting.h
bbi.o: bbi.h
dsmclusteringchromosome.o: dsmclusteringchromosome.h chromosome.h global.h
dsmclusteringchromosome.o: myrand.h mt19937ar.h bitwisedistance.h
dsmclusteringchromosome.o: tablelookup.h bbi.h twodarray.h statistics.h
dsmclusteringchromosome.o: triMatrix.h fastcounting.h
global.o: global.h myrand.h mt19937ar.h bitwisedistance.h tablelookup.h
global.o: statistics.h
bitwisedistance.o: bitwisedistance.h
dsmga.o: dsmga.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga.o: tablelookup.h chromosome.h bbi.h statistics.h
dsmga.o: dsmclusteringchromosome.h twodarray.h triMatrix.h
dsmgaMain.o: statistics.h dsmga.h global.h myrand.h mt19937ar.h
dsmgaMain.o: bitwisedistance.h tablelookup.h chromosome.h bbi.h
dsmgaMain.o: dsmclusteringchromosome.h twodarray.h triMatrix.h
chromosome.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
chromosome.o: tablelookup.h
myrand.o: myrand.h mt19937ar.h
tablelookup.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
tablelookup.o: tablelookup.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h tablelookup.h
fastcounting.o: fastcounting.h
bisection.o: statistics.h dsmga.h global.h myrand.h mt19937ar.h
bisection.o: bitwisedistance.h tablelookup.h chromosome.h bbi.h
bisection.o: dsmclusteringchromosome.h twodarray.h triMatrix.h
bbi.o: bbi.h
dsmclusteringchromosome.o: dsmclusteringchromosome.h chromosome.h global.h
dsmclusteringchromosome.o: myrand.h mt19937ar.h bitwisedistance.h
dsmclusteringchromosome.o: tablelookup.h bbi.h twodarray.h statistics.h
dsmclusteringchromosome.o: triMatrix.h fastcounting.h
global.o: global.h myrand.h mt19937ar.h bitwisedistance.h tablelookup.h
global.o: statistics.h
bitwisedistance.o: bitwisedistance.h
dsmga.o: dsmga.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga.o: tablelookup.h chromosome.h bbi.h statistics.h
dsmga.o: dsmclusteringchromosome.h twodarray.h triMatrix.h
chromosome.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
chromosome.o: tablelookup.h
myrand.o: myrand.h mt19937ar.h
tablelookup.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
tablelookup.o: tablelookup.h
