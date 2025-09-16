INCS := `root-config --cflags`
LIBS := `root-config --libs`
CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS= Functions.o RandomGen.o

esercizio10.1: esercizio10.1.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $<

Functions.o: Functions.cpp Functions.hpp
RandomGen.o: RandomGen.hpp RandomGen.cpp Functions.hpp

clean: 
	@rm -f esercizio10.1
	@rm -rf *.dSYM
	@rm -f *.o