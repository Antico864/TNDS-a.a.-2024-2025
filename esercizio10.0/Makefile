CXXFLAGS = -g3 -Wall --pedantic -std=c++11
INCS := `root-config --cflags`
LIBS := `root-config --libs`

myLIBS= Functions.o RandomGen.o

esercizio10.0: esercizio10.0.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.hpp
	g++ $(CXXFLAGS) -c -o $@ $<

Functions.o: functions.hpp Functions.cpp
RandomGen.o: RandomGen.hpp RandomGen.cpp Functions.hpp

clean: 
	@rm -f esercizio10.0
	@rm -rf *.dSYM
	@rm -f *.o