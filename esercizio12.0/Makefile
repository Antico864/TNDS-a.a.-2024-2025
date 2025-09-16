LIBS := `root-config --libs`
INCS := `root-config --cflags`
CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS = Prism.o RandomGen.o functions.o

esercizio12.0: esercizio12.0.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp
	g++ $(CXXFLAGS) -c -o $@ $<

Prism.o: Prism.cpp Prism.hpp
RandomGen.o: RandomGen.cpp RandomGen.hpp
functions.o: functions.cpp functions.hpp

clean:
	@rm -rf *.dSYM
	@rm -f esercizio12.0
	@rm -f *.o