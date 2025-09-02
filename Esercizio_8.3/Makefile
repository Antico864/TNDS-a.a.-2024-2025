CXXFLAGS = -g3 -Wall --pedantic -std=c++11

LIBS := `root-config --libs`
INCS := `root-config --cflags`

myLIBS= VectFunction.o

esercizio8.3: esercizio8.3.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $<

VectFunction.o: VectFunction.cpp VectFunction.h VectAlg.h

clean: 
	@rm -f esercizio8.3
	@rm -rf *.dSYM
	@rm -f *.o