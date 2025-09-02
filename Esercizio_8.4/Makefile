CXXFLAGS = -g3 -Wall --pedantic -std=c++11

LIBS := `root-config --libs`
INCS := `root-config --cflags`

myLIBS= VectFunction.o sign.o

esercizio8.4: esercizio8.4.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $<

VectFunction.o: VectFunction.cpp VectFunction.h VectAlg.h sign.h
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio8.4
	@rm -rf *.dSYM
	@rm -f *.o