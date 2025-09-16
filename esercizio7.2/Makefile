CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS = sign.o Integral.o

LIBS := `root-config --libs`
INCS := `root-config --cflags`

esercizio7.2: esercizio7.2.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c $<

sign.o: sign.h sign.cpp
Integral.o: Integral.h Integral.cpp funzioni.h sign.h

clean: 
	@rm -f esercizio7.2
	@rm -rf *.dSYM
	@rm -f *.o