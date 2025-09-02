CXXFLAGS = -g3 -Wall --pedantic -std=c++11

LIBS :=`root-config --libs`
INCS :=`root-config --cflags`

myLIBS = Integral.o sign.o

esercizio7.1: esercizio7.1.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c $<

sign.o: sign.cpp sign.h
Integral.o: Integral.cpp Integral.h funzioni.h sign.h

clean: 
	@rm -f esercizio7.1
	@rm -rf *.dSYM
	@rm -f *.o