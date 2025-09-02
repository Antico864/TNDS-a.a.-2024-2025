CXXFLAGS = -g3 -Wall --pedantic -std=c++11

LIBS := `root-config --libs`
INCS := `root-config --cflags`

myLIBS = Integral.o sign.o

esercizio7.0: esercizio7.0.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c $<

Integral.o: Integral.h Integral.cpp sign.h funzioni.h
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio7.0
	@rm -rf *.dSYM
	@rm -f *.o