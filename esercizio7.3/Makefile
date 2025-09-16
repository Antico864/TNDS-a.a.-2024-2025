CXXFLAGS := -g3 -Wall --pedantic -std=c++11

myLIBS:= Integral.o sign.o

INCS := `root-config --cflags`
LIBS := `root-config --libs`

esercizio7.3: esercizio7.3.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c $<

Integral.o: Integral.h Integral.cpp funzioni.h sign.h
funzioni.o: funzioni.h funzioni.cpp
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio7.3
	@rm -rf *.dSYM
	@rm -f *.o