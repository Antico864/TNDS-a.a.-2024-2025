CXXFLAGS = -g3 -Wall --pedantic -std=c++11
myLIBS=  Solutore.o sign.o

esercizio6.2: esercizio6.2.cpp ${myLIBS}
	g++ $(CXXFLAGS) -Wall -o $@ $^

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -Wall -c $<

Solutore.o: Solutore.h Solutore.cpp sign.h
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio6.2
	@rm -rf *.dSYM
	@rm -f *.o
	@rm -f *.h.gch