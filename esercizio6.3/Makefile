CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS = Solutore.o sign.o

esercizio6.3: esercizio6.3.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c $<

Solutore.o: Solutore.h Solutore.cpp sign.h FunzioneBase.h
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio6.3
	@rm -rf *.dSYM
	@rm -f *.o