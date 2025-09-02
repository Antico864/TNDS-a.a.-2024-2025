CXXFLAGS = -g3 -Wall --pedantic


esercizio8.0: esercizio8.0.cpp VectAlg.h
	g++ -Wall -o $@ $^$(CXXFLAGS)

%.o: %.cpp %.h
	g++ -Wall -c $< $(CXXFLAGS)


clean: 
	@rm -f esercizio8.0
	@rm -rf *.dSYM
	@rm -f *.o