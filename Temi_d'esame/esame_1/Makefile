INCS := `root-config --cflags`
LIBS := `root-config --libs`

CXXFLAGS = -g3 -ggdb3 -Wall --pedantic -std=c++11

MyLibs = functions.o Integral.o RandomGen.o

compito: compito.cpp $(MyLibs)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.hpp
	g++ $(CXXFLAGS) -c $<

## DIPENDENZE!!!!!

clean:
	@rm -f *.o
	@rm -f compito
	@rm -rf *.dSYM