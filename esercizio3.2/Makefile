CXXFLAGS = -g3 -Wall --pedantic -std=c++17

LIBS :=`root-config --libs`
INCS :=`root-config --cflags`

esercizio3.2: esercizio3.2.cpp funzioni.h 
	g++ -o esercizio3.2 esercizio3.2.cpp ${INCS} $(LIBS) $(CXXFLAGS)

clean: 
	@rm -f esercizio3.2
	@rm -rf esercizio3.2.dSYM