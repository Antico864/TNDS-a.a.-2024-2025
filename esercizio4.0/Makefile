CXXFLAGS = -g3 -Wall --pedantic -std=c++17

LIBS :=`root-config --libs`
INCS :=`root-config --cflags`

esercizio4.0: esercizio4.0.cpp funzioni.h
	g++ -o esercizio4.0 esercizio4.0.cpp ${INCS} ${LIBS} $(CXXFLAGS)

clean: 
	@rm -f esercizio4.0
	@rm -rf *.dSYM
	@rm -f Trend.pdf