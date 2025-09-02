CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS= VectAlg.o VectFunction.o sign.o

LIBS :=`root-config --libs`
INCS :=`root-config --cflags`

esercizio8.1: esercizio8.1.cpp ${myLIBS}
	g++ $(CXXFLAGS) -o $@ $^ ${LIBS} ${INCS} 

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $<

VectAlg.o: VectAlg.h VectAlg.cpp
VectFunction.o: VectFunction.h VectFunction.cpp sign.h
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio8.1
	@rm -rf *.dSYM
	@rm -f *.o