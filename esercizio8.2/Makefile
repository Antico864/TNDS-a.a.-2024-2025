CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS= VectFunction.o sign.o

LIBS :=`root-config --libs`
INCS :=`root-config --cflags`

esercizio8.2: esercizio8.2.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${LIBS} ${INCS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $< 

VectFunction.o: VectFunction.h VectFunction.cpp VectAlg.h
sign.o: sign.h sign.cpp

clean: 
	@rm -f esercizio8.2
	@rm -rf *.dSYM
	@rm -f *.o