INCS := `root-config --cflags`
LIBS := `root-config --libs`

CXXFLAGS = -g3 -Wall --pedantic -std=c++11

myLIBS = Functions.o RandomGen.o Integral.o

EXECUTABLES = calcsaveavg calcsavehom analysisavg analysishom

calcsaveavg: calcsaveavg.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}
calcsavehom: calcsavehom.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}
analysisavg: analysisavg.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}
analysishom: analysishom.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp
	g++ $(CXXFLAGS) -c -o $@ $<

Functions.o: Functions.cpp Functions.hpp
RandomGen.o: RandomGen.cpp RandomGen.hpp
Integral.o: Integral.cpp Integral.hpp

calcsaveavg.o: calcsaveavg.cpp Functions.hpp Integral.hpp RandomGen.hpp
calcsavehom.o: calcsavehom.cpp Functions.hpp Integral.hpp RandomGen.hpp
analysisavg.o: analysisavg.cpp Functions.hpp Integral.hpp RandomGen.hpp
analysishom.o: analysishom.cpp Functions.hpp Integral.hpp RandomGen.hpp

clean:
	@rm -rf *.dSYM
	@rm -f calcsaveavg
	@rm -f calcsavehom
	@rm -f analysisavg
	@rm -f analysishom
	@rm -f *.o
	@rm -f *.txt