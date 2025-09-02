CXXFLAGS = -g3 -Wall --pedantic -std=c++11

LIBS := `root-config --libs`
INCS := `root-config --cflags`

myLIBS = Posizione.o PuntoMateriale.o CampoVettoriale.o Particella.o

esercizio5.3: esercizio5.3.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}
%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $<

Posizione.o: Posizione.h Posizione.cpp
PuntoMateriale.o: PuntoMateriale.h PuntoMateriale.cpp
CampoVettoriale.o: CampoVettoriale.h CampoVettoriale.cpp
Particella.o: Particella.h Particella.cpp

clean: 
	@rm -f esercizio5.3
	@rm -f *.o
	@rm -rf *.dSYM