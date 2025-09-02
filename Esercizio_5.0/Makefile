CXXFLAGS = -g3 -Wall --pedantic -std=c++17
myLIBS = Posizione.o

LIBS := `root-config --libs`
INCS := `root-config --cflags`

posizioneprova: posizioneprova.cpp $(myLIBS)
	g++ $(CXXFLAGS) -o $@ $^ ${INCS} ${LIBS}

%.o: %.cpp %.h
	g++ $(CXXFLAGS) -c -o $@ $<

posizioneprova: posizioneprova.cpp Posizione.o
Posizione.o: Posizione.h Posizione.cpp

clean: 
	@rm -rf *.dSYM
	@rm -f posizioneprova
	@rm -f *.o