CXXFLAGS = -g3 -Wall --pedantic -std=c++11

esercizio2.2: esercizio2.2.o Vettore.o funzioni.o
	g++ esercizio2.2.o Vettore.o funzioni.o -o esercizio2.2 $(CXXFLAGS)
esercizio2.2.o: esercizio2.2.cpp Vettore.h funzioni.h
	g++ -c esercizio2.2.cpp -o esercizio2.2.o $(CXXFLAGS)
Vettore.o: Vettore.cpp Vettore.h
	g++ -c Vettore.cpp -o Vettore.o $(CXXFLAGS)
funzioni.o: funzioni.cpp funzioni.h
	g++ -c funzioni.cpp -o funzioni.o $(CXXFLAGS)

clean: 
	@rm -r esercizio2.2
	@rm -f *.o
	@rm -f 19410.txt