INCS := `root-config --cflags`
LIBS := `root-config --libs`

CXXFLAGS = -g3 -ggdb3 -Wall --pedantic -std=c++20

compito: compito.o
	g++ $(CXXFLAGS) -o compito compito.o ${INCS} ${LIBS}

compito.o: compito.cpp myLibrary.hpp
	g++ $(CXXFLAGS) -c compito.cpp -o compito.o ${INCS}

clean:
	@rm -f *.o
	@rm -f compito
	@rm -rf *.dSYM