#include "Vettore.h"
#include "funzioni.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {

  if (argc < 2) {
    cout << "Use of program: " << argv[0] << " <n_data> <filename> " << endl;
    return -1;
  }
  int ndata = atoi(argv[1]);
  char *Filename = argv[2];
  Vettore vett0 = Read(ndata, Filename);

  Print(vett0);

  double media = CalcAvg(vett0);
  cout << "Media = " << media << endl;

  double varianza = CalcVar(vett0);
  cout << "Varianza = " << varianza << endl;

  double mediana = CalcMed(vett0);

  Print("19410.txt", vett0);
  Print(vett0);

  cout << "Mediana = " << mediana << endl;

  Print("19410.txt", "\nMedia = " + to_string(media));
  Print("19410.txt", "Varianza = " + to_string(varianza));
  Print("19410.txt", "Mediana = " + to_string(mediana));

  return 0;
}
