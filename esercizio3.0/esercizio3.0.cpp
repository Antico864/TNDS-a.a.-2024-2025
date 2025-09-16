#include <iostream>

#include "Vettore.h"
#include "funzioni.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 2) {
    cout << "Use of program : " << argv[0] << " <n_data> <filename> " << endl;
    return -1;
  }

  int ndata = atoi(argv[1]);
  char *filename = argv[2];
  Vettore<double> vett0 = Read<double>(ndata, filename);

  Print(vett0);

  double media = CalcAvg(vett0);
  double varianza = CalcVar(vett0);

  cout << "Media = " << media << endl;
  cout << "Varianza = " << varianza << endl;

  Vettore<double> A(vett0);

  double mediana = CalcMed(A);
  Print(A);
  cout << "Mediana = " << mediana << endl;

  Print("19410.txt", A);
  Print("19410.txt", "\nMedia = " + to_string(media));
  Print("19410.txt", "Varianza = " + to_string(varianza));
  Print("19410.txt", "Mediana = " + to_string(mediana));

  return 0;
}
