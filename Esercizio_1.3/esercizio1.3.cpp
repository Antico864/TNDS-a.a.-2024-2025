#include "funzioni.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {

  if (argc < 3) {
    cout << "Uso del programma : " << argv[0] << " <n_data> <filename> "
         << endl;
    return -1;
  }
  int ndata = atoi(argv[1]);
  char *Filename = argv[2];
  double *data = readDataFromFile(Filename, ndata);

  for (int k = 0; k < ndata; k++) {
    cout << data[k] << endl;
  }

  double media = doMedia(data, ndata);
  cout << "Media = " << media << endl;

  double varianza = doVarianza(data, ndata);
  cout << "Varianza = " << varianza << endl;

  double *vcopy = new double[ndata];

  for (int i = 0; i < ndata; i++) {
    vcopy[i] = data[i];
  }

  double mediana = doMediana(vcopy, ndata);

  Print("19410.txt", vcopy, ndata);
  Print(vcopy, ndata);
  cout << "Mediana = " << mediana << endl;

  Print("19410.txt", "\nMediana = " + to_string(mediana));
  Print("19410.txt", "Media = " + to_string(media));
  Print("19410.txt", "Varianza = " + to_string(varianza));

  delete[] vcopy;

  return 0;
}
