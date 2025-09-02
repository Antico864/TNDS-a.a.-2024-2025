#include <cstdlib>
#include <iostream>
#include <vector>

#include "funzioni.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 2) {
    cout << "Uso del programma : " << argv[0] << " <filename> " << endl;
    return -1;
  }

  vector<double> a = Read<double>(argv[1]);

  MergeSort(a, 0, a.size()-1);

  double media = CalcAvg<double>(a);
  double varianza = CalcVar<double>(a);
  double mediana = CalcMed<double>(a);

  Print("19410.txt", a);
  Print("19410.txt", "Media = ", media);
  Print("19410.txt", "Varianza = ", varianza);
  Print("19410.txt", "Mediana = ", mediana);

  Print(a);

  cout << "Media = " << media << endl;
  cout << "Varianza = " << varianza << endl;
  cout << "Mediana = " << mediana << endl;

  return 0;
}
