#include <iomanip>
#include <iostream>

#include "FunzioneBase.h"
#include "Solutore.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 4) {
    cout << "Uso del programma: " << argv[0] << " <a> " << " <b> "
         << " <precisione>" << endl;
    cout << "<precisione> è l'ordine di cifra che si desidera visualizzare nel "
            "valore di ascissa dello zero. "
         << endl;
    exit(2);
  }

  int xm = atof(argv[1]);
  int xM = atof(argv[2]);
  double epsilon = atof(argv[3]);
  int cifre_sign = (int)(-log10(epsilon));

  Parabola mypar(3, 5, -2);

  Bisezione bysector;
  double zero = bysector.CercaZeriRef(xm, xM, mypar, epsilon, 500);
  cout << fixed;
  cout << "x_0 = " << setprecision(cifre_sign) << zero << " +- " << epsilon
       << endl;

  return 0;
}
