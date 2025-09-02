#include <iostream>

#include "vectalg.h"

using namespace std;

int main(int argc, const char **argv) {
  double F_1, F_2, F_3, G_1, G_2, G_3;

  cout << "Insert components of force n° 1 [N]: " << endl;
  cin >> F_1;
  cin >> F_2;
  cin >> F_3;
  cout << "Insert components of force n° 2 [N]: " << endl;
  cin >> G_1;
  cin >> G_2;
  cin >> G_3;

  const vector<double> F{F_1, F_2, F_3};
  const vector<double> G{G_1, G_2, G_3};
  //  Metterli const assicura che venga usato l'operator+ che non sovrascrive
  //  il vettore, permettendoci di fare il prodotto scalare anche
  //  dopo la somma.

  vector<double> sum = F + G;
  double product = F * G;

  cout << "sum: " << endl;
  Print(sum);

  cout << "Scalar product F*G = " << product << endl;

  return 0;
}