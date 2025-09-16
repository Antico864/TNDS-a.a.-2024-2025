#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "FunzioneBase.h"
#include "Solutore.h"
#include "sign.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 2) {
    cout << "Use of program: " << "./esercizio6.3 " << "<precision>" << endl;
    exit(33);
  }

  double prec = atof(argv[1]);
  if (prec > 1E-6)
    prec = 1E-6;
  double epsilon = prec;
  int cifre_sign = (int)(-log10(prec));

  vector<double> v(20);

  Tan_men_x myfunc;

  for (int n = 0; n < 20; n++) {

    double mx = (((double)n))*M_PI;
    double Mx = (((double)n) + 0.5 - 0.8 * epsilon) * M_PI;
    // Potremmo fare 5 - 0.1 * epsilon, per essere sicuri di evitare il rounding...?

    Bisezione bysector;

    double zero = bysector.CercaZeriRef(mx, Mx, myfunc, epsilon, cifre_sign);

    v[n] = zero;

    cout << fixed;
    cout << "Zero n°" << n + 1 << ": " << setprecision(cifre_sign) << v[n]
         << " +- " << epsilon << endl;
  }

  return 0;
}
