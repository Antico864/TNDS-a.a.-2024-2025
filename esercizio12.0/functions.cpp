#include "functions.hpp"

using namespace std;

void PrintWithUncertainty(double value, double uncertainty) {
  // Arrotondare il valore e l'incertezza alla precisione richiesta
  int precision = ceil(-log10(uncertainty)); // Numero di cifre decimali
  cout << fixed << setprecision(precision) << value << " ± " << uncertainty
       << endl;
}

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace std;

string PrintWithUncertaintyScient(double value, double uncertainty) {
  if (uncertainty == 0)
    return "Value without uncertainty";

  // Ordine di grandezza dell'incertezza:
  int order_of_unc = floor(log10(uncertainty));
  // Ordine di grandezza del valore:
  int order_of_val = floor(log10(fabs(value)));

  // Numero di cifre significative per il valore:
  unsigned int sign_fig_val = abs(order_of_val - order_of_unc);

  stringstream result;
  result << scientific << setprecision(sign_fig_val) << value;
  result << " ± " << scientific << setprecision(0) << uncertainty;

  return result.str();
}
