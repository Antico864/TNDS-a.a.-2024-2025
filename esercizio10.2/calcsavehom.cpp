#include "Functions.hpp"
#include "Integral.hpp"
#include "RandomGen.hpp"

using namespace std;

int main() {

  const xsinx f(1., 1., 0.);
  unsigned int N = 100;
  unsigned int i = 5;
  vector<double> integ(10000);
  HoMIntegrator hom1(1);

  while (N <= 100000) {

    RandomGen rand(1);

    for (unsigned int k = 0; k < 10000; k++) {
      integ[k] = hom1.Integrate(rand, f, 0., M_PI / 2, N, M_PI / 2);
    }

    // Stampa i dati del vector nel corretto file di testo:
    const string output = "output" + to_string(N) + "hom.txt";
    Print(output, integ);

    //  Verifica dei valori medi:
    cout << "Value of integral, N = " << N << ": " << Avg(integ) << endl;

    //  Aggiorna N:
    N *= i;
    if (i % 2 != 0) {
      i = 2;
    } else
      i = 5;
  }

  return 0;
}
