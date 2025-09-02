#include "Integral.hpp"
#include "RandomGen.hpp"

using namespace std;

//  Nota bene: fmax non è necessario qui, ma dobbiamo mantenerlo
//  per rispettare il metodo ereditato. Anche in eq diff succedeva
//  qualcosa di simile, con la dipendenza dal tempo negli oscillatori.
double AvgIntegrator::Integrate(RandomGen &ran, const FunzioneBase &f,
                                double inf, double sup, unsigned int points,
                                double fmax) {
  double avg = 0.;
  double delta = 0.;
  for (unsigned int i = 0; i < points; i++) {
    double x = ran.Unif(inf, sup);
    delta = f.Eval(x) - avg;
    avg += delta / (i + 1);
  }
  return (sup - inf)*avg;
}

double HoMIntegrator::Integrate(RandomGen &ran, const FunzioneBase &f,
                                double inf, double sup, unsigned int points,
                                double fmax) {
  unsigned int Hcounter = 0;
  for (unsigned int i = 0; i < points; i++) {
    double x = ran.Unif(inf, sup);
    double y = ran.Unif(0., fmax);
    if (y < f.Eval(x))
      Hcounter++;
  }
  return (sup - inf) * fmax *
         (static_cast<double>(Hcounter) / static_cast<double>(points));
}
