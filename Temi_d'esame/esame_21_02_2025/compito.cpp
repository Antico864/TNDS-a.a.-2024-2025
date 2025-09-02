#include "myLibrary.hpp"

using namespace std;

int main() {
  const unsigned int N = 10;
  double sigmanec = 1E-10;
  // Svolgimento punto 1)
  cout << "Svolgimento punto 1)\n" << endl;
  functiontest myf;
  Simpson mysimp(1., 2.);
  double I = mysimp.Integrate(N, myf);
  cout << "I = " << I << endl;
  cout << "\n\n" << endl;

  // Svolgimento punto 2)
  cout << "Svolgimento punto 2)\n" << endl;
  double errorsimp = mysimp.Error(10, myf);
  cout << "L'errore commesso integrando con Simpson a 10 intervalli è "
       << errorsimp << endl;
  unsigned int N_simp = mysimp.Nnec(sigmanec, N, errorsimp);
  cout << "Per ottenere un errore di " << sigmanec
       << " con Simpson è necessario dividere l'intervallo di integrazione in "
       << N_simp << " sottointervalli. " << endl;
  cout << "\n\n" << endl;

  // Svolgimento punto 3)
  cout << "Svolgimento punto 3)\n" << endl;
  vector<double> state = {0.};
  vectftest diffeq(myf);
  Runge_Kutta myrk;
  double t;
  const double h = 1. / N;
  for (t = 1.; t <= 2.; t += h) {
    state = myrk.Step(t, state, h, diffeq);
  }
  cout << "I = " << state[0] << endl;
  cout << "\n\n" << endl;

  // Svolgimento punto 4)
  cout << "Svolgimento punto 4)\n" << endl;
  state = {0.};
  double errorrk = myrk.errorcomponent(1., 2., state, h, diffeq);
  cout << "L'errore commesso integrando con Runge-Kutta a 10 step è " << errorrk
       << endl;
  unsigned int Nrk = myrk.Nnecess(N, errorrk, sigmanec, diffeq, 1., 2.);
  cout << "Per ottenere un errore di " << sigmanec << " con Runge-Kutta sono "
       << Nrk << " passi. " << endl;
  cout << "\n\n" << endl;

  // Svolgimento punto 5)
  cout << "Svolgimento punto 5)\n" << endl;
  RandomGen myrg(1);
  AvgIntegrator myavg(1);
  double Iavg = myavg.Integrate(myrg, myf, 1., 2., 10, myf.Eval(1.));
  cout << "Metodo della media montecarlo: I = " << Iavg << endl;
  double erroravg = myavg.GetError();
  cout << "L'errore commesso integrando con metodo della media montecarlo a 10 "
          "step è "
       << erroravg << endl;
  cout << "\n\n" << endl;

  // Svolgimento punto 6)
  cout << "Svolgimento punto 6)\n" << endl;
  double Navg = Nnecavg(N, erroravg, sigmanec);
  cout << "Per ottenere un errore di " << sigmanec
       << " con metodo della media montecarlo è necessario estrarre " << Navg
       << " punti. " << endl;
  cout << "\n\n" << endl;
  // Va bene anche double, spero, perché è enorme.

  // Svolgimento punto 7)
  cout << "Svolgimento punto 7)\n" << endl;
  RandomGen myrgimp(1);
  AvgIntegrator myavgimp(1);
  double Iimpsample = myavgimp.ImportanceSample(myrgimp, myf, 1., 2., 10);
  cout << "Metodo della media importance sampling: I = " << Iimpsample << endl;
  double errorimpsample = myavgimp.GetError();
  cout << "L'errore commesso integrando con metodo dell'importance sampling a "
          "10 step è "
       << errorimpsample << endl;
  cout << "L'errore rispetto al metodo della media è migliorato di "
       << erroravg - errorimpsample << endl;
  cout << "\n\n" << endl;
  return 0;
}
