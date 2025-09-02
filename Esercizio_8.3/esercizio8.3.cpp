#include "VectAlg.h"
#include "VectFunction.h"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "Use of program: " << argv[0] << " <resolution_step_size>" << endl;
    exit(-1);
  };

  cout << " " << endl;
  cout << "//=========================================================" << endl;
  cout << "//  The program produces a graph of the period of oscillation"
       << endl;
  cout << "//  for the motion of a simple pendulum" << endl;
  cout << "//  using the Runge-Kutta method. " << endl;
  cout << "//=========================================================" << endl;
  cout << " " << endl;

  TApplication app("app", 0, 0);

  double h = atof(argv[1]);
  double theta0;
  double length = 0.5;

  Pendulum mypen(length);

  Runge_Kutta runge_kutta;

  TGraph per;

  // Calcolo del periodo:
  for (int i = 0; i < 30; i++) {
    theta0 = 0.1 * (i + 1); // Incremento theta-{0} di 0.1;
    double v = 0.;          // Resetto i parametri del sistema;
    double t = 0.;
    vector<double> x{-theta0, v};
    while (x[1] >= 0.) {
      v = x[1]; // Imposta a dx/dt il valore della velocità;
      x = runge_kutta.Step(t, x, h,
                           mypen); //  Calcola la soluzione fino al tempo t;
      t = t + h; //  Aggiorna il valore di t con h passato da riga di comando;
    }
    // x[1] ha cambiato segno!
    t = t - v * h / (x[1] - v); // Interpolazione lineare sul tempo;
    double period = 2 * t;      // obv
    per.SetPoint(i, theta0, period);
  }

  TCanvas c("Armonic oscillator solution", "Armonic oscillator solution");

  c.cd();

  per.SetMarkerSize(0.5);
  per.SetMarkerStyle(20);
  per.SetMarkerColor(kBlue);
  per.SetFillColor(0);
  per.SetLineColor(kBlack);

  string str = "Period of oscillation, h = " + convert(h);

  per.SetTitle(str.c_str());
  per.GetXaxis()->SetTitle("#theta_{0} [rad]");
  per.GetYaxis()->SetTitle("T [s]");

  per.Draw("APL");

  c.Update();

  app.Run();
}