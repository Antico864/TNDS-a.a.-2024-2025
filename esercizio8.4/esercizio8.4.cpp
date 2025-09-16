#include "VectFunction.h"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

int main(int argc, const char **argv) {
  if (argc < 2) {
    cerr << "Use of program: " << argv[0] << " <resolution precision> " << endl;
    exit(-1);
  };

  TApplication app("app", 0, 0);

  double h = atof(argv[1]);
  double omega0 = 10.;
  double omegaF = 9.99;
  double alpha = 1. / 30;

  ostringstream om0s;
  om0s << fixed << setprecision(1) << omega0;
  string om0 = "#omega_{0} = " + om0s.str() + " rad/s";

  ostringstream omFs;
  omFs << fixed << setprecision(1) << omegaF;
  string omF = "#omega_{F} = " + omFs.str() + " rad/s";

  ostringstream alphas;
  alphas << fixed << setprecision(3) << alpha;

  Runge_Kutta rk;
  FDOscillator FDOsc(omega0, alpha, omegaF);
  TGraph lorentzian;
  TGraph solution;

  double t = 0.;
  vector<double> state{0., 0.};

  // Faccio evolvere il sistema e riempio il grafico della soluzione:
  while (t < 300) {
    state = rk.Step(t, state, h, FDOsc);
    solution.SetPoint(30 * t, t, state[0]);
    t += h;
  }

  // Omega va da 9 a 11 con passo di 0.05;
  for (unsigned int i = 0; i <= 40; i++) {
    //  Imposto i parametri;
    omegaF = 9 + i * 0.05;
    double v = 0.;
    double t = 0.;
    double x = 0.;
    state = {x, v};

    FDOscillator myosc(omega0, alpha, omegaF);
    // Faccio evolvere il sistema fino a t = 1/alpha;
    while (t < 300) {
      state = rk.Step(t, state, h, myosc);
      t += h;
    }

    // Faccio evolvere il sistema per poco più di un semiperiodo
    // e aggiorno il valore dell'ampiezza:
    double T = 2 * M_PI / omegaF;
    double amp = fabs(state[0]);
    while (t < 300 + 1.1 * T / 2) {
      state = rk.Step(t, state, h, myosc);
      t += h;
      if (fabs(state[0]) >= amp) {
        amp = fabs(state[0]);
      }
    }

    lorentzian.SetPoint(i, omegaF, amp);

    //  E ripeti.
  }

  TCanvas csol("Solution", "Resonance solution");

  solution.SetMarkerSize(0.5);
  solution.SetMarkerStyle(20);
  solution.SetMarkerColor(kBlue);
  solution.SetFillColor(0);
  solution.SetLineColor(kBlack);

  string str = "Forced damped oscillator solution";

  solution.SetTitle(str.c_str());
  solution.GetXaxis()->SetTitle("t [s]");
  solution.GetYaxis()->SetTitle("x [m]");

  solution.Draw("APC");

  TLegend legend(0.75, 0.75, 0.9, 0.9);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry((TObject *)0, om0.c_str(), "");
  legend.AddEntry((TObject *)0, omF.c_str(), "");

  legend.Draw();

  csol.Update();
  csol.SaveAs("Images/solution.pdf");

  TCanvas c("Lorentzian curve", "Lorentzian");

  c.cd();

  lorentzian.SetMarkerSize(0.5);
  lorentzian.SetMarkerStyle(20);
  lorentzian.SetMarkerColor(kBlue);
  lorentzian.SetFillColor(0);
  lorentzian.SetLineColor(kBlack);

  string str1 = "Lorentzian curve, #omega_{0} = " + om0s.str() +
                " rad/s, #alpha = " + alphas.str() + " s^{-1}";

  lorentzian.SetTitle(str1.c_str());
  lorentzian.GetXaxis()->SetTitle("#omega_{F} [rad/s]");
  lorentzian.GetYaxis()->SetTitle("A [m]");

  lorentzian.Draw("APC");

  c.Update();

  c.SaveAs("Images/lorentz.pdf");

  app.Run();
}
