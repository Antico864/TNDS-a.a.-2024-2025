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
  cout << "//=====================================================" << endl;
  cout << "//  The program produces a graph of the solution" << endl;
  cout << "//  for an armonic oscillator differential equation" << endl;
  cout << "//  between times [0, 70] s using the Runge-Kutta method.  " << endl;
  cout << "//=====================================================" << endl;
  cout << " " << endl;

  TApplication app("app", 0, 0);

  double h = atof(argv[1]);
  double tmax = 70.;
  int ind = int(tmax / h + 0.5);
  vector<double> x{0., 1.};
  double t = 0.;

  Runge_Kutta runge_kutta;

  ArmonicOscillator myosc(1.);

  TGraph SolGraph;
  TGraph error_graph;

  for (unsigned int steps = 0; steps < ind; steps++) {
    SolGraph.SetPoint(steps, t, x[0]);
    x = runge_kutta.Step(t, x, h, myosc);
    t += h;
  }

  //  Riempimento del grafico degli errori:
  for (int i = 0; i < 10; i++) {
    double h = 0.1 * pow(0.5, i); //  Ogni volta dimezza il passo;
    vector<double> x{0., 1.};
    t = 0.;
    int step_max = int(tmax / h + 0.5);

    for (int step = 0; step < step_max; step++) {
      x = runge_kutta.Step(t, x, h, myosc);
      t = t + h;
    }

    double err = fabs(x[0] - sin(t));
    error_graph.SetPoint(i, h, err);
  }

  TCanvas c("Armonic oscillator solution", "Armonic oscillator solution");
  TCanvas c_1("Runge_kutta error", "Runge_Kutta error");

  c.cd();
  c.SetGridx();
  c.SetGridy();

  SolGraph.SetMarkerSize(0.5);
  SolGraph.SetMarkerStyle(20);
  SolGraph.SetFillColor(0);
  SolGraph.SetLineColor(kBlack);

  string str =
      "Armonic oscillator solution (Runge-Kutta method, h = " + convert(h) +
      " s)";

  SolGraph.SetTitle(str.c_str());
  SolGraph.GetXaxis()->SetTitle("t [s]");
  SolGraph.GetYaxis()->SetTitle("x [m]");

  SolGraph.Draw("AL");

  c.Update();

  c_1.cd();
  c_1.SetGridx();
  c_1.SetGridy();
  c_1.SetLogx();
  c_1.SetLogy();

  error_graph.SetMarkerSize(0.5);
  error_graph.SetMarkerStyle(20);
  error_graph.SetFillColor(0);
  error_graph.SetLineColor(kBlack);

  string strg = "Runge-Kutta error, h in [0.001, 0.1] s";

  error_graph.SetTitle(strg.c_str());
  error_graph.GetXaxis()->SetTitle(" h [s] ");
  error_graph.GetYaxis()->SetTitle(" Error [m] ");

  error_graph.Draw("APL");

  c_1.Update();

  app.Run();
}