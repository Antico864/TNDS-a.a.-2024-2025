#include <cstdlib>
#include <iomanip>

#include "Integral.h"
#include "funzioni.h"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

int main(int argc, char **argv) {

  TApplication app("app", 0, 0);

  if (argc != 2) {
    cerr << "Use of program: " << argv[0] << " <number of steps>" << endl;
    exit(-1);
  };

  unsigned int nstep = atoi(argv[1]);

  if (nstep % 2 != 0) {
    cerr << "Number of steps must be even. " << endl;
    exit(-1);
  };

  xsinx myf;

  Simpson mysimp(0., M_PI/2);

  double I = mysimp.Integrate(nstep, myf);

  cout << "Numero di passi = " << nstep << "; " << endl;
  cout << "Integrale = "<< I << "; " << endl;

  TGraph error_graph;

  // Riempimento del canvas
  // con i valori dell'errore:

  double Iv = 1;
  double I_value = 0.;
  double err = 0.;
  double nsteps = 4;
  double h = M_PI/8;

  for (unsigned int i = 0; i < 40; i += 2) {
    I_value = mysimp.Integrate(nsteps, myf);
    err = fabs(Iv - I_value);
    error_graph.SetPoint(i, h, err);
    nsteps *= 2;
    h /= 2;
  };

  TCanvas mycanv("Error", "Error");

  mycanv.cd();
  mycanv.SetGridx();
  mycanv.SetGridy();
  mycanv.SetLogx();
  mycanv.SetLogy();

  error_graph.SetMarkerSize(0.5);
  error_graph.SetMarkerStyle(20);
  error_graph.SetMarkerColor(kBlack);
  error_graph.SetLineColor(kRed);
  error_graph.SetFillColor(0);

  error_graph.GetXaxis()->SetLimits(5E-07, 0.4);
  error_graph.GetYaxis()->SetRangeUser(1E-16, 5E-04);

  error_graph.SetTitle("Errore di integrazione (metodo Simpson)");
  error_graph.GetXaxis()->SetTitle("Passo di integrazione");
  error_graph.GetYaxis()->SetTitle("Errore");
  error_graph.Draw("AP");

  mycanv.Update();

  app.Run();
}