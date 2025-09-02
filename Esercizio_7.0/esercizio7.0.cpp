#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "Integral.h"
#include "funzioni.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

int main(int argc, char **argv) {

  TApplication app("app", 0, 0);

  if (argc < 2) {
    cerr << "Use of program: " << argv[0] << " <number of steps>" << endl;
    exit(-1);
  };

  unsigned int nstep = atoi(argv[1]);

  xsinx myf;

  Midpoint mymidp(0., M_PI / 2);

  const double I_value = mymidp.Integrate(nstep, myf);

  cout << "Passi= " << setw(11) << nstep << "\nI= " << setw(15)<< I_value << endl;

  // Riempimento della tabella
  // e del grafico con i punti di errore:

  TGraph error_graph;
  double I = 0.;
  double err = 0.;
  const double I_true = 1.;
  unsigned int steps = 4;
  double h = M_PI/(2*steps);

  cout << "Divisioni: " << setw(13) << "Errore: "<< endl; 

  for (unsigned int i = 0; i < 20; i++) {
    I = mymidp.Integrate(steps, myf);
    err = fabs(I_true - I);
    error_graph.SetPoint(i, h, err);
    cout << steps << "\t\t" << err << endl;
    steps*=2;
    h = M_PI/(2*steps);
  };

  TCanvas mycanv("Errore", "Errore");

  mycanv.cd();
  mycanv.SetGridx();
  mycanv.SetGridy();
  mycanv.SetLogx();
  mycanv.SetLogy();

  error_graph.SetMarkerSize(0.5);
  error_graph.SetMarkerStyle(20);
  error_graph.SetFillColor(0);
  error_graph.SetLineColor(kBlack);

  error_graph.SetTitle("Andamento dell'errore");
  error_graph.GetXaxis()->SetTitle("Passo di integrazione");
  error_graph.GetYaxis()->SetTitle("Errore sull'integrale");
  error_graph.Draw("APE");

  mycanv.Update();

  app.Run();
}