#include <iomanip>
#include <iostream>

#include "Integral.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "funzioni.h"

using namespace std;

int main(int argc, char **argv) {

  TApplication app("app", 0, 0);

  if (argc < 3) {
    cerr << "Use of program: " << argv[0] << " <numero di dev.std.>"
         << " <precisione>" << endl;
    exit(-1);
  }

  double t = atof(argv[1]);
  double prec = atof(argv[2]);

  if (fabs(t) > 5.) {
    cerr << "Error: t must be between -5 and 5. " << endl;
    exit(-1);
  }

  NormGaussian myGaussian;

  Trapezoidi trapezio(-t, t);

  double prob = trapezio.GaussProb(myGaussian, prec);
  double error = (4 / 3) * fabs(trapezio.GetIn1() - trapezio.GetIn());
  int signfig = abs(getsignfig(error));

  cout << fixed;
  cout << "Gaussian probability for measure to fall in [-t*std.dvt., "
          "t*std.dvt.] = "
       << setprecision(signfig) << prob << "; " << endl;
  cout << "Error = " << setprecision(signfig) << error << "; " << endl;

  TGraph GaussProbability;
  TLegend legend(0.7, 0.7, 0.9, 0.9);
  TCanvas mycanv;

  double I = 0.;
  int i = 0;

  for (double x = 0.; x <= 5; x += 0.08) {
    Trapezoidi trap(-x, x);
    I = trap.GaussProb(myGaussian, prec);
    GaussProbability.SetPoint(i, x, I);
    if (I > 1) {
      cerr << "Error: probability must not exceed 1. " << endl;
      cerr << "Number of iteration = " << i << endl;
    }
    i++;
  }

  ////////////////////////////////////
  //  La distribuzione salta a valori prossimi allo zero vicino a 5 sigma dalla
  //  media. Il valore di x a cui avviene il salto è circa inversamente
  //  proporzionale alla precisione...
  ////////////////////////////////////

  mycanv.cd();
  mycanv.SetGridx();
  mycanv.SetGridy();

  GaussProbability.SetMarkerSize(0.5);
  GaussProbability.SetMarkerStyle(20);
  GaussProbability.SetMarkerColor(kBlack);
  GaussProbability.SetLineColor(kBlue);
  GaussProbability.SetFillColor(0);

  GaussProbability.SetTitle("Gaussian probability");
  GaussProbability.GetXaxis()->SetTitle(
      "Half interval width [units of std.dvt.]");
  GaussProbability.GetYaxis()->SetTitle("Probability");
  GaussProbability.Draw("APL");

  app.Run();

  return 0;
}