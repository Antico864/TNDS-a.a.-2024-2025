#include <cstdlib>
#include <vector>

#include "Functions.hpp"
#include "Integral.hpp"
#include "RandomGen.hpp"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"

using namespace std;

int main() {

  TApplication app("app", 0, 0);

  unsigned int N = 500;
  unsigned int i = 2;
  vector<double> var(6);

  //  Canvas con gli istogrammi:
  TCanvas c("Convergence", "Convergenza di un integrale");
  c.Divide(3, 2); // Divido il canvas in 6 parti;

  //  Canvas della varianza:
  TCanvas c1("Variance", "Variance trend");

  TGraph vardist;
  vector<TH1F *> histovec;
  double minx, maxx;

  for (unsigned int counter = 0; counter < 6; counter++) {
    const string output = "output" + to_string(N) + "avg.txt";
    vector<double> data = Read<double>(output);
    double min = *min_element(data.begin(), data.end());
    double max = *max_element(data.begin(), data.end());

    if (counter == 0) {
      minx = min;
      maxx = max;
    }

    //  Calcolo la varianza e la metto nel suo vettore:
    var[counter] = Var(data);

    string hisname = to_string(N) + " points";
    string histitle = "N = " + to_string(N);

    TH1F *his = new TH1F(hisname.c_str(), histitle.c_str(), 100, minx, maxx);

    //  Riempio l'istogramma:
    for (unsigned int j = 0; j < data.size(); j++) {
      his->Fill(data[j]);
    }
    histovec.push_back(his);

    N *= i;
    if (i % 2 != 0) {
      i = 2;
    } else
      i = 5;
  }

  // Calcolo della costante di proporzionalità
  // tra varianza e N:
  double k = var[0] * 500;
  unsigned int N_min = ceil(k / (0.001 * 0.001));
  cout << "\nMinimum number of points for an error <= 0.001: N = " << N_min
       << "\n"
       << endl;

  for (unsigned int k = 0; k < 6; k++) {
    c.cd(k + 1);
    histovec[k]->Draw("HIST");

    histovec[k]->SetLineColor(kBlue);
    histovec[k]->SetFillStyle(0);

    c.Update();
  }

  c1.cd();

  N = 500;
  i = 2;

  for (unsigned int k = 0; k < 6; k++) {
    vardist.SetPoint(k, N, var[k]);

    N *= i;
    if (i % 2 != 0) {
      i = 2;
    } else
      i = 5;
  }

  c1.SetGridx();
  c1.SetGridy();

  c1.SetLogx();
  c1.SetLogy();

  c1.Update();

  gPad->SetLeftMargin(0.15);

  vardist.SetMarkerStyle(20);
  vardist.SetMarkerColor(kBlack);
  vardist.SetTitle("#sigma^{2} in funzione di N (metodo della media)");

  vardist.GetXaxis()->SetTitle("Numero di punti N");
  vardist.GetYaxis()->SetTitle("Varianza #sigma^{2}");

  vardist.Draw("APL");

  c.SaveAs("Images/convergenceavg.pdf");
  c1.SaveAs("Images/varianceavg.pdf");

  app.Run();
  for (unsigned int i = 0; i < histovec.size(); i++) {
    delete histovec[i];
  }
}
