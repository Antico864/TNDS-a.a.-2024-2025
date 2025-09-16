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
  vector<double> variance(6);

  //  Canvas con gli istogrammi:
  TCanvas c("Convergence", "Convergence of an integral");
  c.Divide(3, 2); // Divido il canvas in 6 parti;

  //  Canvas della varianza:
  TCanvas c1("Variance", "Variance trend");
  //  Faccio un vector di puntatori a TH1F*,
  //  che riempio nel ciclo for.

  TGraph vardist;
  vector<TH1F *> histovec;
  double minx, maxx;

  for (unsigned int counter = 0; counter < 6; counter++) {
    const string output = "output" + to_string(N) + "hom.txt";
    vector<double> data = Read<double>(output);
    double min = *min_element(data.begin(), data.end());
    double max = *max_element(data.begin(), data.end());

    if (counter == 0) {
      minx = min;
      maxx = max;
    }

    //  Calcolo la varianza e la metto nel suo vettore:
    variance[counter] = Var(data);

    string hisname = to_string(N) + " points";
    string histitle = "N = " + to_string(N);

    TH1F *his = new TH1F(hisname.c_str(), histitle.c_str(), 84, minx, maxx);

    //  Riempio l'istogramma counter-esimo
    //  e lo metto nel vettore;
    for (unsigned int j = 0; j < data.size(); j++) {
      his->Fill(data[j]);
    }
    histovec.push_back(his);

    //  Aggiorno N;
    N *= i;
    if (i % 2 != 0) {
      i = 2;
    } else
      i = 5;
  }

  // Calcolo della costante di proporzionalità
  // tra varianza e N:
  double k = variance[0] * 500;
  unsigned int N_min = ceil(k / (0.001 * 0.001));
  cout << "\nMinimum number of points for an error of 0.001: N = " << N_min
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
    vardist.SetPoint(k, N, variance[k]);

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

  //  Per far vedere il titolo dell'asse y;
  gPad->SetLeftMargin(0.15);

  vardist.SetMarkerStyle(20);
  vardist.SetMarkerColor(kBlack);
  vardist.SetTitle("#sigma^{2} in funzione di N (metodo Hit-or-Miss)");

  vardist.GetXaxis()->SetTitle("Numero di punti N");
  vardist.GetYaxis()->SetTitle("varianza #sigma^{2}");

  vardist.Draw("APL");

  c.SaveAs("Images/convergencehom.pdf");
  c1.SaveAs("Images/variancehom.pdf");

  app.Run();

  for (unsigned int i = 0; i < histovec.size(); i++) {
    delete histovec[i];
  }
}
