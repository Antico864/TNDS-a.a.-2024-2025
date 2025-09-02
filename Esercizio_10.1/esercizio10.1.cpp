#include <cstdlib>
#include <vector>

#include "RandomGen.hpp"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"

using namespace std;

int main() {

  TApplication myapp("myapp", 0, 0);

  RandomGen myran(1);

  unsigned int Nentries = 200000;
  unsigned int NMax = 12;

  vector<double> entries(Nentries);

  for (unsigned int j = 0; j < Nentries; j++) {
    double rand = myran.Rand();
    entries[j] = rand;
  }

  TGraph deviation;
  vector<double> varianza(NMax);

  TCanvas c("Central limit", "Verification of central limit theorem");
  c.Divide(4, 3); // Divido il canvas in 12 parti;

  for (unsigned int N = 1; N < NMax + 1; N++) {

    string name = "histo n#circ " + to_string(N);
    TH1F *histo = new TH1F(name.c_str(), "", 200, 0., 1.);

    histo->SetMinimum(0.0);

    unsigned int sumSize = trunc((Nentries) / N);

    vector<double> sum(sumSize);

    // Riempio gli istogrammi, dividendo i dati in gruppi di N numeri:
    for (unsigned int k = 0; k < sumSize; k++) {
      for (unsigned int i = k * N; i < (k + 1) * N && i < Nentries; i++) {
        sum[k] += entries[i];
      }
      histo->Fill(sum[k] / N);
    }

    histo->SetLineColor(kBlue);
    histo->SetFillStyle(0);

    c.cd(N); // Disegna nel posto giusto:
    histo->Draw();
    c.Update(); // Aggiorna il canvas per ogni disegno:

    // Calcola la varianza:
    double var = Var(sum);
    if (var != var) { // Controlla se la varianza è 'nan'
      cerr << "Errore nel calcolo della varianza per N = " << N << endl;
      var = 0;
    }
    // Riempi correttamente il vector della varianza:
    varianza[N - 1] = var;

  }

  c.Update();
  c.SaveAs("Images/central_limit.pdf");

  TCanvas c1("Variance", "Variance distribution");

  for (unsigned int i = 0; i < 12; i++) {
    deviation.SetPoint(i, i + 1, varianza[i]);
  }

  c1.SetGrid();

  deviation.SetMarkerStyle(25);
  deviation.SetMarkerSize(1.0);
  deviation.SetMarkerColor(kBlack);
  deviation.SetTitle("#sigma^{2} vs N");
  deviation.GetXaxis()->SetTitle("N");
  deviation.GetYaxis()->SetTitle("#sigma^{2}");

  deviation.Draw("APL");

  c1.Update();
  c1.SaveAs("Images/variance_graph.pdf");

  myapp.Run();
}