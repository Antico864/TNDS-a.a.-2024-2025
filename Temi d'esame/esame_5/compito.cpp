#include "myLibrary.hpp"
#include "TROOT.h"

using namespace std;

int main() {
  // Svolgimento punto 1)
  const double err_rel = 0.03;
  ExpResults res = RunExperiment(10000, 1, err_rel, err_rel, err_rel, err_rel);
  double C_avg = res.C_final;
  double err_rel_C = res.sigmarel_C;
  cout << "C misurata: C = " << C_avg << " F" << endl;
  cout << "Errore relativo su C: err = " << fixed << setprecision(3)
       << err_rel_C * 100 << "%" << endl;





  // Svolgimento punto 2)


  // Solo errore su V0
  ExpResults res1 = RunExperiment(10000, 1, err_rel, 0., 0., 0.);
  // Solo errore su V1
  ExpResults res2 = RunExperiment(10000, 1, 0., err_rel, 0., 0.);
  // Solo errore su R
  ExpResults res3 = RunExperiment(10000, 1, 0., 0., err_rel, 0.);
  // Solo errore su t
  ExpResults res4 = RunExperiment(10000, 1, 0., 0., 0., err_rel);
  vector<ExpResults> vecres = {res1, res2, res3, res4};
  vector<char *> vars = {"V0", "V1", "R", "t"};
  double error = res1.sigmarel_C;
  unsigned int counter = 0;
  for (unsigned int i = 0; i < 4; i++) {
    if (vecres[i].sigmarel_C >= error) {
      error = vecres[i].sigmarel_C;
      counter = i;
    }
    cout << vars[i] << "  " << vecres[i].sigmarel_C * 100 << endl;
  }
  cout << "Il maggior contributo all'errore è portato dalla grandezza "
       << vars[counter] << endl;
  cout << "Errore percentuale a C con contributo solo di " << vars[counter]
       << " è err = " << vecres[counter].sigmarel_C * 100 << "%" << endl;

  // Svolgimento punto 3)

  vector<double> errorvec;
  vector<double> sigmavec;
  FillWithError(sigmavec, errorvec, err_rel);

  TApplication app ("myapp", 0, 0);
  
  TGraph newgraph;
  for(unsigned int i = 0; i < sigmavec.size(); i++) {
    newgraph.SetPoint(i, 100*sigmavec[i], 100*errorvec[i]);
  }

  TCanvas c3 ("Errore", "Errore");

  c3.cd();
  newgraph.SetMarkerSize(0.5);
  newgraph.SetMarkerStyle(20);
  newgraph.SetMarkerColor(kBlue);
  newgraph.SetFillColor(0);
  newgraph.SetLineColor(kBlack);

  string str = "#sigma_{C}% in dipendenza da #sigma_{V}% ";

  newgraph.SetTitle(str.c_str());
  newgraph.GetXaxis()->SetTitle("#sigma_{V}% [V]");
  newgraph.GetYaxis()->SetTitle("#sigma_{C}% [F]");

  newgraph.Draw("APL");
  





  // Punto 4)
  double m = (sigmavec[7] - sigmavec[6]) / (errorvec[7] - errorvec[6]);
  double sig_nec = sigmavec[6] + (0.07 - errorvec[6]) * m;
  cout << "Per ottenere un errore del 7% su C è necessaria un'incertezza "
          "relativa su V dell' "
       << sig_nec*100 << " %. " << endl;
  // È circa del 5%. 


  c3.Update();
  c3.SaveAs("Images/Errors.pdf");

  app.Run();

}
