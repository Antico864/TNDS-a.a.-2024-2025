#include "CampoVettoriale.h"
#include "Particella.h"
#include "Posizione.h"
#include "PuntoMateriale.h"
#include "funzioni.h"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"

using namespace std;

int main(int argc, char **argv) {
  // Ricezione:
  if (argc < 4) {
    cout << "Use of program: " << argv[0] << " <x_punto> <y_punto> <z_punto>"
         << endl;
    exit(33);
  }
  double x = atoi(argv[1]);
  double y = atoi(argv[2]);
  double z = atoi(argv[3]);
  const double delta = 1E-10;
  const Posizione p(x, y, z);

  // Costruzione del dipolo:
  Elettrone e_p;
  Protone p_p;
  const Posizione pos_pr(delta / 2, 0, 0);
  const Posizione pos_el(-delta / 2, 0, 0);
  PuntoMateriale proton(p_p, pos_pr);
  PuntoMateriale electron(e_p, pos_el);

  // Calcolo del campo elettrico:
  const CampoVettoriale E_p = proton.CampoElettrico(p);
  const CampoVettoriale E_e = electron.CampoElettrico(p);
  CampoVettoriale E_tot = E_p + E_e;

  cout << " Il campo risultante in P (" << x << ", " << y << ", " << z
       << ") ha intensità E_tot = " << E_tot.Modulo()
       << ".\n Le componenti del campo sono: " << endl;
  cout << " Ex = " << E_tot.getEx() << "\n Ey = " << E_tot.getEy()
       << "\n Ez = " << E_tot.getEz() << endl;

  // Calcolo della costante di potenza:
  double alpha = 0.;
  Posizione test(1, 1, 1);
  if (p != test) {
    CampoVettoriale E1 = proton.CampoElettrico(p) + electron.CampoElettrico(p);
    CampoVettoriale E2 =
        proton.CampoElettrico(test) + electron.CampoElettrico(test);
    alpha = log(E1.Modulo() / E2.Modulo()) / log(p.getR() / test.getR());
  } else {
    Posizione test2(2, 2, 2);
    CampoVettoriale E1 = proton.CampoElettrico(p) + electron.CampoElettrico(p);
    CampoVettoriale E2 =
        proton.CampoElettrico(test2) + electron.CampoElettrico(test2);
    alpha = log(E1.Modulo() / E2.Modulo()) / log(p.getR() / test2.getR());
  }
  cout << "E ~ k•R^(" << alpha << ") for R -> +∞ " << endl;

  // Grafico del campo lungo l'asse:
  TApplication app("app", 0, 0);

  TGraph *graph = new TGraph;

  for (int i = 100; i <= 1000; i++) {
    Posizione p(i * delta, 0, 0);
    CampoVettoriale E_sum =
        electron.CampoElettrico(p) + proton.CampoElettrico(p);

    graph->SetPoint(i - 100, i * delta, E_sum.Modulo());
  }

  TCanvas *c1 = new TCanvas("c1", "Campo elettrico di dipolo lungo l'asse");

  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  // c1->SetLogx();
  c1->SetLogy();

  graph->SetMarkerSize(0.5);
  graph->SetMarkerStyle(20);
  graph->SetFillColor(0);
  graph->SetLineColor(kBlack);

  graph->SetTitle("Campo elettrico di dipolo lungo l'asse");
  graph->GetXaxis()->SetTitle("x [m]");
  graph->GetYaxis()->SetTitle("E_{tot} [N/C]");
  graph->GetXaxis()->SetLimits(50 * delta, 1050 * delta);
  graph->Draw("APE");

  c1->Update();

  c1->SaveAs("Campo_dipolo.pdf");

  app.Run();
}
