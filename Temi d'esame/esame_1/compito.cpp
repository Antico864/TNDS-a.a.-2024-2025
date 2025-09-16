#include <cstdlib>

#include "Integral.hpp"
#include "RandomGen.hpp"
#include "functions.hpp"

// In realtà questo non serve, si può anche fare una tabella.
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

int main() {

  TApplication app("app", 0, 0);
  // Valore vero dell'integrale
  const double real_value = 3. * M_E * M_E / 16.;

  // Punto 1)

  Midpoint mymidp(0., sqrt(M_E));
  Midright mymidr(0., sqrt(M_E));
  firstfunction f(3., 2.);
  double integral = 0.;

  vector<double> x(11);
  vector<double> err(11);
  vector<double> x_r(11);
  vector<double> err_r(11);

  TGraph error_graph;
  TGraph error_graph_r;

  cout << "\nReal value of integral is " << real_value << endl;
  cout << "\n=======================================" << endl;
  cout << "\tMidpoint method" << endl;
  cout << "=======================================\n" << endl;
  cout << "Number of\tIntegration\tError" << endl;
  cout << "steps\t\tstep\n" << endl;

  // Integra con midpoint:
  for (unsigned int i = 0; i < 10; i++) {
    unsigned int nsteps = pow(2, i);
    integral = mymidp.Integrate(nsteps, f);

    x[i] = sqrt(M_E) / nsteps; // Definizione di passo di integrazione
    err[i] = fabs(real_value - integral);

    cout << nsteps << "\t\t" << x[i] << " \t" << err[i] << endl;

    error_graph.SetPoint(i, x[i], err[i]);
  }

  cout << "Value of integral (1024 points): I = " << integral << endl;

  // Punto 2)

  // Così funzionano. Che bastardata ziopera...
  double k2 = log(err[0] / err[7]) / log(x[0] / x[7]);
  double k1 = pow(M_E, log(err[5]) - k2 * log(x[5]));

  cout << "error = (k_1)*h^(k_2), with k_1 = " << k1 << ", k_2 = " << k2
       << endl;

  // cout << "\nError with parameters:\n" << endl;

  // // Visualizzazione della verosimiglianza:
  // for(unsigned int i = 1; i <= 10; i++) {

  //     cout << pow(2, i) << "\t\t" << x[i-1] << " \t" << k1*pow(x[i-1], k2) <<
  //     "\t" << err[i-1] << endl;
  // }

  // Punto 3)

  cout << "\nReal value of integral is " << real_value << endl;
  cout << "\n=======================================" << endl;
  cout << "\tMidright method" << endl;
  cout << "=======================================\n" << endl;
  cout << "Number of\tValue of\tError" << endl;
  cout << "steps\t\tintegral\n" << endl;

  // Integra con midright:
  for (unsigned int i = 0; i < 10; i++) {
    unsigned int nsteps = pow(2, i + 1);
    integral = mymidr.Integrate(nsteps, f);

    x_r[i] = sqrt(M_E) / (2 * nsteps);
    err_r[i] = fabs(real_value - integral);

    cout << nsteps << "\t\t" << integral << "  \t" << err_r[i] << endl;

    error_graph_r.SetPoint(i, x_r[i], err_r[i]);
  }

  // Punto 4)
  double k2_r = log(err_r[0] / err_r[7]) / log(x_r[0] / x_r[7]);
  double k1_r = pow(M_E, log(err_r[5]) - k2 * log(x_r[5]));
  cout << log(err_r[0]) - k2 * log(x_r[0]) << endl;
  // Attenzione: gli argomenti dei logaritmi devono essere positivi!

  cout << "error = (k_1)*h^(k_2), with k_1 = " << k1_r << ", k_2 = " << k2_r
       << endl;

  // // Punti 5) e 6)

  // // Uso un integratore Montecarlo con metodo della media:
  // RandomGen myran(1);
  // AvgIntegrator myavg(1);
  // vector<double> integavg (1000);
  // // Integra con la media a 16 punti, 1000 volte:
  // for(unsigned int j = 0; j < 1000; j++) {
  //     integavg[j] = myavg.Integrate(myran, f, 0., sqrt(M_E), 16);
  // }
  // double fin_integ = DoAvg(integavg);

  // // Calcolo l'errore del midpoint a 16 punti: ce l'ho già:
  // double sig_test = err[4];

  // // Calcolo il numero minimo di punti:
  // double k = Var(integavg)*1000;
  // double avg_error = sqrt(k/1000);
  // unsigned int min_points = ceil(k/sig_test*sig_test);

  // cout << "Value of integral, average-integration: I = " << fin_integ << "±"
  // << avg_error << endl; cout << min_points << " points are needed to get an
  // error <= " << sig_test << " with average-integration. " << endl;

  // // Punto 7)
  // // In base all'andamento dell'errore userei il metodo midpoint, che sembra
  // più preciso.

  // secondfunction f2;
  // Midpoint midp2 (0., 2.);
  // double integral2 = 0.;

  // vector<double> err2(10);
  // vector<double> x2(10);

  // cout << "\nIntegral of second function " << endl;
  // cout << "\n=======================================" << endl;
  // cout << "\tMidpoint method" << endl;
  // cout << "=======================================\n" << endl;
  // cout << "Number of\tValue of\tError" << endl;
  // cout << "steps\t\tintegral\n" << endl;

  // // Integra con midpoint:
  // for(unsigned int i = 1; i <= 10; i++) {
  //     unsigned int nsteps = pow(2, i);
  //     double h = 2./nsteps;
  //     integral2 = midp2.Integrate(nsteps, f2);

  //     x2[i] = sqrt(M_E)/nsteps;
  //     // err2[i] = ??? Non sappiamo quanto vale davvero l'integrale...

  //     cout << nsteps << "\t\t" << integral2 << " \t" <<
  //     // err2[i-1] <<
  //     endl;
  // }

  // Non riesco a trovare i coefficienti giusti per il midright
  // e non riesco a fare il punto 7.

  TCanvas c;

  // Disegna error_graph su c
  c.cd();
  // c.SetLogx();
  // c.SetLogy();
  error_graph.SetMarkerStyle(20);
  error_graph.SetMarkerColor(kBlack);
  error_graph.SetMarkerSize(1);
  error_graph.SetTitle("Error of integral, midpoint");
  error_graph.GetXaxis()->SetTitle("Integration step h");
  error_graph.GetYaxis()->SetTitle("Error");
  error_graph.Draw("APC");

  c.SetTitle("Integral error");
  c.SaveAs("Images/error_midp.pdf");

  TCanvas cr;

  error_graph_r.SetTitle("Error of integral, midright");
  error_graph_r.SetMarkerStyle(20);
  error_graph_r.SetMarkerColor(kBlue);
  error_graph_r.SetMarkerSize(1);
  error_graph_r.GetXaxis()->SetTitle("Integration step h");
  error_graph_r.GetYaxis()->SetTitle("Error");
  error_graph_r.Draw("APC");

  cr.SetTitle("Integral error");
  cr.SaveAs("Images/error_midr.pdf");

  app.Run();
}