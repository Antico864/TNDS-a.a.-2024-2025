#include <iomanip>

#include "Prism.hpp"
#include "RandomGen.hpp"
#include "functions.hpp"

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

int main() {

  TApplication app("app", 0, 0);

  unsigned int nmeasures = 10000;

  PrismExperiment myprism(1);

  TH1F *th0_mis =
      new TH1F("#theta_{0} measured", "#theta_{0} parameter distribution", 100,
               1.5695, 1.5725);
  th0_mis->SetCanExtend(TH1::kXaxis);
  TH1F *th1_mis =
      new TH1F("#theta_{1} measured", "#theta_{1} parameter distribution", 100,
               2.548, 2.551);
  th1_mis->SetCanExtend(TH1::kXaxis);
  TH1F *th2_mis =
      new TH1F("#theta_{2} measured", "#theta_{2} parameter distribution", 100,
               2.655, 2.658);
  th2_mis->SetCanExtend(TH1::kXaxis);
  TH1F *A_mis =
      new TH1F("A measured", "A parameter distribution", 100, 2.694, 2.706);
  A_mis->SetCanExtend(TH1::kXaxis);
  TH1F *B_mis =
      new TH1F("B measured", "B parameter distribution", 100, 59E-15, 61E-15);
  B_mis->SetCanExtend(TH1::kXaxis);
  TH1F *dm1_mis =
      new TH1F("#delta_{m, 1} measured", "#delta_{m, 1} parameter distribution",
               100, 0.977, 0.9803);
  dm1_mis->SetCanExtend(TH1::kXaxis);
  TH1F *dm2_mis =
      new TH1F("#delta_{m, 2} measured", "#delta_{m, 2} parameter distribution",
               100, 1.084, 1.088);
  dm2_mis->SetCanExtend(TH1::kXaxis);
  TH1F *n1_mis = new TH1F("n_{1} measured", "n_{1} parameter distribution", 100,
                          1.6957, 1.69758);
  n1_mis->SetCanExtend(TH1::kXaxis);
  TH1F *n2_mis = new TH1F("n_{2} measured", "n_{2} parameter distribution", 100,
                          1.7502, 1.752);
  n2_mis->SetCanExtend(TH1::kXaxis);

  // Istogrammi bidimensionali:
  TH2F *AB_cov =
      new TH2F("A, B covariance", "A, B covariance ", 100, 1., 0., 100, 1., 0.);
  AB_cov->SetCanExtend(TH2::kXaxis);
  AB_cov->SetCanExtend(TH2::kYaxis);
  TH2F *dm_cov = new TH2F("#delta_{m, 1}, #delta_{m, 2} covariance",
                          "#delta_{m, 1}, #delta_{m, 2} covariance ", 100, 1.,
                          0., 100, 1., 0.);
  dm_cov->SetCanExtend(TH2::kXaxis);
  dm_cov->SetCanExtend(TH2::kYaxis);
  TH2F *n_cov = new TH2F("n_{1}, n_{2} covariance", "n_{1}, n_{2} covariance ",
                         100, 1., 0., 100, 1., 0.);
  n_cov->SetCanExtend(TH2::kXaxis);
  n_cov->SetCanExtend(TH2::kYaxis);

  TCanvas angles("Angles", "Angles distribution");
  angles.Divide(3, 1);
  TCanvas letters("A, B coefficients",
                  "A, B coefficients distribution"); // con correlazione;
  letters.Divide(3, 1);
  TCanvas dms("#delta_{m} coefficients",
              "#delta_{m} coefficients distribution"); // con correlazione;
  dms.Divide(3, 1);
  TCanvas ns("n_{i} coefficients",
             "n_{i} coefficients distribution"); // con correlazione;
  ns.Divide(3, 1);

  // Dichiarazione variabili covarianza:
  double A_sum = 0.;
  double A_sum2 = 0.;
  double B_sum = 0.;
  double B_sum2 = 0.;
  double AB_sum = 0.;
  double sig_A = 0.;
  double sig_B = 0.;

  double dm1_sum = 0.;
  double dm1_sum2 = 0.;
  double dm2_sum = 0.;
  double dm2_sum2 = 0.;
  double dm12_sum = 0.;
  double sig_dm1 = 0.;
  double sig_dm2 = 0.;

  double n1_sum = 0.;
  double n1_sum2 = 0.;
  double n2_sum = 0.;
  double n2_sum2 = 0.;
  double n12_sum = 0.;
  double sig_n1 = 0.;
  double sig_n2 = 0.;

  // Riempimento degli istogrammi
  for (unsigned int i = 0; i < nmeasures; i++) {
    myprism.Execute();
    myprism.Analyse();

    //  Riempimento istogrammi:
    // Angoli:
    th0_mis->Fill(myprism.Getth0outp());
    th1_mis->Fill(myprism.Getth1outp());
    th2_mis->Fill(myprism.Getth2outp());

    // A e B:
    A_mis->Fill(myprism.GetAoutp());
    B_mis->Fill(myprism.GetBoutp());

    // delta_m_1 e delta_m_2:
    dm1_mis->Fill(myprism.Getdm1outp());
    dm2_mis->Fill(myprism.Getdm2outp());

    // n_i:
    n1_mis->Fill(myprism.Getn1outp());
    n2_mis->Fill(myprism.Getn2outp());

    // Covarianze:
    AB_cov->Fill(myprism.GetAoutp(), myprism.GetBoutp());
    dm_cov->Fill(myprism.Getdm1outp(), myprism.Getdm2outp());
    n_cov->Fill(myprism.Getn1outp(), myprism.Getn2outp());

    // Calcolo somme correlazione:
    // A e B:
    A_sum += myprism.GetAoutp();
    A_sum2 += pow(myprism.GetAoutp(), 2);
    B_sum += myprism.GetBoutp();
    B_sum2 += pow(myprism.GetBoutp(), 2);
    AB_sum += (myprism.GetAoutp()) * (myprism.GetBoutp());

    // delta_m_1 e delta_m_2:
    dm1_sum += myprism.Getdm1outp();
    dm1_sum2 += pow(myprism.Getdm1outp(), 2);
    dm2_sum += myprism.Getdm2outp();
    dm2_sum2 += pow(myprism.Getdm2outp(), 2);
    dm12_sum += (myprism.Getdm1outp()) * (myprism.Getdm2outp());

    // n_i:
    n1_sum += myprism.Getn1outp();
    n1_sum2 += pow(myprism.Getn1outp(), 2);
    n2_sum += myprism.Getn2outp();
    n2_sum2 += pow(myprism.Getn2outp(), 2);
    n12_sum += (myprism.Getn1outp()) * (myprism.Getn2outp());
  }

  // Recupero media e la dev.std. di ogni istogramma:
  // Angoli:
  th0_mis->StatOverflows(kTRUE);
  double th0_mean = th0_mis->GetMean();
  double th0_dev = th0_mis->GetStdDev();

  th1_mis->StatOverflows(kTRUE);
  double th1_mean = th1_mis->GetMean();
  double th1_dev = th1_mis->GetStdDev();

  th2_mis->StatOverflows(kTRUE);
  double th2_mean = th2_mis->GetMean();
  double th2_dev = th2_mis->GetStdDev();

  // A e B:
  A_mis->StatOverflows(kTRUE);
  double A_mean = A_mis->GetMean();
  double A_dev = A_mis->GetStdDev();

  B_mis->StatOverflows(kTRUE);
  double B_mean = B_mis->GetMean();
  double B_dev = B_mis->GetStdDev();

  // Angoli di deviazione:
  dm1_mis->StatOverflows(kTRUE);
  double dm1_mean = dm1_mis->GetMean();
  double dm1_dev = dm1_mis->GetStdDev();

  dm2_mis->StatOverflows(kTRUE);
  double dm2_mean = dm2_mis->GetMean();
  double dm2_dev = dm2_mis->GetStdDev();

  // Indici di rifrazione:
  n1_mis->StatOverflows(kTRUE);
  double n1_mean = n1_mis->GetMean();
  double n1_dev = n1_mis->GetStdDev();

  n2_mis->StatOverflows(kTRUE);
  double n2_mean = n2_mis->GetMean();
  double n2_dev = n2_mis->GetStdDev();

  // Calcolo degli errori:
  sig_A = sqrt(A_sum2 - pow(A_sum, 2) / nmeasures) / sqrt(nmeasures);
  sig_B = sqrt(B_sum2 - pow(B_sum, 2) / nmeasures) / sqrt(nmeasures);

  sig_dm1 = sqrt(dm1_sum2 - pow(dm1_sum, 2) / nmeasures) / sqrt(nmeasures);
  sig_dm2 = sqrt(dm2_sum2 - pow(dm2_sum, 2) / nmeasures) / sqrt(nmeasures);

  sig_n1 = sqrt(n1_sum2 - pow(n1_sum, 2) / nmeasures) / sqrt(nmeasures);
  sig_n2 = sqrt(n2_sum2 - pow(n2_sum, 2) / nmeasures) / sqrt(nmeasures);

  //  Calcolo delle covarianze:
  double cov_AB =
      ((AB_sum / nmeasures - (A_sum * B_sum) / (nmeasures * nmeasures))) /
      (sig_A * sig_B);
  double cov_dm12 =
      ((dm12_sum / nmeasures - (dm1_sum * dm2_sum) / (nmeasures * nmeasures))) /
      (sig_dm1 * sig_dm2);
  double cov_n12 =
      ((n12_sum / nmeasures - (n1_sum * n2_sum) / (nmeasures * nmeasures))) /
      (sig_n1 * sig_n2);

  //  Stampa dei risultati:
  cout << "#theta_{0} = ";
  cout << PrintWithUncertaintyScient(th0_mean, th0_dev) << " rad" << endl;

  cout << "#theta_{1} = ";
  cout << PrintWithUncertaintyScient(th1_mean, th1_dev) << " rad" << endl;

  cout << "#theta_{2} = ";
  cout << PrintWithUncertaintyScient(th2_mean, th2_dev) << " rad" << endl;

  cout << "A = ";
  cout << PrintWithUncertaintyScient(A_mean, A_dev) << endl;

  cout << "B = ";
  cout << PrintWithUncertaintyScient(B_mean, B_dev) << " m^2" << endl;

  cout << "#delta_{m, 1} = ";
  cout << PrintWithUncertaintyScient(dm1_mean, dm1_dev) << " rad" << endl;

  cout << "#delta_{m, 2} = ";
  cout << PrintWithUncertaintyScient(dm2_mean, dm2_dev) << " rad" << endl;

  cout << "n1 = ";
  cout << PrintWithUncertaintyScient(n1_mean, n1_dev) << endl;

  cout << "n2 = ";
  cout << PrintWithUncertaintyScient(n2_mean, n2_dev) << endl;

  cout << "Covariance between A and B is " << fixed << setprecision(2)
       << cov_AB * 100 << "%" << endl;
  cout << "Covariance between delta m, 1 and delta m, 2 is " << fixed
       << setprecision(2) << cov_dm12 * 100 << "%" << endl;
  cout << "Covariance between n1 and n2 is " << fixed << setprecision(2)
       << cov_n12 * 100 << "%" << endl;

  //  Riempimento canvas:
  angles.cd(1);
  th0_mis->Draw("HIST");

  angles.cd(2);
  th1_mis->Draw("HIST");

  angles.cd(3);
  th2_mis->Draw("HIST");
  angles.Update();

  letters.cd(1);
  A_mis->Draw();

  letters.cd(2);
  B_mis->Draw();

  letters.cd(3);
  AB_cov->Draw("COLZ");
  letters.Update();

  dms.cd(1);
  dm1_mis->Draw();

  dms.cd(2);
  dm2_mis->Draw();

  dms.cd(3);
  dm_cov->Draw("COLZ");
  dms.Update();

  ns.cd(1);
  n1_mis->Draw();

  ns.cd(2);
  n2_mis->Draw();

  ns.cd(3);
  n_cov->Draw("COLZ");
  ns.Update();

  angles.SaveAs("Images/angles.pdf");
  letters.SaveAs("Images/coefficients.pdf");
  dms.SaveAs("Images/min_deviation.pdf");
  ns.SaveAs("Images/rifraction_indexes.pdf");

  app.Run();

  delete th0_mis;
  delete th1_mis;
  delete th2_mis;
  delete dm1_mis;
  delete dm2_mis;
  delete n1_mis;
  delete n2_mis;

  delete &angles;
  delete &letters;
  delete &dms;
  delete &ns;
}
