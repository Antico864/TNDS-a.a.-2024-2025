#include "RandomGen.hpp"

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

int main(int argc, const char **argv) {

  if (argc < 2) {
    cerr << "Use of program: " << argv[0] << " <center of gaussian> "
         << " <sigma> " << endl;
    exit(-1);
  };

  double mu = atof(argv[1]);
  double sigma = atof(argv[2]);

  TApplication myapp("myapp", 0, 0);

  RandomGen myran(1);
  Gauss mygauss(mu, sigma);
  double max = 1. / sqrt(2 * M_PI * sigma);

  TH1F uniform("uniform", "Uniform distribution", 100, 5., 10.);
  TH1F exponential("exponential", "Exponential distribution", 100, 0., 10.);
  TH1F gaussianBM("gaussian Box-Muller", "Gaussian distribution, Box-Muller",
                  100, mu - 4*sigma, mu + 4*sigma);
  TH1F gaussianAR("gaussian accept-reject",
                  "Gaussian distribution, accept-reject", 100, mu - 4*sigma, mu + 4*sigma);

  for (unsigned int i = 0; i < 10000; i++) {
    double x = myran.Unif(5., 10.);
    uniform.Fill(x);
    double y = myran.Exp(1.);
    exponential.Fill(y);
    double z = myran.GaussBM(mu, sigma);
    gaussianBM.Fill(z);
    double t = myran.FunctionAR(mygauss, max, mu - 3 * sigma, mu + 3 * sigma);
    gaussianAR.Fill(t);
  }

  TCanvas c("Uniform", "Uniform distribution");

  c.cd();

  uniform.SetLineColor(kBlue);
  uniform.SetFillStyle(0);

  uniform.SetTitle("Uniform distribution in [5, 10]");

  uniform.Draw();

  c.Update();

  TCanvas c1("Exponential", "Exponential distribution");

  c1.cd();

  exponential.SetLineColor(kBlue);
  exponential.SetFillStyle(0);

  exponential.SetTitle("Exponential distribution");

  exponential.Draw();

  c1.Update();

  TCanvas c2("gaussian accept-reject", "Gaussian distribution, accept-reject");

  c2.cd();

  gaussianBM.SetLineColor(kBlue);
  gaussianBM.SetFillStyle(0);

  gaussianBM.SetTitle("Gaussian distribution, Box-Muller method");

  gaussianBM.Draw();

  c2.Update();

  TCanvas c3("gaussian accept-reject", "Gaussian distribution, accept-reject");

  c3.cd();

  gaussianAR.SetLineColor(kBlue);
  gaussianAR.SetFillStyle(0);

  gaussianAR.SetTitle("Gaussian distribution, accept-reject method");

  gaussianAR.Draw();

  c3.Update();

  myapp.Run();
}