#include <cstdlib>
#include <iostream>
#include <vector>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

#include "funzioni.h"

using namespace std;

int main(int argc, char **argv) {

  TApplication app("app", 0, 0);

  TGraphErrors trend;

  int index{0};

  for (int i = 1941; i < 2024; i++) {

    string filename = to_string(i) + ".txt";

    vector<double> a = Read<double>(filename.c_str());

    double avg = DoAvg(a);
    double error = DoErr(a, avg);

    cout << "Anno " << to_string(i) << ":   Delta medio = " << avg << " +/- "
         << error << endl;

    trend.SetPoint(index, i, avg);
    trend.SetPointError(index, 0, error);

    index++;
  }

  TCanvas c("Temperature trend", "Temperature trend");

  c.cd();
  c.SetGridx();
  c.SetGridy();

  trend.SetMarkerSize(0.5);
  trend.SetMarkerStyle(20);
  trend.SetFillColor(0);
  trend.SetLineColor(kBlack);

  trend.SetTitle("Temperature trend");
  trend.GetXaxis()->SetTitle("Year");
  trend.GetYaxis()->SetTitle("#Delta T [#circ C]");
  trend.Draw("APE");

  c.Update();

  c.SaveAs("Trend.pdf");

  app.Run();
}