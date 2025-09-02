#include <cstdlib>
#include <iostream>
#include <vector>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "funzioni.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 3) {
    cout << "Uso del programma : " << argv[0] << " ndata" << " <filename> "
         << endl;
    return -1;
  }

  TApplication app("app", 0, 0);

  vector<double> a = Read<double>(atoi(argv[1]), argv[2]);
  for (int k = 0; k < a.size(); k++) {
    cout << a[k] << endl;
  }

  double min_val = *min_element(a.begin(), a.end()) - 1;
  double max_val = *max_element(a.begin(), a.end()) + 1;

  TH1F *histo =
      new TH1F("histo", "Andamento temperatura", 100, min_val, max_val);
  histo->StatOverflows(kTRUE);
  for (int k = 0; k < a.size(); k++)
    histo->Fill(a[k]);

  cout << "Media dei valori = " << histo->GetMean() << endl;

  TCanvas *mycanvas = new TCanvas("Histo", "Histo");
  histo->Draw();
  histo->GetXaxis()->SetTitle("#Delta T");
  histo->GetYaxis()->SetTitle("Number of counts");
  mycanvas->Update();

  mycanvas->SaveAs("Delta_T.pdf");

  app.Run();
}