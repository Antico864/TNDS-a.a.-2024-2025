#include "myLibrary.hpp"

using namespace std;

int main() {
    const double d = 100E-06;
    const double L = 1.;
    const double lambda_0 = 500E-09;
    double x = 0.;
    const double lambda_1 = 400E-09;
    const double lambda_2 = 450E-09;
    function4 myf4 (lambda_0, L, d, x);

    // Svolgimento punto 1)

    TApplication app ("app", 0, 0);

    // Dichiarazione integratore e grafico:
    Trapezoidi mytrap (-d/2., d/2.); // con metodo Integrate() a precisione fissata
    TGraph mygraph;

    graph4Fill(mytrap, myf4, mygraph);

    // Svolgimento punto 2)

    // Calcolare il valore più basso di |x| per cui l’ampiezza è nulla 
    // utilizzando il metodo della bisezione con una precisione di 1 μm.
    // Guardiamo il grafico e cerchiamo in un intervallo consono. 
    // Si può fare meglio ovviamente. 

    first_zero(lambda_0, myf4, mytrap);

    // Svolgimento punto 3)
    first_zero(lambda_1, myf4, mytrap);
    first_zero(lambda_2, myf4, mytrap);







    TCanvas c;

    string s = "Figura di interferenza (#lambda = 400 nm)";
    const char* s_c = s.c_str();

    c.cd();
    c.SetGridx();
    c.SetGridy();

    mygraph.SetMarkerStyle(20);
    mygraph.SetMarkerColor(kBlack);
    mygraph.SetMarkerSize(1);
    mygraph.SetTitle(s_c);
    mygraph.GetXaxis()->SetTitle("x [m]");
    mygraph.GetYaxis()->SetTitle("A [m]");
    mygraph.Draw("AC");

    c.SetTitle("Figura di interferenza");
    c.SaveAs("Images/interferenza.pdf");

    app.Run();
}