#include "myLibrary.hpp"

using namespace std;

int main() {
    // Svolgimento punto 1)

    cout << "Punto 1):\n" << endl;
    Particle Ne20 = {unitq, 20 * unitm};
    const double B0 = 0.5;
    const double v0x = 9.15E04;
    const double h = 1E-08;
    vector<double> state {0., 0., v0x, 0.};
    double d = Diameter(0., B0, Ne20, state, h).diam;
    cout << "d = " << d << " m\n\n" << endl;




    // Svolgimento punto 2)
    // Sempre con Ne20
    cout << "Punto 2):\n" << endl;
    state = {0., 0., v0x, 0.};
    double d1 = Diameter(1., B0, Ne20, state, h).diam;
    cout << "d = " << d1 << " m\n\n" << endl;



    // Svolgimento punto 3)
    cout << "Punto 3):\n" << endl;
    state = {0., 0., v0x, 0.};
    Particle Ne22 = {unitq, 22. * unitm};
    double d2 = Diameter(1., B0, Ne22, state, h).diam;
    cout << "d = " << d2 << " m\n" << endl;
    cout << "Il diametro differisce da quello del Ne20 di  delta_d = " << fabs(d1 - d2) << " m\n\n" << endl;



    // Svolgimento punto 4)
    cout << "Punto 4):\n" << endl;
    TApplication app("app", 0, 0);
    const double sigmar = 0.01;
    TH1F dhist20 ("d1 hist", "d1 hist", 100, d1 - 5. * sigmar * d1, d1 + 5. * sigmar * d1);
    TH1F dhist22 ("d2 hist", "d2 hist", 100, d2 - 5. * sigmar * d2, d2 + 5. * sigmar * d2);
    TCanvas c1("d1", "d1");
    TCanvas c2("d2", "d2");
    Results resul = SimulateTraj(B0, dhist20, dhist22, v0x, sigmar, h, Ne20, Ne22);
    // Disegna nel canvas;
    string s1 = "d1";
    draw(dhist20, c1, s1);
    string s2 = "d2";
    draw(dhist22, c2, s2);
    c1.SaveAs("Images/hist1.pdf");
    c2.SaveAs("Images/hist2.pdf");
    cout << "d1 = " << resul.davg1 << " ± " << resul.sig1 << endl;
    cout << "d2 = " << resul.davg2 << " ± " << resul.sig2 << endl;



    // Svolgimento punto 5)
    cout << "Punto 5):\n" << endl;
    if(fabs(resul.davg2 - resul.davg1) > 3 * resul.sig2) cout << "Le righe sono risolvibili;" << endl;
    else cout << "Le righe non sono risolvibili;" << endl;



    // Svolgimento punto 6)
    cout << "Punto 6):\n" << endl;
    // Fare errore runge kutta e verificare che è <= di 1/10 della delta avg
    Runge_Kutta rk;
    state = {0., 0., v0x, 0.};
    double t = Diameter(1., B0, Ne22, state, h).t;
    // Devo trovare il tempo che ci mette!!!
    state = {0., 0., v0x, 0.};
    DiffEqZ myeq (1., B0, Ne22);
    double error = rk.errorcomponent(t, state, h, myeq);
    cout << "error =" << error << endl;
    cout << "0.1 * fabs(resul.davg2 - resul.davg1) =" << 0.1 * fabs(resul.davg2 - resul.davg1) << endl;
    if(error <= 0.1 * fabs(resul.davg2 - resul.davg1)) cout << "Successo!" << endl;
    else cout << "Insuccesso!" << endl;

    app.Run();
    return 0;
}
