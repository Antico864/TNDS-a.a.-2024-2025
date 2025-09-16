#include "myLibrary.hpp"

using namespace std;

int main() {

    TApplication app ("app", 0, 0);
    // Svolgimento punto 1)
    // Costanti: 
    const double omega0 = 1.15;
    const double alpha = 0.01;
    const double x0 = 1.;
    const double v0 = 0.;
    double h = 0.1;
    const double sigmav0 = 0.003;

    // Vettore di stato iniziale:
    vector<double> state {x0, v0}; 
    Runge_Kutta myrk;
    DOscillator myosc (omega0, alpha);
    // Evolvi il sistema:
    double t = 0.;
    for(t = 0.; t <=43; t += h) {
        state = myrk.Step(t, state, h, myosc);
    }
    // Determina la posizione:
    const double fin_x = state[0];
    cout << "x(" << t << " s) = " << fin_x << " m" << endl;



    // Svolgimento punto 2)
    state = {x0, v0};
    double error = myrk.rkError(43, state, h, myosc);
    cout << "Passo di integrazione: h = " << h << endl;
    cout << "Errore = " << error << endl;



    // Svolgimento punto 3)
    const double epsilon = 50E-06;
    // Uso l'errore già calcolato:
    double h_nec = myrk.maxh(epsilon);
    cout << "Per ottenere un errore di 50E-06 m serve un passo di integrazione h <= " << h_nec << " s" << endl;




    // Svolgimento punto 4)
    TH1F myhist ("Distribuzione posizioni", "Posizione, t = 43 s", 100, fin_x - 5*0.001, v0 + 5*0.001);
    DoPosDistr(myhist, h_nec, sigmav0, myrk, myosc);

    TCanvas c ("Posizione", "Distribuzione posizione");
    c.cd();

    myhist.SetLineColor(kBlue);
    myhist.SetFillStyle(0);

    myhist.Draw();

    c.SaveAs("Images/hist_3");




    // Svolgimento punto 5)
    // Quello con sigma = 0.003 m/s è già fatto
    vector<double> devstd; 
    vector<double> sigmas {0.005, 0.008, 0.012, 0.015};
    cout << "Sigma v [m/s]\tsigma x [m]" << endl;
    cout << sigmav0 << "\t" << myhist.GetStdDev() << endl;
    for(unsigned int i = 0; i < 4; i ++) {
        double sigmav = sigmas[i];
        devstd.push_back(DevStd(h_nec, sigmav, myrk, myosc));
        cout << sigmav << "\t" << devstd[i] << endl;
    }

    app.Run();
}