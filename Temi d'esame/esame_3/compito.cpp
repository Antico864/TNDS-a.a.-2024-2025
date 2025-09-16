#include "myLibrary.hpp"

using namespace std;

int main() {
    TApplication app ("app", 0, 0);
    // Svolgimento punto 1)
    // Costanti:
    const double x0 = 1.;
    const double y0 = 0.;
    const double vx0 = 0.;
    const double vy0 = 1.;
    const double R_0 = 1.;
    const double alpha_0 = 0.;
    // Facciamo 5000 divisioni di 10*2pi:
    const double h = 2*M_PI/500.;
    vector<double> state {x0, y0, vx0, vy0};

    // Calcolo del periodo di rivoluzione:

    Runge_Kutta myrk;
    diffeq3 mydiffeq (alpha_0);
    cout << "m_B = " << pow(R_0, - 0.5 * alpha_0) << endl;
    TGraph traject;
    double t = 0.;
    unsigned int counter = 0;
    cout << "y = " << state[1] << endl;
    cout << "h = " << h << endl;
    while(state[1] >= 0.) {

        state = myrk.Step(t, state, h, mydiffeq);
   
        counter++;
        t += h;
    }
    // Adesso interpola e calcola T:
    const double y_curr = state[1];
    const double t_curr = t;
    state = myrk.Step(t, state, -h, mydiffeq);
    const double y_prev = state[1];
    double m = (y_curr - y_prev) / h;
    double t_true = t_curr - y_curr / m;
    double period = 2*t_true;
    // Così si interpola!!!!!!!

    cout << "Integrazione con passo h = " << h << endl;
    cout << "Periodo = " << period << " s" << endl;


    FillTrajectoryGraph(mydiffeq, state, myrk, period, h, traject);

    double delta = 1. - fabs(state[0]);
    // cout << state[0] << endl;
    cout << "x differisce da 1 della quantità delta = " << delta << endl;

    TCanvas c ("Trajectory", "Trajectory");

    draw(traject, c);

    traject.SetTitle("Traiettoria, #alpha = 0");

    c.SetTitle("Traiettoria, #alpha = 0");
    c.SaveAs("Images/Traiettoria.pdf");








    // Svolgimento punto 2;
    const double alpha_1 = 2.;
    const double alpha_2 = -2.;
    // Imposta i parametri iniziali:
    vector<double> state1 = {1.1, 0., 0., 1.};
    diffeq3 eq_1 (alpha_1);

    TGraph traj1;

    FillTrajectoryGraph(eq_1, state1, myrk, M_PI, h, traj1);

    // Disegnalo nel suo canvas;
    TCanvas c1 ("Trajectory 1", "Trajectory 1");
    draw(traj1, c1);
    traj1.SetTitle("Traiettoria, #alpha = 2");

    c1.SetTitle("Traiettoria1");
    c1.SaveAs("Images/Traiettoria1.pdf");
    // La traiettoria diverge!

    diffeq3 eq_2 (alpha_2);
    state1 = {1.1, 0., 0., 1.};
    TGraph traj2;
    FillTrajectoryGraph(eq_2, state1, myrk, M_PI, h, traj2);

    // Disegnalo nel suo canvas;
    TCanvas c2 ("Trajectory 2", "Trajectory 2");
    draw(traj2, c2);
    traj2.SetTitle("Traiettoria, #alpha = -2");

    c2.SetTitle("Traiettoria2");
    c2.SaveAs("Images/Traiettoria2.pdf");
    // Il moto è oscillante attorno all'origine!







    // Svolgimento punto 3)
    diffeq3_3 neweq(alpha_1);
    state1 = {1.1, 0., 0., 1.};
    TGraph traj3;
    FillTrajectoryGraph(neweq, state1, myrk, M_PI, h, traj3);

    // Disegnalo nel suo canvas;
    TCanvas c3 ("Trajectory 3", "Trajectory 3");
    draw(traj3, c3);
    traj3.SetTitle("Traiettoria con nuovo B, #alpha = -2");

    c3.SetTitle("Traiettoria3");
    c3.SaveAs("Images/Traiettoria3.pdf");
    // Bisogna disegnare un andamento di r con il tempo!! 
    // Nulla di più facile...

    state1 = {1.1, 0., 0., 1.};
    TGraph radius;
    FillRadiusGraph(neweq, state1, myrk, M_PI, h, radius);

    // Disegnalo nel suo canvas;
    TCanvas c4 ("Raggio", "Raggio");
    draw(radius, c4);
    radius.SetTitle("Raggio, #alpha = -2");
    radius.GetXaxis()->SetTitle("t [s]");
    radius.GetYaxis()->SetTitle("R [m]");

    c4.SetTitle("Raggio");
    c4.SaveAs("Images/Raggio.pdf");


    app.Run();
}