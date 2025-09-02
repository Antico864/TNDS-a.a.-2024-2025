#pragma once

#include "VectAlg.h"
#include "sign.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

// ===========================================================================
// Classe astratta BasicVectorialFunction, restituisce il vettore df/dx(x):
// ===========================================================================
class BasicVectorialFunction {

public:
  virtual vector<double> Eval(double t, const vector<double> &x) const = 0;
};

class ArmonicOscillator : public BasicVectorialFunction {

public:
  ArmonicOscillator(double omega0) { m_omega0 = omega0; };

  //  Questo Eval non dipende da t.
  //  Questa dipendenza esiste per oscillazioni forzate, ad es...
  //  Eulero è lineare in h...
  virtual vector<double> Eval(double t, const vector<double> &x) const override;

private:
  double m_omega0;
};

class Pendulum : public BasicVectorialFunction {

public:
  Pendulum(double l) { m_l = l; };

  virtual vector<double> Eval(double t, const vector<double> &x) const override;

private:
  double m_l;
  double m_gravity = 9.08665;
};

class FDOscillator : public BasicVectorialFunction {

public:
  FDOscillator() {
    m_omega0 = 10;
    m_alpha = 1 / 30;
    m_omegaF = 5;
  };
  FDOscillator(double omega0, double alpha, double omegaF) {
    if (omega0 < alpha) {
      throw runtime_error(
          "Error: gamma must not be greater than omega. Exiting. ");
    };
    m_omega0 = omega0;
    m_alpha = alpha;
    m_omegaF = omegaF;
  };
  vector<double> Eval(double t, const vector<double> &x) const override;

private:
  double m_alpha;
  double m_omega0;
  double m_omegaF;
};

//==========================
//  Classe astratta:
//==========================

class BasicDifferentialEquation {

public:
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const = 0;
};

//==========================
//  Integratore di Eulero:
//==========================

class Euler : public BasicDifferentialEquation {

public:
  //  Valuta l'estrapolazione dopo h.
  //  Funziona per vettori di stato di taglia qualsiasi.
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const override;
};

//==============================
//  Integratore di Runge-Kutta:
//==============================

class Runge_Kutta : public BasicDifferentialEquation {

public:
  //  Valuta l'estrapolazione dopo h, sviluppato all'ordine 4.
  //  Funziona per vettori di stato di taglia qualsiasi.
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const override;
};

string convert(double h);
