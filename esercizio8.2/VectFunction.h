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

//==========================
//  Classe virtuale:
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

class Runge_Kutta : public BasicDifferentialEquation {

public:
  //  Valuta l'estrapolazione dopo h, sviluppato all'ordine 4.
  //  Funziona per vettori di stato di taglia qualsiasi.
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const override;
};

string convert(double h);
