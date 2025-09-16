#pragma once

#include "VectAlg.h"
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

class Pendulum : public BasicVectorialFunction {

public:
  Pendulum(double l) { m_l = l; };

  virtual vector<double> Eval(double t, const vector<double> &x) const override;

private:
  double m_l;
  double m_gravity = 9.08665;
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
//  Integratore Runge-Kutta:
//==========================

class Runge_Kutta : public BasicDifferentialEquation {

public:
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const override;
};

string convert(double h);
