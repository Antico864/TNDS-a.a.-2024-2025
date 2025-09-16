#pragma once

#include <cmath>
#include <iostream>

using namespace std;

class FunzioneBase {

public:
  virtual double Eval(double x) const = 0;
};

class Gauss : public FunzioneBase {

public:
  Gauss();
  Gauss(double, double);
  virtual double Eval(double x) const override;

private:
  double m_mu, m_sigma;
};
