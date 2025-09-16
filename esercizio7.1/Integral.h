#pragma once

#include "funzioni.h"
#include "sign.h"

using namespace std;

class Integral {

public:
  Integral(double a, double b);

  virtual double Integrate(unsigned int nstep, const FunzioneBase &f) = 0;

protected:
  void CheckInterval(double a, double b);

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

class Simpson : public Integral {

public:
  Simpson(double a, double b) : Integral(a, b) {;};

  double Integrate(unsigned int nstep, const FunzioneBase &f) override;
};
