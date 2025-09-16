#pragma once

#include "funzioni.h"
#include "sign.h"

using namespace std;

class Integral {

public:
  Integral(double a, double b) {
    CheckInterval(a, b);
    m_nstep = 0;
    m_sum = 0;
    m_integral = 0;
    m_h = 0;
  }

  virtual double Integrate(unsigned int nstep, const FunzioneBase &) = 0;

protected:
  void CheckInterval(double a, double b) {
    m_a = min(a,b);
    m_b = max(a,b);
    if(sign(b-a) == 1) m_sign = 1;
    else m_sign = -1;
  };

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

class Midpoint : public Integral {

public:
  Midpoint(double a, double b) : Integral(a, b) { ; };

  double Integrate(unsigned int nstep, const FunzioneBase &f) override;
};
