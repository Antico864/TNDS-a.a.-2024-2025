#pragma once

#include "funzioni.h"
#include "sign.h"
#include <cassert>

using namespace std;

class Integral {

public:
  Integral(double a, double b);

  virtual double Integrate(unsigned int nstep, const FunzioneBase &f) = 0;
  virtual double Integrate(double prec, const FunzioneBase &f) = 0;

  // Metodi Get:

  double GetH() { return m_h; };
  double GetN() { return m_nstep; };

protected:
  void CheckInterval(double a, double b);

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

class Trapezoidi : public Integral {

public:
  Trapezoidi(double a, double b) : Integral(a, b) { ; };

  ///////////////
  //  Metodi Get:
  ///////////////

  double GetIn() { return I_n; };
  double GetIn1() { return I_n1; };

  ///////////////
  //  Metodi Integrate:
  ///////////////

  double Integrate(unsigned int nstep, const FunzioneBase &f) override {
    cerr << "Method not implemented for Trapezoidi." << endl;
    return 0;
  }

  double Integrate(double prec, const FunzioneBase &f) override;

  double GaussProb(const NormGaussian &h, double prec) {
    double prob = (Trapezoidi::Integrate(prec, h));
    return prob;
  }

private:
  double I_n, I_n1;
  // Sono in più rispetto ai data membri di Integral
};

bool are_close(double a, double b, double eps = 1E-7) {
  return fabs(a - b) < eps;
};
