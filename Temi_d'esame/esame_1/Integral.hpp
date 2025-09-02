#pragma once

#include "RandomGen.hpp"
#include "functions.hpp"

using namespace std;

class Integral {
public:
  Integral(double a, double b);
  virtual double Integrate(unsigned int nstep, const BasicFunction &f) = 0;

protected:
  void CheckInterval(double a, double b);

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

class Midpoint : public Integral {
public:
  Midpoint(double a, double b) : Integral(a, b) { ; };
  double Integrate(unsigned int nstep, const BasicFunction &f) override;
};

class Midright : public Integral {
public:
  Midright(double a, double b) : Integral(a, b) { ; };
  double Integrate(unsigned int nstep, const BasicFunction &f) override;
};

class IntegralMC {

public:
  IntegralMC(unsigned int seed) : m_gen(seed) {
    m_error = 0.;
    m_npoints = 0;
  };
  virtual double Integrate(RandomGen &ran, const BasicFunction &f, double inf,
                           double sup, unsigned int points) = 0;
  unsigned int GetNPoints() { return m_npoints; };
  double GetError() { return m_error; };

protected:
  RandomGen m_gen;
  double m_error;
  unsigned int m_npoints;
};

//  Metodo della media:
class AvgIntegrator : public IntegralMC {

public:
  AvgIntegrator(unsigned int seed) : IntegralMC(seed) { ; };
  virtual double Integrate(RandomGen &ran, const BasicFunction &f, double inf,
                           double sup, unsigned int points);
};