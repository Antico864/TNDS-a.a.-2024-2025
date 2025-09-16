#pragma once

#include <cmath>
#include <iostream>

#include "Functions.hpp"
#include "RandomGen.hpp"

using namespace std;

//  Classe astratta:
class IntegralMC {

public:
  IntegralMC(unsigned int seed) : m_gen(seed) {
    m_error = 0.;
  };
  virtual double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                           double sup, unsigned int points, double fmax) = 0;
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
  virtual double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                           double sup, unsigned int points, double fmax);
};

class HoMIntegrator : public IntegralMC, RandomGen {

public:
  HoMIntegrator(unsigned int seed) : IntegralMC(seed), RandomGen(seed) { ; };
  virtual double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                           double sup, unsigned int points, double fmax);
};
