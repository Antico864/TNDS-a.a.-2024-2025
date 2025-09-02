#pragma once

#include <cmath>
#include <iostream>

#include "Functions.hpp"

using namespace std;

class RandomGen {

public:
  RandomGen(unsigned int);

  double Rand();
  double Unif(double, double);
  double Exp(double);
  double GaussBM(double, double);
  double FunctionAR(const FunzioneBase &, double, double, double);

protected:
  unsigned int m_a, m_c, m_m;
  unsigned int m_seed;
};
