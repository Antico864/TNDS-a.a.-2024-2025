#pragma once

#include <cmath>
#include <iostream>

#include "functions.hpp"

using namespace std;

class RandomGen {

public:
  RandomGen(unsigned int);

  double Rand();
  double Unif(double a, double b);

protected:
  unsigned int m_a, m_c, m_m;
  unsigned int m_seed;
};
