#include "RandomGen.hpp"

using namespace std;

RandomGen::RandomGen(unsigned int seed) {
  m_seed = seed;
  m_a = 1664525;
  m_c = 1013904223;
  m_m = 1 << 31;
}

double RandomGen::Rand() {
  m_seed = (m_a * m_seed + m_c) % m_m;
  return static_cast<double>(m_seed) / m_m;
}

double RandomGen::Unif(double a, double b) { return a + (b - a) * Rand(); }

double RandomGen::Exp(double lambda) { return -(1 / lambda) * log(1 - Rand()); }

double RandomGen::GaussBM(double avg, double sigma) {
  double s = Rand();
  double t = Rand();
  return avg + sigma * sqrt(-2 * log(1 - s)) * cos(2 * M_PI * t);
}

//  f deve essere sempre positiva e limitata!
double RandomGen::FunctionAR(const FunzioneBase &f, double max, double a,
                             double b) {
  double s = Rand();
  double t = Rand();
  double x = a - (b - a) * s;
  double y = max * t;
  while (y > f.Eval(x)) {
    s = Rand();
    t = Rand();
    x = a + (b - a) * s;
    y = max * t;
  }
  return x;
}
