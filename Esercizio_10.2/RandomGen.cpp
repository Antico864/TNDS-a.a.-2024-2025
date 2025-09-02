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
