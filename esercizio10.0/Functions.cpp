#include "Functions.hpp"

using namespace std;

Gauss::Gauss() {
  m_mu = 0.;
  m_sigma = 1.;
}

Gauss::Gauss(double mu, double sigma) {
  m_mu = mu;
  m_sigma = sigma;
}

double Gauss::Eval(double x) const {
  return (1 / sqrt(2 * M_PI * m_sigma)) *
         pow(M_E, -pow(x - m_mu, 2) / (2 * m_sigma));
}