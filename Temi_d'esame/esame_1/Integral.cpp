#include "Integral.hpp"

using namespace std;

Integral::Integral(double a, double b) {
  CheckInterval(a, b);
  m_nstep = 0;
  m_sum = 0;
  m_integral = 0;
  m_h = 0;
};

void Integral::CheckInterval(double a, double b) {
  m_a = min(a, b);
  m_b = max(a, b);
  if (a > b)
    m_sign = -1;
  else
    m_sign = 1;
};

double Midpoint::Integrate(unsigned int nstep, const BasicFunction &f) {
  if (nstep <= 0) {
    cerr << "Error: number of steps must be positive" << endl;
    return 1;
  }
  m_nstep = nstep;
  m_h = (m_b - m_a) / m_nstep;
  m_sum = 0;

  for (unsigned int i = 0; i < m_nstep; i++) {
    m_sum += f.Eval(m_a + (i + 0.5) * m_h);
  };
  m_integral = m_sign * m_sum * m_h;
  return m_integral;
};

double Midright::Integrate(unsigned int nstep, const BasicFunction &f) {
  if (nstep <= 0) {
    cerr << "Error: number of steps must be positive" << endl;
    return 1;
  }
  m_nstep = nstep;
  m_h = (m_b - m_a) / m_nstep;
  m_sum = 0;

  for (unsigned int i = 0; i < m_nstep; i++) {
    m_sum += f.Eval(m_a + i * m_h);
  };
  m_integral = m_sign * m_sum * m_h;
  return m_integral;
};

double AvgIntegrator::Integrate(RandomGen &ran, const BasicFunction &f,
                                double inf, double sup, unsigned int points) {
  double avg = 0.;
  double delta = 0.;
  for (unsigned int i = 0; i < points; i++) {
    double x = ran.Unif(inf, sup);
    delta = f.Eval(x) - avg;
    avg += delta / (i + 1);
  }
  return (sup - inf) * avg;
}