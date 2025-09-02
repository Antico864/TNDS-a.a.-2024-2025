#include "Integral.h"
#include "funzioni.h"
#include "sign.h"

using namespace std;

Integral::Integral(double a, double b) {
  CheckInterval(a, b);
  m_nstep = 0;
  m_sum = 0;
  m_integral = 0;
  m_h = 0;
};

void Integral::CheckInterval(double a, double b) {
  m_a = min(a,b);
  m_b = max(a,b);
  if(sign(b-a)==1) m_sign = 1;
  else m_sign = -1;
};

double Simpson::Integrate(unsigned int nstep, const FunzioneBase &f) {

  if (nstep <= 0 || nstep % 2 != 0) {
    throw runtime_error("Error: number of steps must be positive and even. ");
  }
  m_nstep = nstep;
  m_h = (m_b - m_a)/m_nstep;

  double sum_odd = 0., sum_even = 0.;
  for (unsigned int i = 1; i < m_nstep; i += 2) {
    sum_odd += f.Eval(m_a + i*m_h);
  }
  for (unsigned int i = 2; i < m_nstep; i += 2) {
    sum_even += f.Eval(m_a + i * m_h);
  }
  m_integral = m_sign*(m_h/3)*(f.Eval(m_a) + 4*sum_odd + 2*sum_even + f.Eval(m_b));
  return m_integral;
};
