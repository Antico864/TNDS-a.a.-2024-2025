#include "Integral.h"
#include "funzioni.h"

using namespace std;

Integral::Integral(double a, double b) {
  CheckInterval(a, b);
  m_nstep = 0;
  m_sum = 0;
  m_integral = 0;
  m_h = 0;
}

void Integral::CheckInterval(double a, double b) {
  m_a = min(a, b);
  m_b = max(a, b);
  if (sign(b - a) == 1) m_sign = 1;
  else m_sign = 1;
};

double Trapezoidi::Integrate(double prec, const FunzioneBase &f) {
  m_sum = (f.Eval(m_a) + f.Eval(m_b)) / 2.0;
  m_nstep = 1;
  m_integral = (m_b - m_a) * m_sum;
  m_h = (m_b - m_a);

  double x = 0.;
  I_n = 0.;
  I_n1 = m_h * (f.Eval(m_a) + f.Eval(m_b)) / 2.0;

  while (fabs(I_n1 - I_n) >= prec) {

    I_n = I_n1; // Aggiorno I_n;
    m_sum = 0.; // Resetto la somma;

    // Calcolo la somma sulla partizione:
    for (unsigned int i = 0; i < m_nstep; i++) {
      x = m_a + (i + 0.5) * m_h;
      m_sum += f.Eval(x);
    };

    I_n1 = (I_n + m_h * m_sum) / 2; // Aggiorno  I_n1;
    m_h /= 2.; // Aggiorno h:
    m_nstep *= 2;
    // Debugging:
    // cout << "m_h = " << m_h << endl;
    // cout << "m_sum = " << m_sum << endl;
  };

  m_integral = m_sign * I_n1;
  return m_integral;
};
