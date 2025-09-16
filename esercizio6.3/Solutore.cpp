#include <cmath>

#include "Solutore.h"
#include "sign.h"

using namespace std;

double Bisezione::CercaZeriRef(double xmin, double xmax, const FunzioneBase &f,
                               double prec, unsigned int nmax) {
  m_niterations = 0;
  m_prec = prec;
  m_nmax = nmax;

  if (m_prec <= 0 || m_nmax <= 0) {
    cerr << "Errore: precisione e numero massimo di iterazioni devono essere "
            "positivi."
         << endl;
    exit(1);
  }

  m_a = xmin;
  m_b = xmax;

  double f_a = f.Eval(m_a);
  double f_b = f.Eval(m_b);

  if (f_a == 0)
    return m_a;
  if (f_b == 0)
    return m_b;

  if (sign(f_a) * sign(f_b) > 0) {
    cerr << "Errore: f(xmin) e f(xmax) devono avere segni opposti." << endl;
    exit(1);
  }

  double c = 0.0;
  while (fabs(m_a - m_b) > m_prec) {
    c = 0.5 * (m_a + m_b);
    double f_c = f.Eval(c);

    if (m_niterations >= m_nmax)
      break;
    m_niterations++;

    if (f_c == 0)
      return c;

    if (sign(f_a) * sign(f_c) < 0) {
      m_b = c;
      f_b = f_c;
    } else {
      m_a = c;
      f_a = f_c;
    }
  }

  return 0.5 * (m_a + m_b);
}