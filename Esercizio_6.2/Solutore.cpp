#include <cmath>

#include "Solutore.h"
#include "sign.h"

using namespace std;

double Bisezione::CercaZeriRef(double xmin, double xmax, const FunzioneBase &f,
                               double prec, unsigned int nmax) {
  m_niterations = 0;
  m_prec = prec;
  m_nmax = nmax;

  if (xmin < xmax) {
    m_a = xmin;
    m_b = xmax;
  } else {
    m_a = xmax;
    m_b = xmin;
  };

  double f_a = f.Eval(m_a);
  double f_b = f.Eval(m_b);

  while (fabs(m_a - m_b) > m_prec) {
    double c = 0.5 * (m_b + m_a);
    double f_c = f.Eval(c);
    if (m_niterations > m_nmax)
      break;
    m_niterations++;
    if (sign(f_a) * sign(f_c) <= 0) {
      m_b = c;
      f_b = f_c;
    } else if (sign(f_b) * sign(f_c) <= 0) {
      m_a = c;
      f_a = f_c;
    } else
      return 0.;
  };
  return 0.5 * (m_b + m_a);
};
