#include "Integral.h"
#include "funzioni.h"

using namespace std;

double Midpoint::Integrate(unsigned int nstep, const FunzioneBase &f) {
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
