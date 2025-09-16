#include "VectFunction.h"
#include "sign.h"

using namespace std;

vector<double> ArmonicOscillator::Eval(double t,
                                       const vector<double> &x) const {

  vector<double> v{x[1], -m_omega0 * m_omega0 * x[0]};
  return v;
}

vector<double> Euler::Step(double t, const vector<double> x, double h,
                           const BasicVectorialFunction &f) const {
  if (sign(h) != 1) {
    cerr << "Integration step must be positive. " << endl;
    exit(-1);
  };
  vector<double> solut = x + h * f.Eval(t, x);
  return solut;
}

string convert(double h) {
  int cifre_significative = -log10(h);
  ostringstream streamObj3;
  streamObj3 << fixed;
  streamObj3 << setprecision(cifre_significative);
  streamObj3 << h;
  string strObj3 = streamObj3.str();
  return strObj3;
}
