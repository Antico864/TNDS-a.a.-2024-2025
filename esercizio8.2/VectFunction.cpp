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

vector<double> Runge_Kutta::Step(double t, const vector<double> x, double h,
                                 const BasicVectorialFunction &f) const {
  if (sign(h) != 1) {
    cerr << "Integration step must be positive. " << endl;
    exit(-1);
  };
  vector<double> k1 = f.Eval(t, x);
  vector<double> k2 = f.Eval(t + h / 2, x + h / 2 * k1);
  vector<double> k3 = f.Eval(t + h / 2, x + h / 2 * k2);
  vector<double> k4 = f.Eval(t + h, x + h * k3);
  vector<double> solution = x + (k1 + 2. * k2 + 2. * k3 + k4) * (h / 6);
  return solution;
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
