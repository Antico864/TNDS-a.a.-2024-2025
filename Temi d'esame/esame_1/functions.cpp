#include "functions.hpp"

using namespace std;

double firstfunction::Eval(double x) const {
  return pow(x, m_exp1) * log(sqrt(M_E + pow(x, m_exp2)));
};

double secondfunction::Eval(double x) const { return 1 / (sqrt(4 - x * x)); };

double Var(vector<double> v) {
  double sum = 0.;
  double sumq = 0.;

  if (v.size() <= 1) {
    __throw_runtime_error(
        "Errore: vettore troppo piccolo per calcolarne la varianza.");
  }

  for (unsigned int i = 0; i < v.size(); i++) {
    sum += v[i];
    sumq += v[i] * v[i];
  }

  double meanq = sumq / v.size();
  double sumofq = (sum / v.size()) * (sum / v.size());

  return meanq - sumofq;
}