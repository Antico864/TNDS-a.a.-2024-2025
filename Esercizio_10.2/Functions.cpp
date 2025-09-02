#include "Functions.hpp"

using namespace std;

double Avg(const vector<double> &v) {
  if (v.size() == 1) {
    return v[0];
  }
  double avg = 0.;
  for (unsigned int i = 0; i < v.size(); i++) {
    double delta = v[i] - avg;
    avg += delta / (i + 1);
  }
  return avg;
}

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
