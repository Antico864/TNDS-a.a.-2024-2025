#include "Functions.hpp"

using namespace std;

Gauss::Gauss() {
  m_mu = 0.;
  m_sigma = 1.;
}

Gauss::Gauss(double mu, double sigma) {
  m_mu = mu;
  m_sigma = sigma;
}

double Gauss::Eval(double x) const {
  return (1 / sqrt(2 * M_PI * m_sigma)) *
         pow(M_E, -pow(x - m_mu, 2) / (2 * m_sigma));
}

double Avg(vector<double> v) {
  if (v.size() < 2) {
    cerr << "Errore: vettore troppo piccolo per calcolare la varianza." << endl;
    exit(-1);
  }
  double avg = 0.;
  for (unsigned int i = 0; i < v.size(); i++) {
    avg += v[i];
  }
  return avg / v.size();
}

double Var(vector<double> v) {
  if (v.size() < 2) {
    cerr << "Errore: vettore troppo piccolo per calcolare la varianza." << endl;
    exit(-1);
  }

  double sumofq = 0.;
  double sum = 0.;

  for (unsigned int i = 0; i < v.size(); i++) {
    sumofq += v[i] * v[i];
    sum += v[i];
  }

  // double avg = Avg(v);
  // double deltaq = 0.;
  // for(unsigned int i = 0; i < v.size(); i++) {
  //     deltaq += (v[i] - avg)*(v[i] - avg);
  // }

  return (v.size() * sumofq - sum * sum) / (v.size() * v.size() - v.size());
}
