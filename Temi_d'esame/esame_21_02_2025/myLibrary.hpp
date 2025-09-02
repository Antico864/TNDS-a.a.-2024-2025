#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"

using namespace std;

// Funzioni matematiche:

inline double sign(double x) { return (x == 0. ? 0. : (x > 0 ? 1. : -1)); };

inline int getsignfig(const double number) {
  double dummy = number;
  int j = 0;
  while (dummy < 1 || dummy >= 10) {
    if (dummy < 1) {
      dummy *= 10;
      j--;
      // cout << "Meno uno; j = " << j << endl;
    } else {
      dummy /= 10;
      j++;
      // cout << "Più uno; j = " << j << endl;
    }
  }
  cout << "j = " << j << endl;
  return j;
}

// FunzioneBase:

class FunzioneBase {
public:
  virtual double Eval(double x) const = 0;
  virtual ~FunzioneBase() { ; };
};

class functiontest : public FunzioneBase {
public:
  functiontest() { ; };
  inline double Eval(double x) const override {
    return (12. * pow(M_E, -(x + 1))) / (pow(x, 3.) + 1);
  }
  inline double Evalimp(double t) const {
    return Eval(t) / ((M_E * M_E) / (M_E - 1) * pow(M_E, -t));
  }
};

// Integrali:

class Integral {

public:
  Integral(double a, double b) {
    CheckInterval(a, b);
    m_nstep = 0;
    m_sum = 0;
    m_integral = 0;
    m_h = 0;
  }

  virtual double Integrate(unsigned int nstep, const FunzioneBase &) = 0;

protected:
  void CheckInterval(double a, double b) {
    m_a = min(a, b);
    m_b = max(a, b);
    if (a > b)
      m_sign = -1;
    else
      m_sign = 1;
  };

  unsigned int m_nstep;
  double m_a, m_b;
  double m_sum, m_integral, m_h;
  int m_sign;
};

class Simpson : public Integral {

public:
  Simpson(double a, double b) : Integral(a, b) { ; };

  inline double Integrate(unsigned int nstep, const FunzioneBase &f) override {

    if (nstep <= 0 || nstep % 2 != 0) {
      throw runtime_error("Error: number of steps must be positive and even. ");
    }

    m_nstep = nstep;
    m_h = (m_b - m_a) / m_nstep;

    double sum_odd = 0., sum_even = 0.;
    for (unsigned int i = 1; i < m_nstep; i += 2) {
      sum_odd += f.Eval(m_a + i * m_h);
    }
    for (unsigned int i = 2; i < m_nstep; i += 2) {
      sum_even += f.Eval(m_a + i * m_h);
    }
    m_integral = m_sign * (m_h / 3) *
                 (f.Eval(m_a) + 4 * sum_odd + 2 * sum_even + f.Eval(m_b));
    return m_integral;
  };
  double Error(unsigned int nstep, const FunzioneBase &f) {
    double In = 0.;
    double I2n = 0.;
    In = Integrate(nstep, f);
    I2n = Integrate(2 * nstep, f);
    return 16 / 15 * fabs(In - I2n);
  }
  inline unsigned int Nnec(double sigma, unsigned int nstep, double errnstep) {
    double hN = fabs((m_b - m_a) / nstep);
    double h = sqrt(sqrt(sigma * (pow(hN, 4) / errnstep)));
    return ceil((m_b - m_a) / h);
  }
};

// Generatori di numeri casuali:

class RandomGen {

public:
  RandomGen(unsigned int seed) {
    m_seed = seed;
    m_a = 1664525;
    m_c = 1013904223;
    m_m = 1 << 31;
  };

  inline double Rand() {
    m_seed = (m_a * m_seed + m_c) % m_m;
    return static_cast<double>(m_seed) / m_m;
  };
  inline double Unif(double a, double b) { return a + (b - a) * Rand(); };
  inline double ImpSample(double a, double b) {
    return (M_E * M_E) / (M_E - 1) * pow(M_E, -Unif(a, b));
  };

protected:
  unsigned int m_a, m_c, m_m;
  unsigned int m_seed;
};

// Algebra vettoriale:

// Somma:
template <typename T>
inline vector<T> operator+(const vector<T> &a, const vector<T> &b) {
  if (a.size() != b.size()) {
    cerr << "Vector sizes in operator+ must be equal. " << endl;
    exit(-1);
  }
  vector<T> sum(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    sum[i] = a[i] + b[i];
  }
  return sum;
};

// Differenza:
template <typename T>
inline vector<T> operator-(const vector<T> &a, const vector<T> &b) {
  if (a.size() != b.size()) {
    cerr << "Vector sizes in operator- must be equal. " << endl;
    exit(-1);
  }
  vector<T> diff(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    diff[i] = a[i] - b[i];
  }
  return diff;
};

// Prodotto scalare canonico:
template <typename T>
inline T operator*(const vector<T> &a, const vector<T> &b) {
  if (a.size() != b.size()) {
    cerr << "Vector sizes in operator* must be equal. " << endl;
    exit(-1);
  }
  T prod = a[0] * b[0];
  for (unsigned int i = 1; i < static_cast<int>(a.size()); i++) {
    prod = prod + a[i] * b[i];
  }
  return prod;
};

// Scalare per vector:
template <typename T>
inline vector<T> operator*(const T k, const vector<T> &a) {
  vector<T> prod(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    prod[i] = k * a[i];
  }
  return prod;
};

// Vector per scalare:
template <typename T>
inline vector<T> operator*(const vector<T> &a, const T k) {
  vector<T> prod(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    prod[i] = k * a[i];
  }
  return prod;
};

// vector/scalare:
template <typename T>
inline vector<T> operator/(const vector<T> &a, const T k) {
  vector<T> div(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    div[i] = (1 / k) * a[i];
  }
  return div;
};

// Funzioni vettoriali, equazioni differenziali:
class BasicVectorialFunction {

public:
  virtual vector<double> Eval(double t, const vector<double> &x) const = 0;
};

// Equazioni differenziali:
class BasicDifferentialEquation {

public:
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const = 0;
};

class Runge_Kutta : public BasicDifferentialEquation {

public:
  //  Valuta l'estrapolazione dopo h, sviluppato all'ordine 4.
  //  Funziona per vettori di stato di taglia qualsiasi.
  virtual vector<double> Step(double t, const vector<double> x, double h,
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
  // Errore runtime:
  double errorcomponent(double tstart, double tfinal, const vector<double> x,
                        double h, const BasicVectorialFunction &f) {
    Runge_Kutta rk_test;
    vector<double> sol1 = x; // Passo h
    vector<double> sol2 = x; // Passo h/2
    unsigned int counter1 = 0;
    unsigned int counter2 = 0;
    // Per controllare che avanzi di h ogni due volte che avanza di h/2:
    unsigned int ntimesh = 0;
    double t_test;
    for (t_test = tstart; t_test <= tfinal; t_test += h / 2) {
      if (counter1 % 2 == 0) {
        sol1 = rk_test.Step(t_test, sol1, h, f);
        ntimesh++;
      }
      counter1++;
      sol2 = rk_test.Step(t_test, sol2, h / 2, f);
      counter2++;
    }
    if (counter1 != counter2) {
      sol2 = rk_test.Step(t_test, sol2, h / 2, f);
      sol1 = rk_test.Step(t_test, sol1, h, f);
      counter1++;
    }
    double delta_N = fabs(sol1[0] - sol2[0]);
    m_error = 16. / 15. * delta_N; // tanto per essere sicuri,
    // ma non serve un data membro privato;
    return 16. / 15. * delta_N;
    // Ricorda: 15/16*k*h^4 = delta_N.
  }
  inline unsigned int Nnecess(unsigned int N0, double sigN0, double sigmanec,
                              const BasicVectorialFunction &f, double t0,
                              double tf) {
    double hN = fabs(1. / N0);
    double h = sqrt(sqrt(sigmanec * (pow(hN, 4) / sigN0)));
    // Numero minimo che garantisce errore = sigmanec:
    return ceil(1. / h);
  }

private:
  double m_error;
};

// Integrali MonteCarlo:
class IntegralMC {

public:
  IntegralMC(unsigned int seed) : m_gen(seed) {
    m_error = 0.;
    m_npoints = 0;
  };
  virtual double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                           double sup, unsigned int points, double fmax) = 0;
  unsigned int GetNPoints() { return m_npoints; };
  double GetError() { return m_error; };

protected:
  RandomGen m_gen;
  double m_error;
  unsigned int m_npoints;
};

class AvgIntegrator : public IntegralMC {

public:
  AvgIntegrator(unsigned int seed) : IntegralMC(seed) { ; };
  inline double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                          double sup, unsigned int points,
                          double fmax) override {
    double avg = 0.;
    double delta = 0.;
    double M2 = 0.;
    double sum = 0.;
    for (unsigned int i = 0; i < points; i++) {
      double x = ran.Unif(inf, sup);
      delta = f.Eval(x) - avg;
      avg += delta / (i + 1);
      M2 += delta * delta;
    }
    // Stima di sigma_N con un singolo integrale:
    double sig_f = sqrt(M2 / (points - 1));
    m_error = (sig_f * (sup - inf)) / sqrt(points);
    return (sup - inf) * avg;
    // In questo modo posso stimare sigma_n invocando il
    // metodo GetError() della classe madre, facendo un solo integrale.
  }
  inline double ImportanceSample(RandomGen &ran, const functiontest &f,
                                 double a, double b, unsigned int points) {
    double avg = 0.;
    double delta = 0.;
    double M2 = 0.;
    for (unsigned int i = 0; i < points; i++) {
      double x = ran.ImpSample(a, b);
      delta = f.Evalimp(x) - avg;
      avg += delta / (i + 1);
      M2 += delta * delta;
    }
    // Stima di sigma_N con un singolo integrale:
    double sig_f = sqrt(M2 / (points - 1));
    m_error = (sig_f * (b - a)) / sqrt(points);
    return (b - a) * avg;
    // In questo modo posso stimare sigma_n invocando il
    // metodo GetError() della classe madre, facendo un solo integrale.
  }
};

// ESAME VERO!!!!!

class vectftest : public BasicVectorialFunction, functiontest {
public:
  vectftest(const functiontest &f) { m_f = f; };
  inline vector<double> Eval(double t, const vector<double> &x) const override {
    // Ritorna la derivata del vettore:
    vector<double> xeval = {m_f.Eval(t)};
    return xeval;
  }

private:
  functiontest m_f;
};

inline double Nnecavg(unsigned int points0, double sig0, double sigmanec) {
  // cout << "sig0 * sqrt(points0) = " << sig0 * sqrt(points0) << endl;
  // cout << "(sigmanec * sigmanec) = " << (sigmanec * sigmanec) << endl;
  return sig0 * sqrt(points0) / (sigmanec * sigmanec);
}
