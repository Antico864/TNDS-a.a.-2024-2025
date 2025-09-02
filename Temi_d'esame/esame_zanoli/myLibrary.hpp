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
// PRIMA DI CONSEGNARE LO SCRITTO!!!!
// clang-format -i *.cpp *.hpp

// Per impostare i limiti dinamicamente:
// TH1F histo; // con dettagli
// histo.SetCanExtend(TH1F::kAllAxes);

// histo.SetStatOverflows(kTrue); per far edev std e media con tutti i dati
// dell'istogramma, anche quelli fuori limite dell'asse.

// Conversione da double a stringa:
inline string convert(double h) {
  int cifre_significative = -log10(h);
  ostringstream streamObj3;
  streamObj3 << fixed;
  streamObj3 << setprecision(cifre_significative);
  streamObj3 << h;
  string strObj3 = streamObj3.str();
  return strObj3;
}

// Funzioni con vector:

// Ricorda: c'è un metodo che si chiama sort(v.begin(), v.end()).
// USALO!!! Non serve fare i fighi con il MergeSort!!

template <typename T> vector<T> Read(const char *filename) {
  vector<T> v;
  T val;
  ifstream inputFile(filename);
  if (!inputFile) {
    cout << "Cannot open file " << filename << endl;
    exit(11);
  } else {
    while (!inputFile.eof()) {
      inputFile >> val;
      v.push_back(val);
    };
  };
  inputFile.close();
  return v;
};

template <typename T> double DoAvg(vector<T> a) {
  double avg = 0;
  for (int k = 0; k < a.size(); k++) {
    avg = (avg * (k) + a[k + 1]) / (k + 1);
  }
  return avg;
};

template <typename T> double DevStd(vector<T> a) {
  double delta2sum = 0;
  unsigned int counter = 0;
  double avg = DoAvg(a);
  for (unsigned int k = 0; k < a.size(); k = k + 7) {
    delta2sum += pow(a[k] - avg, 2);
    counter++;
  }
  double devstdavg = sqrt(delta2sum / (counter - 1)) / sqrt(counter);
  return devstdavg;
};

template <typename T> double CalcMed(vector<T> a) {
  double mediana = 0;
  int mN = a.size();
  sort(a.begin(), a.end());
  if (mN % 2 != 0) {
    mediana = a[(mN + 1) / 2];
  } else {
    mediana = (a[mN / 2 - 1] + a[mN / 2]) / 2;
  }
  return mediana;
};

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

class Parabola : public FunzioneBase {
public:
  Parabola() {
    m_a = 0;
    m_b = 0;
    m_c = 0;
  };
  Parabola(double a, double b, double c) {
    m_a = a;
    m_b = b;
    m_c = c;
  };
  ~Parabola() { ; };
  inline double Eval(double x) const override {
    return m_a * x * x + m_b * x + m_c;
  };

  ///////////////////////
  //  Metodi Set:
  ///////////////////////

  void SetA(double a) { m_a = a; };
  void SetB(double b) { m_b = b; };
  void SetC(double c) { m_c = c; };

  ///////////////////////
  //  Metodi Get:
  ///////////////////////

  double GetA() { return m_a; };
  double GetB() { return m_b; };
  double GetC() { return m_c; };

  double GetVertex() const { return -m_b / (2 * m_a); };

private:
  double m_a, m_b, m_c;
};

class Tan_men_x : public FunzioneBase {
public:
  Tan_men_x() {
    m_A = 1.;
    m_omega = 1.;
    m_phi = 0.;
  }

  inline double Eval(double x) const override { return sin(x) - x * cos(x); }

private:
  double m_A, m_omega, m_phi;
};

class xsinx : public FunzioneBase {

public:
  xsinx() {
    m_A = 1.;
    m_omega = 1.;
    m_phi = 0.;
  };
  xsinx(double A, double omega, double phi) {
    m_A = A;
    m_omega = omega;
    m_phi = phi;
  };

  inline double Eval(double x) const override { return x * sin(x); };

private:
  double m_A, m_omega, m_phi;
};

// Solutore e bisezione per zeri:

class Solutore {

public:
  ///////////////////
  //  Costruttori
  ///////////////////

  Solutore() { ; };
  Solutore(double prec) { m_prec = prec; };
  Solutore(unsigned int n_max) { m_nmax = n_max; };

  virtual ~Solutore() { ; };

  ///////////////////
  //  CercaZeri:
  ///////////////////

  virtual double CercaZeriRef(double xmin, double xmax, const FunzioneBase &f,
                              double prec, unsigned int nmax) = 0;

  ///////////////////
  //  Metodi Get:
  ///////////////////

  unsigned int GetNMaxIt() { return m_nmax; };
  unsigned int GetNIt() { return m_niterations; };
  double GetPrec() { return m_prec; };

  ///////////////////
  //  Metodi Set:
  ///////////////////

  void SetNMaxIt(unsigned int nint) { m_nmax = nint; };
  void SetPrec(double epsilon) { m_prec = epsilon; };

protected:
  double m_a, m_b; //  Estremi intervallo di ricerca
  double m_prec;
  unsigned int m_nmax;        //  Massime
  unsigned int m_niterations; //  Effettuate
};

class Bisezione : public Solutore {

public:
  ///////////////////
  //  Costruttori:
  ///////////////////

  Bisezione() { ; };
  Bisezione(double prec) { m_prec = prec; };
  Bisezione(unsigned int n) { m_nmax = n; };
  virtual ~Bisezione() { ; };

  ///////////////////
  //  CercaZeri:
  ///////////////////

  inline double CercaZeriRef(double xmin, double xmax, const FunzioneBase &f,
                             double prec, unsigned int nmax) override {
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
};

// Integrali, con midpoint e Simpson:

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

class Midpoint : public Integral {

public:
  Midpoint(double a, double b) : Integral(a, b) { ; };

  inline virtual double Integrate(unsigned int nstep,
                                  const FunzioneBase &f) override {
    if (nstep <= 0) {
      cerr << "Error: number of steps must be positive" << endl;
      return 1;
    }
    m_nstep = nstep;
    m_h = (m_b - m_a) / m_nstep;
    m_sum = 0;

    for (unsigned int i = 0; i < m_nstep; i++) {
      m_sum += f.Eval(m_a + (i + 0.5) * m_h);
      // cout << "Sum = " << m_sum << endl;
    };
    m_integral = m_sign * m_sum * m_h;
    return m_integral;
  };
};

// Attenzione a questo eh
class Simpson : public Integral {

public:
  Simpson(double a, double b) : Integral(a, b) { ; };

  inline double Integrate(unsigned int nstep, const FunzioneBase &f) override {

    if (nstep <= 0 || nstep % 2 != 0) {
      throw runtime_error("Error: number of steps must be positive and even. ");
    }

    m_nstep = nstep;
    m_h = (m_b - m_a) / m_nstep;
    m_sum = 0.;
    double c_k = 0.;
    double x = 0.;

    for (unsigned int i = 1; i < m_nstep + 1; i++) {
      c_k = 4. - 2. * (i % 2);
      x = m_a + i * m_h;
      m_sum += (c_k * f.Eval(x));
      // cout << "Step n° " << i << "; Coefficient value: " << c_k << "; x
      // value: " << x << "; Value of f in x: " << f.Eval(x) << "; Value of sum:
      // " << m_sum << endl;
    };

    m_integral = m_sign*(m_h / 3) * (f.Eval(m_a) + m_sum + f.Eval(m_b));

  //   double sum_odd = 0., sum_even = 0.;
  //   for (unsigned int i = 1; i < m_nstep; i += 2) {
  //     sum_odd += f.Eval(m_a + i*m_h);
  //   }
  //   for (unsigned int i = 2; i < m_nstep; i += 2) {
  //     sum_even += f.Eval(m_a + i * m_h);
  //   }
  //   m_integral = m_sign*(m_h/3)*(f.Eval(m_a) + 4*sum_odd + 2*sum_even + f.Eval(m_b));
    return m_integral;
  };
};

class Trapezoidi : public Integral {

public:
  Trapezoidi(double a, double b) : Integral(a, b) { ; };
  Trapezoidi(Trapezoidi &trap) : Integral(trap) {
    I_n = trap.GetIn();
    I_n1 = trap.GetIn1();
  }

  double GetIn() { return I_n; };
  double GetIn1() { return I_n1; };

  double Integrate(unsigned int nstep, const FunzioneBase &f) override {
    cerr << "Method not implemented for Trapezoidi." << endl;
    return 0;
  }
  inline double Integrate(double prec, const FunzioneBase &f) {
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
      m_h /= 2.;                      // Aggiorno h:
      m_nstep *= 2;
      // Debugging:
      // cout << "m_h = " << m_h << endl;
      // cout << "m_sum = " << m_sum << endl;
    };

    m_integral = m_sign * I_n1;
    return m_integral;
  };

private:
  double I_n, I_n1;
};

// Posizione:

class Posizione {

public:
  Posizione() {
    m_x = 0;
    m_y = 0;
    m_z = 0;
  };
  Posizione(double x, double y, double z) {
    m_x = x;
    m_y = y;
    m_z = z;
  };

  ~Posizione() { ; };

  double getX() const { return m_x; };
  double getY() const { return m_y; };
  double getZ() const { return m_z; };

  inline double getR() const {
    return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));
  };
  inline double getPhi() const { return atan2(m_x, m_y); };
  inline double getTheta() const { return atan2(m_y, m_x); };
  inline double getRho() const {
    Posizione p(m_x, m_y, m_z);
    double r = p.getR();
    if (r == 0) {
      cout << "Errore: divisione per zero nel calcolo di rho." << endl;
      exit(-1);
    }
    double cosRho = (m_z - p.getZ()) / r;
    if (cosRho > 1)
      cosRho = 1; // Da rivedere...
    if (cosRho < -1)
      cosRho = -1;
    return acos(cosRho);
  };

  Posizione &operator=(const Posizione &p) {
    if (this != &p) { // Evita l'auto-assegnazione
      this->m_x = p.getX();
      this->m_y = p.getY();
      this->m_z = p.getZ();
    }
    return *this;
  };
  inline bool operator!=(const Posizione &p) const {
    return m_x != p.getX() || m_y != p.getY() || m_z != p.getZ();
  };
  inline bool operator==(const Posizione &p) const {
    return m_x == p.getX() && m_y == p.getY() && m_z == p.getZ();
  };
  inline Posizione operator+(const Posizione &p) const {
    return Posizione(m_x + p.getX(), m_y + p.getY(), m_z + p.getZ());
  };
  inline Posizione operator-(const Posizione &p) const {
    return Posizione(m_x - p.getX(), m_y - p.getY(), m_z - p.getZ());
  };

  inline Posizione pos(const Posizione &a) const {
    Posizione diff(getX() - a.getX(), getY() - a.getY(), getZ() - a.getZ());
    return diff;
  };
  inline double getDist(const Posizione &a) const {
    Posizione diff = a - *this;
    return diff.getR();
  };

private:
  double m_x, m_y, m_z;
};

// Campo Vettoriale:

class CampoVettoriale : public Posizione {
public:
  // Costruttori
  CampoVettoriale(const Posizione &p) : Posizione(p) {
    m_Ex = 0;
    m_Ey = 0;
    m_Ez = 0;
  };
  CampoVettoriale(const Posizione &p, double Ex, double Ey, double Ez)
      : Posizione(p) {
    m_Ex = Ex;
    m_Ey = Ey;
    m_Ez = Ez;
  };
  CampoVettoriale(double x, double y, double z, double Ex, double Ey, double Ez)
      : Posizione(x, y, z) {
    m_Ex = Ex;
    m_Ey = Ey;
    m_Ez = Ez;
  };

  // Operazioni
  inline CampoVettoriale &operator+=(const CampoVettoriale &V) {
    return (*this) = (*this) + V;
  };
  inline CampoVettoriale operator+(const CampoVettoriale &V) const {
    if (V.getX() != getX() || V.getY() != getY() || V.getZ() != getZ()) {
      throw invalid_argument(
          "La somma di campi vettoriali deve avvenire nello stesso punto.");
    }

    CampoVettoriale sum(Posizione(getX(), getY(), getZ()));
    sum.setEx(getEx() + V.getEx());
    sum.setEy(getEy() + V.getEy());
    sum.setEz(getEz() + V.getEz());
    return sum;
  };

  // Metodi
  double getEx() const { return m_Ex; };
  double getEy() const { return m_Ey; };
  double getEz() const { return m_Ez; };

  void setEx(double Ex) { m_Ex = Ex; };
  void setEy(double Ey) { m_Ey = Ey; };
  void setEz(double Ez) { m_Ez = Ez; };

  inline double Modulo() {
    double modulo = pow(getEx(), 2) + pow(getEy(), 2) + pow(getEz(), 2);
    return sqrt(modulo);
  };

private:
  double m_Ex, m_Ey, m_Ez;
};

// Particella:

class Particella {
public:
  // Costrutori
  Particella() {
    m_massa = 0;
    m_carica = 0;
  };
  Particella(double m, double q) {
    m_massa = m;
    m_carica = q;
  };

  // Distruttore
  ~Particella() { ; };

  // Metodi:
  double getMassa() const { return m_massa; };
  double getCarica() const { return m_carica; };
  inline void Print() const {
    cout << "Particella: m = " << getMassa() << "\tq = " << getCarica() << endl;
  };

protected:
  double m_massa;
  double m_carica;
};

class Elettrone : public Particella {
public:
  // costruttore
  Elettrone() : Particella(9.1093826E-31, -1.60217653E-19) {;};
  // distruttore
  ~Elettrone() { ; };
  //
  inline void Print() const {
    cout << "Elettrone: m = " << m_massa << "    q = " << m_carica << endl;
  };
};

class Protone : public Particella {
public:
  // costruttore
  Protone();
  // distruttore
  ~Protone();
  //
  inline void Print() const {
    cout << "Protone: m = " << m_massa << "\tq = " << m_carica << endl;
  };
};

// Punto Materiale:

class PuntoMateriale : public Posizione, public Particella {
public:
  // Costruttori:
  PuntoMateriale(double m, double q, const Posizione &p)
      : Posizione(p), Particella(m, q) {
    ;
  };
  PuntoMateriale(double m, double q, double x, double y, double z)
      : Posizione(x, y, z), Particella(m, q) {
    ;
  };
  PuntoMateriale(const Particella &part, const Posizione &p)
      : Posizione(p), Particella(part) {
    ;
  };
  PuntoMateriale(const Particella &part, double x, double y, double z)
      : Posizione(x, y, z), Particella(part) {
    ;
  };

  // Fisica:
  inline CampoVettoriale CampoElettrico(const Posizione &p) const {

    CampoVettoriale E(p);

    const double E_0 =
        getCarica() / (4 * (M_PI) * 8.8541878188e-12 * pow(getDist(p), 2));

    E.setEx(E_0 * pos(p).getX());
    E.setEy(E_0 * pos(p).getY());
    E.setEz(E_0 * pos(p).getZ());

    return E;
  };
  inline CampoVettoriale CampoGravitazionale(const Posizione &p) const {

    CampoVettoriale G(p);

    const double G_0 = 6.67430e-11 * getMassa() / pow(getDist(p), 2);

    G.setEx(G_0 * pos(p).getX());
    G.setEy(G_0 * pos(p).getY());
    G.setEz(G_0 * pos(p).getZ());

    return G;
  };
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
  inline double Exp(double lambda) { return -(1 / lambda) * log(1 - Rand()); };
  inline double GaussBM(double avg, double sigma) {
    double s = Rand();
    double t = Rand();
    return avg + sigma * sqrt(-2 * log(1 - s)) * cos(2 * M_PI * t);
  }

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

// Stampa:
template <typename T> inline void Print(const vector<T> &a) {
  cout << "Printing vector: " << endl;
  for (auto i = a.begin(); i != a.end(); i++)
    cout << *i << " " << endl;
  cout << "Vector printed. " << endl;
};

// Funzioni vettoriali, equazioni differenziali:
class BasicVectorialFunction {

public:
  virtual vector<double> Eval(double t, const vector<double> &x) const = 0;
};

class ArmonicOscillator : public BasicVectorialFunction {

public:
  ArmonicOscillator(double omega0) { m_omega0 = omega0; };

  //  Questo Eval non dipende da t.
  //  Questa dipendenza esiste per oscillazioni forzate, ad es...
  //  Eulero è lineare in h...
  inline vector<double> Eval(double t, const vector<double> &x) const override {

    vector<double> v{x[1], -m_omega0 * m_omega0 * x[0]};
    return v;
  };

private:
  double m_omega0;
};

class Pendulum : public BasicVectorialFunction {

public:
  Pendulum(double l) { m_l = l; };

  inline vector<double> Eval(double t, const vector<double> &x) const override {
    vector<double> v{x[1], -(m_gravity / m_l) * sin(x[0])};
    return v;
  };

private:
  double m_l;
  double m_gravity = 9.08665;
};

class FDOscillator : public BasicVectorialFunction {

public:
  FDOscillator() {
    m_omega0 = 10;
    m_alpha = 1 / 30;
    m_omegaF = 5;
  };
  FDOscillator(double omega0, double alpha, double omegaF) {
    if (omega0 < alpha) {
      throw runtime_error(
          "Error: gamma must not be greater than omega. Exiting. ");
    };
    m_omega0 = omega0;
    m_alpha = alpha;
    m_omegaF = omegaF;
  };
  inline vector<double> Eval(double t, const vector<double> &x) const override {
    vector<double> v{x[1], -(m_omega0 * m_omega0 * x[0]) - (m_alpha * x[1]) +
                               sin(m_omegaF * t)};
    return v;
  };

private:
  double m_alpha;
  double m_omega0;
  double m_omegaF;
};

// Equazioni differenziali:
class BasicDifferentialEquation {

public:
  virtual vector<double> Step(double t, const vector<double> x, double h,
                              const BasicVectorialFunction &f) const = 0;
};

class Euler : public BasicDifferentialEquation {

public:
  //  Valuta l'estrapolazione dopo h. Per creare un grafico bisogna implementare
  //  un ciclo for nel main oppure definire una funzione/metodo di riempimento.
  //  Funziona per vettori di stato di taglia qualsiasi.
  inline vector<double> Step(double t, const vector<double> x, double h,
                             const BasicVectorialFunction &f) const override {
    if (sign(h) != 1) {
      cerr << "Integration step must be positive. " << endl;
      exit(-1);
    };
    vector<double> solut = x + h * f.Eval(t, x);
    return solut;
  };
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
  // Definiamo una funzione qua dentro, che calcola l'errore fino al tempo t. 
  // ATTENZIONE ALL'ERRORE SU CHE COMPONENTE!!!!
  double errorcomponent(double t, const vector<double> x, double h,
    const BasicVectorialFunction &f) {
      Runge_Kutta rk_test;
      vector<double> sol1 = x; // Passo h
      vector<double> sol2 = x; // Passo h/2
      // Attenzione al <=
      unsigned int counter1 = 0;
      unsigned int counter2 = 0;
      double t_test = 0.;
      for(t_test = 0.; t_test <= t; t_test += h / 2) {
        if(counter1%2 != 0) {
          sol1 = rk_test.Step(t_test, sol1, h, f);
        }
        sol2 = rk_test.Step(t_test, sol2, h/2, f);
        counter2++;
        counter1++;
      }
      if(counter1 != counter2) {
        sol2 = rk_test.Step(t_test, sol2, h/2, f);
        counter2++;
        sol1 = rk_test.Step(t_test, sol1, h, f);
        counter1++;
      }
      // Controllo: devono essere esattamente uno il doppio dell'altro. 
      cout << "counter1 = " << counter1 << "\tcounter2 = " << counter2 << endl;
      double delta_N = fabs(sol1[1] - sol2[1]);
      m_error = 16. / 15. * delta_N; // tanto per essere sicuri;
      return 16. / 15. * delta_N;
      // Ricorda: 15/16*k*h^4 = delta_N. 
      // QUI L'ERRORE `E SULLE Y!!!
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
                          double sup, unsigned int points, double fmax) override {
    double avg = 0.;
    double delta = 0.;
    double M2 = 0.;
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
    // In questo modo posso stimare sigma_n invogando il
    // metodo GetError() della classe madre, facendo un solo integrale.
  }
};
// Nota: per trovare sigma basta fare la deviazione standard degli M valori
// ottenuti. Si può trovare k facendo un'interpolazione con 2 o più set di
// estrazioni a N diversi. ±.

class HoMIntegrator : public IntegralMC, RandomGen {

public:
  HoMIntegrator(unsigned int seed) : IntegralMC(seed), RandomGen(seed) { ; };
  inline double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                          double sup, unsigned int points,
                          double fmax) override {
    unsigned int Hcounter = 0;
    for (unsigned int i = 0; i < points; i++) {
      double x = ran.Unif(inf, sup);
      double y = ran.Unif(0., fmax);
      if (y < f.Eval(x))
        Hcounter++;
    }
    return (sup - inf) * fmax *
           (static_cast<double>(Hcounter) / static_cast<double>(points));
  };
};

// Simulazione di apparati sperimentali:
class PrismExperiment {

public:
  PrismExperiment(unsigned int seed)
      : m_rgen(seed), m_lambda1(579.1E-9), m_lambda2(404.7E-9),
        m_alpha(60. * M_PI / 180.), m_sigmat(0.3E-3), m_A_input(2.7),
        m_B_input(60000E-18) {
    // Indici di rifrazione;
    m_n1_input = sqrt(m_A_input + m_B_input / (m_lambda1 * m_lambda1));
    m_n2_input = sqrt(m_A_input + m_B_input / (m_lambda2 * m_lambda2));

    // Imposto theta0 (completamente arbitrario);
    m_th0_input = M_PI / 2;

    // Calcolo le deviazioni minime
    // e theta1, theta2:
    m_dm1_input = 2. * asin(m_n1_input * sin(0.5 * m_alpha)) - m_alpha;
    m_th1_input = m_th0_input + m_dm1_input;
    m_dm2_input = 2. * asin(m_n2_input * sin(0.5 * m_alpha)) - m_alpha;
    m_th2_input = m_th0_input + m_dm2_input;
  }
  ~PrismExperiment() { ; };

  //  Metodi invocati in successione per ogni pseudomisura:
  inline void Execute() {
    // Ottengo le pseudomisure degli angoli.
    // Ripercorro i passi che farebbe lo sperimentatore,
    // distribuendo le costanti in modo gaussiano.

    //  "Misuro" gli angoli perturbando i dati iniziali:
    m_th0_meas = m_rgen.GaussBM(m_th0_input, m_sigmat);
    m_th1_meas = m_rgen.GaussBM(m_th1_input, m_sigmat);
    m_th2_meas = m_rgen.GaussBM(m_th2_input, m_sigmat);
  }
  inline void Analyse() {
    // Determino le deviazioni minime:
    m_dm1_meas = m_th1_meas - m_th0_meas;
    m_dm2_meas = m_th2_meas - m_th0_meas;

    // Inverto per trovare gli indici di rifrazione:
    m_n1_meas = sin((m_dm1_meas + m_alpha) / 2) / sin(m_alpha / 2);
    m_n2_meas = sin((m_dm2_meas + m_alpha) / 2) / sin(m_alpha / 2);

    // Inverto la relazione di Cauchy per trovare A e B:
    m_A_meas = (pow(m_lambda2 * m_n2_meas, 2) - pow(m_lambda1 * m_n1_meas, 2)) /
               (m_lambda2 * m_lambda2 - m_lambda1 * m_lambda1);
    m_B_meas = (pow(m_n2_meas, 2) - pow(m_n1_meas, 2)) /
               (1 / (m_lambda2 * m_lambda2) - 1 / (m_lambda1 * m_lambda1));
  }

  //  Eventuali metodi Get per accedere ai data members.
  //  Se ci serviranno li scriviamo.

  double Getth0inp() { return m_th0_input; };
  double Getth0outp() { return m_th0_meas; };

  double Getth1inp() { return m_th1_input; };
  double Getth1outp() { return m_th1_meas; };

  double Getth2inp() { return m_th2_input; };
  double Getth2outp() { return m_th2_meas; };

  double GetAinp() { return m_A_input; };
  double GetAoutp() { return m_A_meas; };

  double GetBinp() { return m_B_input; };
  double GetBoutp() { return m_B_meas; };

  double Getdm1inp() { return m_dm1_input; };
  double Getdm1outp() { return m_dm1_meas; };

  double Getdm2inp() { return m_dm2_input; };
  double Getdm2outp() { return m_dm2_meas; };

  double Getn1inp() { return m_n1_input; };
  double Getn1outp() { return m_n1_meas; };

  double Getn2inp() { return m_n2_input; };
  double Getn2outp() { return m_n2_meas; };

private:
  // generatore di numeri casuali

  RandomGen m_rgen;

  // parametri dell'apparato sperimentale:

  double m_lambda1, m_lambda2, m_alpha, m_sigmat;

  // valori delle quantita' misurabili :
  // input  : valori assunti come ipotesi nella simulazione
  // meas : valore dopo la simulazione di misura

  // A e B:
  double m_A_input, m_A_meas;
  double m_B_input, m_B_meas;

  // n1 e n2:
  double m_n1_input, m_n1_meas;
  double m_n2_input, m_n2_meas;

  // deltam1 e deltam2:
  double m_dm1_input, m_dm1_meas;
  double m_dm2_input, m_dm2_meas;

  // theta0, theta1 e theta2:
  double m_th0_input, m_th0_meas;
  double m_th1_input, m_th1_meas;
  double m_th2_input, m_th2_meas;

  //  In breve: utti i valori, "veri" e misurati,
  //  di tutte le grandezze coinvolte nel processo di simulazione.
};

// Esame 2:
// Calcola l'errore con passo h al tempo tf
inline void DoPosDistr(TH1F &hist, double h, double sigma,
                BasicDifferentialEquation &res, BasicVectorialFunction &f) {
  RandomGen ran(1);
  const double x0 = 0.;
  const double v0 = 1.;

  for (unsigned int i = 0.; i < 10000; i++) {
    // Perturba v0;
    double v = ran.GaussBM(v0, sigma);
    vector<double> x = {x0, v};
    // Calcola la soluzione a 43 secondi;
    double t = 0.;
    for (t = 0.; t <= 43; t += h) {
      x = res.Step(t, x, h, f);
    }
    // Plotta la distribuzione.
    hist.Fill(x[0]);
  }
};

inline double DevStd(double h, double sigma_v, BasicDifferentialEquation &res,
              BasicVectorialFunction &f) {
  RandomGen ran(1);
  const double x0 = 0.;
  const double v0 = 1.;

  double delta2sum = 0;
  unsigned int counter = 0;
  double avg = 0.;
  double prev_avg = 0.;

  for (unsigned int i = 0.; i < 10000; i++) {
    // Perturba v0;
    double v = ran.GaussBM(v0, sigma_v);
    vector<double> x = {x0, v};
    // Calcola la soluzione a 43 secondi;
    double t = 0.;
    for (t = 0.; t <= 43; t += h) {
      x = res.Step(t, x, h, f);
    }
    // Conta il dato di x[0] per il calcolo della varianza:
    counter++;
    prev_avg = avg;
    avg += (x[0] - avg) / counter;
    delta2sum += (x[0] - prev_avg) * (x[0] - avg);
  }
  return sqrt(delta2sum / (counter - 1));
};

// Esame 3:

inline void FillTrajectoryGraph(BasicVectorialFunction &diff, vector<double> &state,
                         Runge_Kutta &rk, double T, double h, TGraph &graph) {
  unsigned int counter = 0;
  for (double t = 0.; t <= 10 * T; t += h) {
    // Disegna il punto (x, y) in R^2
    graph.SetPoint(counter, state[0], state[1]);
    // E evolvi:
    state = rk.Step(t, state, h, diff);
    counter++;
  }
};

inline double radius(double x, double y) { return (x * x + y * y); }

inline void FillRadiusGraph(BasicVectorialFunction &diff, vector<double> &state,
                     Runge_Kutta &rk, double T, double h, TGraph &graph) {
  unsigned int counter = 0;
  double r = 0.;
  for (double t = 0.; t <= 10 * T; t += h) {
    // Disegna il raggio nel suo grafico:
    r = radius(state[0], state[1]);
    graph.SetPoint(counter, t, r);
    // E evolvi:
    state = rk.Step(t, state, h, diff);
    counter++;
  }
};

inline void draw(TGraph &graph, TCanvas &c) {
  c.cd();
  graph.SetMarkerStyle(20);
  graph.SetMarkerColor(kBlue);
  graph.SetMarkerSize(0.2);
  graph.GetXaxis()->SetTitle("x [m]");
  graph.GetYaxis()->SetTitle("y [m]");
  graph.Draw("APL");
};

// Esame 4:

class function4 : public FunzioneBase {
public:
  function4() { ; };
  function4(const double lambda, const double L, const double d, double x) {
    m_lambda = lambda;
    m_d = d;
    m_L = L;
    m_x = x;
  };
  // Copy constructor:
  function4(function4 &f) {
    m_lambda = f.GetLambda();
    m_L = f.GetL();
    m_d = f.Getd();
    m_x = f.Getx();
  }
  double GetLambda() { return m_L; };
  double GetL() { return m_L; };
  double Getd() { return m_d; };
  double Getx() { return m_x; };
  void Setx(double x) { m_x = x; };
  void SetLambda(double lambda) { m_lambda = lambda; };
  inline double Eval(double t) const override {
    return (1 / m_d) * cos((sqrt(m_L * m_L + (m_x - t) * (m_x - t)) -
                            sqrt(m_L * m_L + m_x * m_x)) /
                           m_lambda);
  }

private:
  double m_lambda, m_L, m_d, m_x;
  // Potremmo dichiararli const qui e non nel main per poi costruire...
};

class Ampiezza : public FunzioneBase {
public:
  // Costruttore:
  Ampiezza(function4 &f4, Trapezoidi &trap) : m_f(f4), m_trap(trap) {};
  // Restituisce l'integrale di m_f:
  inline double Eval(double x) const override {
    m_f.Setx(x);
    return m_trap.Integrate(1E-04, m_f);
  };

private:
  function4 &m_f;
  Trapezoidi &m_trap;
};

inline void graph4Fill(Trapezoidi &trap, function4 &f4, TGraph &graph) {
  unsigned int counter = 0;
  double I = 0.;
  for (double x = -0.1; x <= 0.1; x += 0.002) {
    f4.Setx(x);
    // Integra la funzione rispetto a t e stampa x, A(x):
    I = trap.Integrate(1E-04, f4);
    // cout << x << "\t" << I << endl;
    graph.SetPoint(counter, x, I);
    // cout << "counter    " << counter << endl;
    counter++;
  }
};

inline void first_zero(const double lambda, function4 f4, Trapezoidi trap) {
  f4.SetLambda(lambda);
  const double prec = 1E-06;
  Bisezione mybisect(prec);
  double b = 0.05; // Esempio
  unsigned int nmax = 500;
  Ampiezza myA(f4, trap);

  double zero = mybisect.CercaZeriRef(0., b, myA, prec, nmax);

  cout << "Lambda = " << lambda << " m;   \tZero: x = " << zero << " m" << endl;
};

// Esame 5;

// Usiamo gli errori relativi;
// Parametri iniziali: V0, V1, R, C

class CapacitorExperiment {

public:
  CapacitorExperiment(unsigned int seed, double sigmarelV0, double sigmarelV1,
                      double sigmarelR, double sigmarelt)
      : m_rgen(seed), m_sigmarel_R(sigmarelR), m_sigmarel_V0(sigmarelV0),
        m_sigmarel_V1(sigmarelV1), m_sigmarel_t(sigmarelt), m_C_input(2E-06),
        m_R_input(100E03), m_V0_input(12.), m_V1_input(3.),
        m_t_input(m_C_input * m_R_input * log(m_V0_input / m_V1_input)) {
    ;
  }
  // Per il tempo ho messo un valore a caso, non serve per l'errore.

  // Così posso fare l'esperimento con una sola delle fonti di errore
  // e vedere quale influisce di più!
  ~CapacitorExperiment() { ; };

  //  Metodi invocati in successione per ogni pseudomisura:
  inline void Execute() {
    // Esame:
    // Perturba gaussianamente i valori di input e genera il tempo che serve
    // per ottenere un valore plausibile di C;
    m_R_meas = m_rgen.GaussBM(m_R_input, m_R_input * m_sigmarel_R);
    m_V0_meas = m_rgen.GaussBM(m_V0_input, m_V0_input * m_sigmarel_V0);
    m_V1_meas = m_rgen.GaussBM(m_V1_input, m_V1_input * m_sigmarel_V1);
    m_t_meas = m_rgen.GaussBM(m_t_input, m_t_input * m_sigmarel_t);
  }
  inline void Analyse() {
    // Esame: calcola C con i valori di output.
    m_C_meas = m_t_meas / (m_R_meas * log(m_V0_meas / m_V1_meas));
  }

  //  Eventuali metodi Get per accedere ai data members.
  //  Se ci serviranno li scriviamo.

  double GetV0inp() { return m_V0_input; };
  double GetV0outp() { return m_V0_meas; };

  double GetV1inp() { return m_V1_input; };
  double GetV1outp() { return m_V1_meas; };

  double GetRinp() { return m_R_input; };
  double GetRoutp() { return m_R_meas; };

  double GetCinp() { return m_C_input; };
  double GetCoutp() { return m_C_meas; };
  double GetSigmaCoutp() { return m_sigmarel_C; };

private:
  // generatore di numeri casuali

  RandomGen m_rgen;

  // parametri dell'apparato sperimentale:

  // valori delle quantita' misurabili :
  // input  : valori assunti come ipotesi nella simulazione
  // meas : valore dopo la simulazione di misura

  // V0, V1:
  double m_V0_input, m_V0_meas;
  double m_sigmarel_V0;
  double m_V1_input, m_V1_meas;
  double m_sigmarel_V1;

  // R:
  double m_R_input, m_R_meas;
  double m_sigmarel_R;

  // C:
  double m_C_input, m_C_meas;
  double m_sigmarel_C;

  // t;
  double m_t_input, m_t_meas;
  double m_sigmarel_t;

  //  In breve: utti i valori, "veri" e misurati,
  //  di tutte le grandezze coinvolte nel processo di simulazione.
};

struct ExpResults {
  double C_final, sigmarel_C;
};

// Per trovare i contributi all'errore
inline ExpResults RunExperiment(unsigned int N, unsigned int seed, double sigmarelV0,
                         double sigmarelV1, double sigmarelR,
                         double sigmarelt) {
  CapacitorExperiment exp(seed, sigmarelV0, sigmarelV1, sigmarelR, sigmarelt);
  vector<double> C_out(N);
  for (unsigned int i = 0; i < N; i++) {
    exp.Execute();
    exp.Analyse();
    // Salva i valori per calcolare
    // media e dev std:
    C_out[i] = exp.GetCoutp();
  }
  double C_final = DoAvg(C_out);
  double err_rel_C = DevStd(C_out) / C_final;
  ExpResults results{C_final, err_rel_C};
  return results;
};

inline void FillWithError(vector<double> &sigV, vector<double> &sigC,
                   const double err) {
  cout << "Errore % su V: \tErrore % su C:" << endl;
  for (double sigmaV = 0.02; sigmaV <= 0.07; sigmaV += 0.005) {
    ExpResults res3 = RunExperiment(10000, 1, sigmaV, sigmaV, err, err);
    sigC.push_back(res3.sigmarel_C);
    sigV.push_back(sigmaV);
    cout << sigmaV * 100 << "%\t\t" << res3.sigmarel_C * 100 << "%" << endl;
  }
};

// Esame Zanoli:
const double unitm = 1.66E-27;
const double unitq = 1.60E-19;

struct Particle {
  double carica, massa;
};

class DiffEqZ : public BasicVectorialFunction {
public:
  DiffEqZ(double r, double B0, Particle p) {
    m_r = r;
    m_B0z = B0;
    m_p.carica = p.carica;
    m_p.massa = p.massa;
  }
  inline vector<double> Eval(double t, const vector<double>& state) const override {
    double coeff = 0.;
    if (m_r != 0) coeff = (m_p.carica / m_p.massa) * (m_B0z * pow(M_E, - m_r));
    else coeff = (m_p.carica / m_p.massa) * (m_B0z);
    vector<double> stateeval {state[2], state[3], - coeff * state[3], coeff * state[2]};
    return stateeval;
  }
private:
  double m_r, m_B0z;
  Particle m_p;
};

struct diam {
  double diam, t;
};

inline diam Diameter(double R, double B0, Particle p, vector<double>state, double h) {
  DiffEqZ myeq (R, B0, p);
  Runge_Kutta myrk;double t = 0.;
  while(state[0] >= 0.) {
      state = myrk.Step(t, state, h, myeq);
      // cout << state[0] << endl;
      t += h;
  }
  diam d;
  d.diam = state[1];
  d.t = t;
  return d;
}

struct Results {
  double davg1, davg2;
  double sig1, sig2;
};

inline Results SimulateTraj(const double B0, TH1F& hist1, TH1F& hist2, const double v0x, 
  const double sigrel, double h, Particle p1, Particle p2) {
  unsigned int N = 5000;
  double d1test = 0.;
  double d2test = 0.;
  double v0x1 = v0x;
  RandomGen myran(1);
  for(unsigned int i = 0; i < N; i++) {
      // Fai evolvere e metti nei grafici;
      v0x1 = myran.GaussBM(v0x, v0x * sigrel);
      vector<double> state1 {0., 0., v0x1, 0.};
      vector<double> state2 {0., 0., v0x1, 0.};
      d1test = Diameter(1., B0, p1, state1, h).diam;
      d2test = Diameter(1., B0, p2, state1, h).diam;
      hist1.Fill(d1test);
      hist2.Fill(d2test);
  }
  Results rel;
  rel.davg1 = hist1.GetMean();
  rel.davg2 = hist2.GetMean();
  rel.sig1 = hist1.GetStdDev();
  rel.sig2 = hist2.GetStdDev();
  return rel;
}

inline void draw(TH1F& newgraph, TCanvas& c3, string Title) {

  c3.cd();
  // newgraph.SetMarkerSize(0.5);
  // newgraph.SetMarkerStyle(20);
  newgraph.SetLineColor(kBlue);
  newgraph.SetFillColor(0);

  newgraph.SetTitle(Title.c_str());

  newgraph.Draw("HIST");
  c3.Update();
}
