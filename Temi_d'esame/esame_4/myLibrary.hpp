#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

// Conversione da double a stringa:
string convert(double h) {
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

template <typename T> double DevStd(vector<T> a, double avg) {
    double delta2 = 0;
    double delta2sum = 0;
    unsigned int counter = 0;
    for (unsigned int k = 0; k < a.size(); k = k + 7) {
        delta2 = pow(a[k] - avg, 2);
        delta2sum = delta2sum + delta2;
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

inline double sign(double x) { return (x == 0. ? 0. : (x > 0 ? 1. : -1));};

int getsignfig(const double number) {
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
        virtual ~FunzioneBase() {;};
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
        ~Parabola() {;};
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

        double GetVertex() const { return -m_b/(2 * m_a); };

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

        inline double Eval(double x) const override {
            return sin(x) - x * cos(x);
        }
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

        inline double Eval(double x) const override { return x*sin(x); };
    private:
        double m_A, m_omega, m_phi;
};





// const inline double EvalA(double x, function4 & f4, Trapezoidi & trap) {
//     // Ritorna l'integrale in funzione di x!!!!!
//     f4.Setx(x);
//     return trap.Integrate(1., f4);
// };

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

        Bisezione() {;};
        Bisezione(double prec) { m_prec = prec; };
        Bisezione(unsigned int n) { m_nmax = n; };
        virtual ~Bisezione() {;};

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
                                                if (m_niterations > m_nmax) break;
                                                m_niterations++;
                                                if (sign(f_a) * sign(f_c) <= 0) {
                                                    m_b = c;
                                                    f_b = f_c;
                                                } else if (sign(f_b) * sign(f_c) <= 0) {
                                                    m_a = c;
                                                    f_a = f_c;
                                                } else return NAN;
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
    
        inline virtual double Integrate(unsigned int nstep, const FunzioneBase &) = 0;
    
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
        Midpoint(double a, double b) : Integral(a, b) {;};

        inline virtual double Integrate(unsigned int nstep, const FunzioneBase &f) override {
            if (nstep <= 0) {
                cerr << "Error: number of steps must be positive" << endl;
                return 1;
            }
            m_nstep = nstep;
            m_h = (m_b - m_a) / m_nstep;
            m_sum = 0;

            for (unsigned int i = 0; i < m_nstep; i++) {
                m_sum += f.Eval(m_a + (i + 0.5)*m_h);
                // cout << "Sum = " << m_sum << endl;
            };
            m_integral = m_sign * m_sum * m_h;
            return m_integral;
        };
};

class Simpson : public Integral {

    public:
        Simpson(double a, double b) : Integral(a, b) { ; };
    
        inline virtual double Integrate(unsigned int nstep, const FunzioneBase &f) {

            if (nstep <= 0 || nstep % 2 != 0) {
                __throw_runtime_error("Error: number of steps must be positive and even. ");
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
                // cout << "Step n° " << i << "; Coefficient value: " << c_k << "; x value:
                // " << x << "; Value of f in x: " << f.Eval(x) << "; Value of sum: " <<
                // m_sum << endl;
            };
          
            m_integral = (m_h / 3) * (f.Eval(m_a) + m_sum + f.Eval(m_b));
            return m_integral;
        };
};

class Trapezoidi : public Integral {

    public:
      Trapezoidi(double a, double b) : Integral(a, b) { ; };
      Trapezoidi(Trapezoidi& trap) : Integral(trap) {
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
          m_h /= 2.; // Aggiorno h:
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

        ~Posizione() {;};

        inline double getX() const { return m_x; };
        inline double getY() const { return m_y; };
        inline double getZ() const { return m_z; };

        inline double getR() const {
            return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));
        };
        inline double getPhi() const {return atan2(m_x, m_y);};
        inline double getTheta() const {return atan2(m_y, m_x);};
        inline double getRho() const {
            Posizione p(m_x, m_y, m_z);
            double r = p.getR();
            if (r == 0) {
              cout << "Errore: divisione per zero nel calcolo di rho." << endl;
              exit(-1);
            }
            double cosRho = (m_z - p.getZ()) / r;
            if (cosRho > 1)
              cosRho = 1; // Correzione numerica
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
        Posizione operator+(const Posizione &p) const {
            return Posizione(m_x + p.getX(), m_y + p.getY(), m_z + p.getZ());
        };
        Posizione operator-(const Posizione &p) const {
            return Posizione(m_x - p.getX(), m_y - p.getY(), m_z - p.getZ());
        };

        Posizione pos(const Posizione & a) const {
            Posizione diff(getX() - a.getX(), getY() - a.getY(), getZ() - a.getZ());
            return diff;
        };
        inline double getDist(const Posizione & a) const {
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
        CampoVettoriale(const Posizione &p, double Ex, double Ey,
            double Ez): Posizione(p) {
                                    m_Ex = Ex;
                                    m_Ey = Ey;
                                    m_Ez = Ez;
                                };
                                CampoVettoriale(double x, double y, double z, double Ex,
                                    double Ey, double Ez) : Posizione(x, y, z) {
                                        m_Ex = Ex;
                                        m_Ey = Ey;
                                        m_Ez = Ez;
                                    };

        // Operazioni
        CampoVettoriale &operator+=(const CampoVettoriale &V) {
            return (*this) = (*this) + V;
        };
        CampoVettoriale operator+(const CampoVettoriale &V) const {
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
        inline double getEx() const { return m_Ex; };
        inline double getEy() const { return m_Ey; };
        inline double getEz() const { return m_Ez; };

        inline void setEx(double Ex) { m_Ex = Ex; };
        inline void setEy(double Ey) { m_Ey = Ey; };
        inline void setEz(double Ez) { m_Ez = Ez; };

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
        ~Particella() {;};

        // Metodi:
        inline double getMassa() const { return m_massa; };
        inline double getCarica() const { return m_carica; };
        void Print() const {
            cout << "Particella: m = " << getMassa() << "\tq = " << getCarica() << endl;
        };

    protected:
        double m_massa;
        double m_carica;
};

class Elettrone : public Particella {
    public:
        // costruttore
        Elettrone() : Particella(9.1093826E-31, -1.60217653E-19) {};
        // distruttore
        ~Elettrone() {;};
        //
        void Print() const {
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
        void Print() const {
            cout << "Protone: m = " << m_massa << "\tq = " << m_carica << endl;
        };
};

// Punto Materiale:

class PuntoMateriale : public Posizione, public Particella {
    public:
        // Costruttori:
        PuntoMateriale(double m, double q, const Posizione &p)
        : Posizione(p), Particella(m, q) {;};
        PuntoMateriale(double m, double q, double x, double y, double z)
        : Posizione(x, y, z), Particella(m, q) {;};
        PuntoMateriale(const Particella &part, const Posizione &p)
        : Posizione(p), Particella(part) {;};
        PuntoMateriale(const Particella &part, double x, double y, double z)
        : Posizione(x, y, z), Particella(part) {;};

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
  };
  vector<T> diff(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    diff[i] = a[i] - b[i];
  };
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
        inline vector<double> Eval(double t, const vector<double> &x) const override{

            vector<double> v{x[1], -m_omega0 * m_omega0 * x[0]};
            return v;
        };

    private:
        double m_omega0;
};

class Pendulum : public BasicVectorialFunction {

    public:
      Pendulum(double l) { m_l = l; };
    
      inline vector<double> Eval(double t, const vector<double> &x) const override{
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
      virtual vector<double> Step(double t, const vector<double> x, double h, const BasicVectorialFunction &f) 
            const {
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
      inline double Integrate(RandomGen &ran, const FunzioneBase &f,
        double inf, double sup, unsigned int points,
        double fmax) {
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
                        return (sup - inf)*avg;
                        // In questo modo posso stimare sigma_n invogando il 
                        // metodo GetError() della classe madre, facendo un solo integrale. 
                    }
};

class HoMIntegrator : public IntegralMC, RandomGen {
    
    public:
      HoMIntegrator(unsigned int seed) : IntegralMC(seed), RandomGen(seed) { ; };
      inline double Integrate(RandomGen &ran, const FunzioneBase &f, double inf,
                               double sup, unsigned int points, double fmax) override {
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
    : m_rgen(seed), m_lambda1(579.1E-9), 
    m_lambda2(404.7E-9),
    m_alpha(60. * M_PI / 180.), 
    m_sigmat(0.3E-3), 
    m_A_input(2.7),
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
      void Execute() {
        // Ottengo le pseudomisure degli angoli.
        // Ripercorro i passi che farebbe lo sperimentatore,
        // distribuendo le costanti in modo gaussiano.
      
        //  "Misuro" gli angoli perturbando i dati iniziali:
        m_th0_meas = m_rgen.GaussBM(m_th0_input, m_sigmat);
        m_th1_meas = m_rgen.GaussBM(m_th1_input, m_sigmat);
        m_th2_meas = m_rgen.GaussBM(m_th2_input, m_sigmat);
      }
      void Analyse() {
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













// Esame 4:

class function4 : public FunzioneBase {
    public: 
        function4() {;};
        function4(const double lambda, 
                  const double L, 
                  const double d, 
                  double x) 
                            {
                                m_lambda = lambda;
                                m_d = d;
                                m_L = L;
                                m_x = x;
                            };
        // Copy constructor:
        function4(function4& f) {
            m_lambda = f.GetLambda();
            m_L = f.GetL();
            m_d = f.Getd();
            m_x = f.Getx();
        }
        double GetLambda() {return m_L;};
        double GetL() {return m_L;};
        double Getd() {return m_d;};
        double Getx() {return m_x;};
        void Setx(double x) {m_x = x;};
        void SetLambda(double lambda) {m_lambda = lambda;};
        inline double Eval(double t) const override {
            return (1/m_d) * cos((sqrt(m_L * m_L + (m_x - t) * (m_x - t)) - sqrt(m_L * m_L + m_x * m_x)) / m_lambda);
        }
    private: 
        double m_lambda, m_L, m_d, m_x;
        // Potremmo dichiararli const qui e non nel main per poi costruire...
};

class Ampiezza : public FunzioneBase {
    public:
        // Costruttore:
        Ampiezza (function4 & f4, Trapezoidi & trap) : 
            m_f(f4), m_trap(trap) {};
        // Restituisce l'integrale di m_f:
        inline double Eval(double x) const override {
            m_f.Setx(x);
            return m_trap.Integrate(1E-04, m_f);
        };
    private: 
        function4 & m_f;
        Trapezoidi & m_trap;
};

void graph4Fill(Trapezoidi & trap, function4 & f4, TGraph & graph) {
    unsigned int counter = 0;
    double I = 0.;
    for(double x = -0.1; x <= 0.1; x += 0.002) {
        f4.Setx(x);
        // Integra la funzione rispetto a t e stampa x, A(x): 
        I = trap.Integrate(1E-04, f4);
        // cout << x << "\t" << I << endl;
        graph.SetPoint(counter, x, I);
        // cout << "counter    " << counter << endl;
        counter++;
    }
};

void first_zero(const double lambda, function4 f4, Trapezoidi trap) {
    f4.SetLambda(lambda);
    const double prec = 1E-06;
    Bisezione mybisect (prec);
    double b = 0.05; // Esempio
    unsigned int nmax = 500;
    Ampiezza myA (f4, trap);

    double zero = mybisect.CercaZeriRef(0., b, myA, prec, nmax);

    cout << "Lambda = " << lambda << " m;   \tZero: x = " << zero << " m" << endl;
};


