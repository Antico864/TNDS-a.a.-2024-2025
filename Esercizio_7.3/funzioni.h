#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <tgmath.h>

using namespace std;

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
  virtual double Eval(double x) const override {
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

class FunzioneTrigonometrica : public FunzioneBase {

public:
  void SetA(double A) { m_A = A; };
  void SetOmega(double omega) { m_omega = omega; };
  void SetPhi(double phi) { m_phi = phi; };

  double GetA() { return m_A; };
  double GetOmega() { return m_omega; };
  double GetPhi() { return m_phi; };

  double GetPeriod() {
    double T = (2 * M_PI) / m_omega;
    return T;
  };

  virtual double EvalSin(double x) const {
    return m_A * sin(m_omega * x + m_phi);
  };
  virtual double EvalCos(double x) const {
    return m_A * cos(m_omega * x + m_phi);
  };

protected:
  double m_A, m_omega, m_phi;
};

class Tan_men_x : public FunzioneTrigonometrica {
public:
  Tan_men_x() {
    m_A = 1.;
    m_omega = 1.;
    m_phi = 0.;
  }

  virtual double Eval(double x) const override {
    return EvalSin(x) - x * EvalCos(x);
  }
};

class xsinx : public FunzioneTrigonometrica {

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

  //    virtual double Eval(double x) const override { return x*EvalSin(x); };
  virtual double Eval(double x) const override { return x * sin(x); };
};

class NormGaussian : public FunzioneBase {

public:
  NormGaussian() {
    m_mu = 0;
    m_sigma = 1;
  };
  NormGaussian(double mu, double sigma) {
    m_mu = mu;
    m_sigma = sigma;
  };
  ~NormGaussian() { ; };

  double GetMu() const { return m_mu; };
  double GetSigma() const { return m_sigma; };

  virtual double Eval(double x) const override {
    return (1 / sqrt(2.0 * M_PI)) * pow(M_E, ((-x * x) / 2.0));
  };

private:
  double m_mu, m_sigma;
};

int getsignfig(const double number);
