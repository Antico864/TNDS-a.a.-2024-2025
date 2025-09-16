#pragma once

#include <iostream>

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
