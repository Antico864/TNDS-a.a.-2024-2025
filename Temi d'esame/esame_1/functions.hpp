#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <tgmath.h>

using namespace std;

class BasicFunction {

public:
  virtual double Eval(double x) const = 0;
  virtual ~BasicFunction() { ; };
};

class firstfunction : public BasicFunction {
public:
  firstfunction(double exp1, double exp2) {
    m_exp1 = exp1;
    m_exp2 = exp2;
  };
  ~firstfunction() { ; };
  virtual double Eval(double x) const override;

  // Non sono necessari data membri privati...
private:
  double m_exp1, m_exp2;
};

class secondfunction : public BasicFunction {
public:
  secondfunction() { ; };
  ~secondfunction() { ; };
  virtual double Eval(double x) const override;
};

double Var(vector<double> v);

template <typename T> double DoAvg(vector<T> a) {
  double avg = 0;
  for (int k = 0; k < a.size(); k++) {
    avg = (avg * (k) + a[k + 1]) / (k + 1);
  }
  return avg;
};
