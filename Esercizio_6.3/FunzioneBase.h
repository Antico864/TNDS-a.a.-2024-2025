#pragma once

#include <iostream>
#include <limits>
#include <tgmath.h>

using namespace std;

class FunzioneBase {

public:
  virtual double Eval(double x) const = 0;
  virtual ~FunzioneBase() { ; };
};

class Tan_men_x : public FunzioneBase {
public:
  Tan_men_x() { ; };

  double Eval(double x) const override {
    return tan(x) - x;
  }
};
