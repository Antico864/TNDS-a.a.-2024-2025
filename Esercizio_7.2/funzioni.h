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

class xsinx : public FunzioneBase {

public:
  xsinx() {;};

  double Eval(double x) const override { return x * sin(x); };
};
