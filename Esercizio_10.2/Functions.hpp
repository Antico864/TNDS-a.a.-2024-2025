#pragma once

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

//=============================================
//  Classi:
//=============================================

class FunzioneBase {

public:
  virtual double Eval(double x) const = 0;
};

//  Questa è la funzione che dobbiamo integrare:
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

  virtual double Eval(double x) const override { return x * sin(x); };

private:
  double m_A, m_omega, m_phi;
};

//=============================================
//  Altre funzioni:
//=============================================

double Avg(const vector<double> &v);

double Var(vector<double> v);

template <typename T> vector<T> Read(const string filename) {
  vector<T> v;
  ifstream inputFile(filename);
  if (!inputFile) {
    throw ios_base::failure("Cannot open file " + filename);
  } else {
    while (!inputFile.eof()) {
      T val;
      inputFile >> val;
      v.push_back(val);
    };
  };
  inputFile.close();
  return v;
};

template <typename T> void Print(const string outputFilename, vector<T> a) {
  ofstream outputFile(outputFilename);
  if (!outputFile) {
    __throw_runtime_error("Error opening output file: exiting ");
  }
  unsigned int mN = a.size();
  for (unsigned int i = 0; i < mN - 1; i++) {
    outputFile << a[i] << endl;
  }
  outputFile << a[mN - 1];
  outputFile.close();
};
