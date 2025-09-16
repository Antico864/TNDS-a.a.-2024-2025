#pragma once

#include <iostream>

using namespace std;

class Vettore {
public:
  Vettore();
  Vettore(int N);
  ~Vettore();

  void SetComponent(int, double);
  double GetComponent(int) const;
  int GetSize();
  unsigned int GetN() const { return m_N; };

  void Scambia(int, int);

  Vettore(const Vettore &);

  Vettore &operator=(const Vettore &);
  double &operator[](int);

private:
  unsigned int m_N;
  double *m_v;
};
