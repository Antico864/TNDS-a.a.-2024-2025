#include "Vettore.h"
#include <iostream>

using namespace std;

Vettore::Vettore() {
  m_N = 0;
  m_v = NULL;
}

Vettore::Vettore(int N) {
  m_N = N;
  m_v = new double[N];
  for (int k = 0; k < N; k++) {
    m_v[k] = 0;
  }
}

Vettore::~Vettore() { delete[] m_v; }

void Vettore::SetComponent(int k, double val) { m_v[k] = val; }

double Vettore::GetComponent(int k) const {
  if (k > m_N || k < 0) {
    string s = "Dimensione non valida: indice = " + to_string(k) +
               ", dimensione = " + to_string(m_N);
    __throw_runtime_error(s.c_str());
  }
  return m_v[k];
}

int Vettore::GetSize() {
  int size = 0;
  for (int k = 0; k < m_N; k++) {
    if (m_v[k] != 0) {
      size++;
    }
  }
  return size;
}

void Vettore::Scambia(int primo, int secondo) {
  double temp = GetComponent(primo);
  SetComponent(primo, GetComponent(secondo));
  SetComponent(secondo, temp);
}

Vettore::Vettore(const Vettore &A) {
  m_N = A.GetN();
  m_v = new double[m_N];
  for (int i = 0; i < m_N; i++) {
    m_v[i] = A.GetComponent(i);
  }
}

Vettore &Vettore::operator=(const Vettore &A) {
  m_N = A.GetN();
  if (m_v) {
    delete[] m_v;
  };
  m_v = new double[m_N];
  for (int i = 0; i < m_N; i++)
    m_v[i] = A.GetComponent(i);
  return *this;
}

double &Vettore::operator[](int k) {
  if (k >= m_N || k < 0) {
    string s =
        "Errore: indice = " + to_string(k) + ", dimensione = " + to_string(m_N);
    __throw_runtime_error(s.c_str());
  }
  return m_v[k];
}
