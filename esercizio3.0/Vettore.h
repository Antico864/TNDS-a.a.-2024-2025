#pragma once

#include <fstream>
#include <iostream>

using namespace std;

template <typename T> class Vettore {
public:
  Vettore();

  Vettore(int);

  Vettore(const Vettore &);

  T &operator=(const Vettore &);

  ~Vettore();

  unsigned int GetN() const { return m_N; };

  void SetComponent(int k, T val);

  T GetComponent(int k) const;

  void Scambia(int primo, int secondo);

  unsigned int GetSize();

  T &operator[](int);

private:
  unsigned int m_N;
  T *m_v;
};

template <typename T> Vettore<T>::Vettore() {
  m_N = 0;
  m_v = NULL;
}

template <typename T> Vettore<T>::Vettore(int N) {
  m_N = N;
  m_v = new double[N];
  for (int k = 0; k < N; k++) {
    m_v[k] = 0;
  }
}

template <typename T> Vettore<T>::Vettore(const Vettore &A) {
  m_N = A.GetN();
  m_v = new double[m_N];
  for (int i = 0; i < m_N; i++) {
    m_v[i] = A.GetComponent(i);
  }
}

template <typename T> T &Vettore<T>::operator=(const Vettore &A) {
  m_N = A.GetN();
  if (m_v) {
    delete[] m_v;
  };
  m_v = new double[m_N];
  for (int i = 0; i < m_N; i++)
    m_v[i] = A.GetComponent(i);
  return *this;
}

template <typename T> Vettore<T>::~Vettore() { delete[] m_v; }

template <typename T> void Vettore<T>::SetComponent(int k, T val) {
  m_v[k] = val;
}

template <typename T> T Vettore<T>::GetComponent(int k) const {
  if (k > m_N || k < 0) {
    cout << "Dimensione non valida: indice = " << k << ", dimensione = " << m_N
         << endl;
    exit(2);
  } else {
    return m_v[k];
  }
}

template <typename T> void Vettore<T>::Scambia(int primo, int secondo) {
  double temp = GetComponent(primo);
  SetComponent(primo, GetComponent(secondo));
  SetComponent(secondo, temp);
}

template <typename T> unsigned int Vettore<T>::GetSize() {
  int size = 0;
  for (int k = 0; k < m_N; k++) {
    while (m_v[k] != 0) {
      size = size++;
    }
  }
  return size;
}

template <typename T> T &Vettore<T>::operator[](int i) {
  if (0 < i < m_N)
    return m_v[i];
  else {
    cout << "Errore: indice = " << i << ", dimensione = " << m_N << endl;
    exit(3);
  }
}
