#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

//======================================================================
//  Somma di due vector componente per componente:
//======================================================================

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

//======================================================================
//  Differenza di due vector componente per componente:
//======================================================================

template <typename T>
inline vector<T> operator-(const vector<T> &a, const vector<T> &b) {
  if (a.size() != b.size()) {
    cerr << "Vector sizes in operator- must be equal. " << endl;
    exit(-1);
  }
  vector<T> diff(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    diff[i] = a[i] - b[i];
  }
  return diff;
};

//======================================================================
//  Prodotto scalare canonico tra due vector:
//======================================================================

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

//======================================================================
//  Prodotto scalare*vector:
//======================================================================

template <typename T>
inline vector<T> operator*(const T k, const vector<T> &a) {
  vector<T> prod(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    prod[i] = k * a[i];
  }
  return prod;
};

//======================================================================
//  Prodotto vector*scalare:
//======================================================================

template <typename T>
inline vector<T> operator*(const vector<T> &a, const T k) {
  vector<T> prod(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    prod[i] = k * a[i];
  }
  return prod;
};

// ===============================================================================
// Divisione tra vector/scalare:
// ===============================================================================

template <typename T>
inline vector<T> operator/(const vector<T> &a, const T k) {
  vector<T> div(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    div[i] = (1 / k) * a[i];
  }
  return div;
};

// ===============================================================================
// Somma a+b con risultato messo in a:
// ===============================================================================

template <typename T>
inline vector<T> operator+(vector<T> &a, const vector<T> &b) {
  if (a.size() != b.size()) {
    cerr << "Vector sizes in operator+ must be equal. " << endl;
    exit(-1);
  }
  vector<T> sum(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    a[i] += b[i];
  }
  return a;
};

// ===============================================================================
// Differenza a-b con risultato messo in a:
// ===============================================================================

template <typename T>
inline vector<T> operator-(vector<T> &a, const vector<T> &b) {
  if (a.size() != b.size()) {
    cerr << "Vector sizes in operator+ must be equal. " << endl;
    exit(-1);
  }
  vector<T> sum(a.size());
  for (unsigned int i = 0; i < static_cast<int>(a.size()); i++) {
    a[i] -= b[i];
  }
  return a;
};

// ===============================================================================
// Stampa di un vector:
// ===============================================================================

template <typename T> inline void Print(const vector<T> &a) {
  cout << "Printing vector: " << endl;
  for (auto i = a.begin(); i != a.end(); i++)
    cout << *i << " " << endl;
  cout << "Vector printed. " << endl;
};
