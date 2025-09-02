#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

template <typename T> vector<T> Read(char *filename) {
  vector<T> v;
  ifstream inputFile(filename);
  if (!inputFile) {
    cout << "Cannot open file " << filename << endl;
    exit(11);
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

template <typename T> double CalcAvg(vector<T> &a) {
  double sum = 0;
  int nM = a.size();
  for (unsigned int k = 0; k < nM; k++) {
    sum += a[k];
  }
  return sum / nM;
};

template <typename T> double CalcVar(vector<T> &a) {
  double var = 0;
  int mN = a.size();
  double avg = CalcAvg<double>(a);
  for (unsigned int k = 0; k < mN; k++) {
    var += pow(a[k] - avg, 2);
  }
  return var / (a.size() - 1);
};

template <typename T> double CalcMed(vector<T> &a) {
  double mediana = 0;
  int mN = a.size();
  if (mN % 2 != 0) {
    mediana = a[(mN + 1) / 2];
  } else {
    mediana = abs((a[mN / 2 + 1] - a[mN / 2]) / 2);
  }
  return mediana;
};

template <typename T> void Merge(vector<T> &v, int left, int right) {
  int mid = left + (right - left) / 2;
  int n1 = mid - left + 1;
  int n2 = right - mid;

  vector<T> leftVect(n1);
  vector<T> rightVect(n2);

  for (int i = 0; i < n1; i++)
    leftVect[i] = v[left + i];
  for (int j = 0; j < n2; j++)
    rightVect[j] = v[mid + 1 + j];

  int i = 0, j = 0, k = left;

  // Fusione ordinata
  while (i < n1 && j < n2) {
    if (leftVect[i] <= rightVect[j]) {
      v[k] = leftVect[i];
      i++;
    } else {
      v[k] = rightVect[j];
      j++;
    }
    k++;
  }

  // Copia degli elementi rimanenti
  while (i < n1) {
    v[k] = leftVect[i];
    i++;
    k++;
  }

  while (j < n2) {
    v[k] = rightVect[j];
    j++;
    k++;
  }
};

// DA AGGIUSTARE!!!!!
template <typename T> void MergeSort(vector<T> &v, int left, int right) {
  if (left < right) {
    int mid = left + (right - left) / 2;

    MergeSort(v, left, mid);
    MergeSort(v, mid + 1, right);

    Merge(v, left, right);
  }
};

template <typename T> void Print(const char *outputFilename, vector<T> &a) {
  ofstream outputFile(outputFilename);
  if (!outputFile) {
    cout << "Error opening output file: exiting " << endl;
    exit(33);
  }
  int mN = a.size();
  for (unsigned int i = 0; i < mN; i++) {
    outputFile << a[i] << endl;
  }
  outputFile.close();
};

void Print(const char *outputFilename, string text, double a) {
  ofstream outputFile(outputFilename, ios::app);
  outputFile << text << a << endl;
  outputFile.close();
};

template <typename T> void Print(vector<T> &a) {
  unsigned int mN = a.size();
  for (unsigned int k = 0; k < mN; k++) {
    cout << a[k] << endl;
  }
};
