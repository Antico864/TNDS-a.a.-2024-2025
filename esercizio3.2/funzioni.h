#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

template <typename T> vector<T> Read(int n, char *filename) {
  vector<T> v;
  T val;
  ifstream inputFile(filename);
  if (!inputFile) {
    cout << "Cannot open file " << filename << endl;
    exit(11);
  } else {
    for (int k = 0; k < n; k++) {
      inputFile >> val;
      v.push_back(val);
    };
  };
  inputFile.close();
  return v;
};

template <typename T> double CalcAvg(vector<T> a) {
  double sum = 0;
  int nM = a.size();
  for (unsigned int k = 0; k < nM; k++) {
    sum = sum + a[k];
  }
  double avg = sum / nM;
  return avg;
};

template <typename T> double CalcVar(vector<T> a) {
  double var = 0;
  int mN = a.size();
  double avg = CalcAvg(a);
  for (unsigned int k = 0; k < mN; k++) {
    var += pow(a[k] - avg, 2);
  }
  return var / (a.size() - 1);
};

template <typename T> double CalcMed(vector<T> a) {
  double mediana = 0;
  int mN = a.size();
  sort(a.begin(), a.end());
  if (mN % 2 != 0) {
    mediana = a[(mN + 1) / 2];
  } else {
    mediana = (a[mN / 2 - 1] + a[mN / 2]) / 2;
  }
  return mediana;
};

template <typename T> void Print(const char *outputFilename, vector<T> a) {
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

template <typename T> void Print(vector<T> a) {
  unsigned int mN = a.size();
  for (unsigned int k = 0; k < mN; k++) {
    cout << a[k] << endl;
  }
};
