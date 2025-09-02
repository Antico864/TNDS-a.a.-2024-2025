#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

template <typename T> vector<T> Read(const char *filename) {
  vector<T> v;
  T val;
  ifstream inputFile(filename);
  if (!inputFile) {
    cout << "Cannot open file " << filename << endl;
    exit(11);
  } else {
    while (!inputFile.eof()) {
      inputFile >> val;
      v.push_back(val);
    };
  };
  inputFile.close();
  return v;
};

template <typename T> double DoAvg(const vector<T>& a) {
  if (a.empty()) return 0.;
  double sum = 0;
  for (unsigned int i = 0; i < a.size(); i++) {
    sum += a[i];
  }
  return sum / a.size();
}

template <typename T> double DoErr(vector<T> a, double avg) {
  double delta2sum = 0;
  unsigned int counter = 0;
  for (unsigned int k = 0; k < a.size(); k = k + 7) {
    delta2sum += pow(a[k] - avg, 2);
    counter++;
  }
  double devstdavg = sqrt(delta2sum / counter) / sqrt(counter);
  return devstdavg;
};
