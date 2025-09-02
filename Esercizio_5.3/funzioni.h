#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
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
}

template <typename T> double DoAvg(vector<T> a) {
  double avg = 0;
  for (int k = 0; k < a.size(); k++) {
    avg = (avg * (k) + a[k + 1]) / (k + 1);
  }
  return avg;
}

template <typename T> double DoErr(vector<T> a, double avg) {
  double delta2 = 0;
  double delta2sum = 0;
  unsigned int counter = 0;
  for (unsigned int k = 0; k < a.size(); k = k + 7) {
    delta2 = pow(a[k] - avg, 2);
    delta2sum = delta2sum + delta2;
    counter++;
  }
  double devstdavg = sqrt(delta2sum / (counter - 1)) / sqrt(counter);
  return devstdavg;
}
