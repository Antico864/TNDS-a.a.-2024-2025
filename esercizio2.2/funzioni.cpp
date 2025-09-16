#include "funzioni.h"
#include "Vettore.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

Vettore Read(int n, const char *filename) {
  ifstream inputFile(filename);
  if (!inputFile) {
    cout << "File name not correct" << endl;
    exit(33);
  } else {
    Vettore vett0(n);
    for (unsigned int k = 0; k < n; k++) {
      inputFile >> vett0[k];
    }
    return vett0;
  }
}

double CalcAvg(Vettore &a) {
  double sum = 0;
  int nM = a.GetN();
  for (unsigned int k = 0; k < nM; k++) {
    sum = sum + a.GetComponent(k);
  }
  double avg = sum / nM;
  return avg;
}

double CalcVar(Vettore &a) {
  double var = 0;
  int mN = a.GetN();
  double avg = CalcAvg(a);
  for (unsigned int k = 0; k < mN; k++) {
    var += pow(a[k] - avg, 2);
  }
  return var / (a.GetN() - 1);
}

void Merge(Vettore &v, int left, int right) {
  int mid = left + (right - left) / 2;
  int n1 = mid - left + 1;
  int n2 = right - mid;

  Vettore leftVect(n1);
  Vettore rightVect(n2);

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
}

void MergeSort(Vettore &v, int left, int right) {
  if (left < right) {
    int mid = left + (right - left) / 2;

    MergeSort(v, left, mid);
    MergeSort(v, mid + 1, right);

    Merge(v, left, right);
  }
}

double CalcMed(Vettore &a) {
  MergeSort(a, 0, a.GetN() - 1);
  double mediana = 0;
  int mN = a.GetN();
  if (mN % 2 != 0) {
    mediana = a[(mN) / 2];
  } else {
    mediana = abs((a[mN / 2 + 1] - a[mN / 2]) / 2);
  }
  return mediana;
}

void Print(const char *outputFilename, Vettore &a) {
  ofstream outputFile(outputFilename);
  if (!outputFile) {
    cout << "Error opening output file: exiting " << endl;
    exit(33);
  }
  int mN = a.GetN();
  for (unsigned int i = 0; i < mN; i++) {
    outputFile << a[i] << endl;
  }
  outputFile.close();
}

void Print(const char *outputFilename, string text) {
  ofstream outputFile(outputFilename, ios::app);
  outputFile << text << endl;
  outputFile.close();
}

void Print(Vettore &a) {
  unsigned int mN = a.GetN();
  for (unsigned int k = 0; k < mN; k++) {
    cout << a[k] << endl;
  }
}
