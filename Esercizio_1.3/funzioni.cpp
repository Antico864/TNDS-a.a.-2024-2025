#include "funzioni.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

double *readDataFromFile(const char *filename, int n) {
  ifstream inputFile(filename);
  if (!inputFile) {
    cout << "File name not correct" << endl;
    exit(33);
  }
  double *a = new double[n];
  for (int k = 0; k < n; k++) {
    inputFile >> a[k];
  }
  inputFile.close();
  return a;
}

double doMedia(double *data, int ndata) {
  double media = 0;
  double datasum = 0;
  for (int i = 0; i < ndata; i++) {
    datasum = datasum + data[i];
  }
  media = datasum / ndata;
  return media;
}

double doVarianza(double *a, int n) {
  double var = 0;
  double avg = doMedia(a, n);
  for (unsigned int k = 0; k < n; k++) {
    var += pow(a[k] - avg, 2);
  }
  return var / (n - 1);
}

void Merge(double *&v, int left, int right) {
  int mid = left + (right - left) / 2;
  int n1 = mid - left + 1;
  int n2 = right - mid;

  double *leftArray = new double[n1];
  double *rightArray = new double[n2];

  for (int i = 0; i < n1; i++)
    leftArray[i] = v[left + i];
  for (int j = 0; j < n2; j++)
    rightArray[j] = v[mid + 1 + j];

  int i = 0, j = 0, k = left;

  // Fusione ordinata
  while (i < n1 && j < n2) {
    if (leftArray[i] <= rightArray[j]) {
      v[k] = leftArray[i];
      i++;
    } else {
      v[k] = rightArray[j];
      j++;
    }
    k++;
  }

  // Copia degli elementi rimanenti
  while (i < n1) {
    v[k] = leftArray[i];
    i++;
    k++;
  }

  while (j < n2) {
    v[k] = rightArray[j];
    j++;
    k++;
  }

  delete[] leftArray;
  delete[] rightArray;
}

// Ordinamento crescente
void MergeSort(double *&v, int left, int right) {
  if (left < right) {
    int mid = left + (right - left) / 2;

    MergeSort(v, left, mid);
    MergeSort(v, mid + 1, right);

    Merge(v, left, right);
  }
}

double doMediana(double *v, int n) {
  MergeSort(v, 0, n - 1);
  double mediana = 0;
  if (n % 2 == 0) {
    mediana = abs((v[n / 2 + 1] - v[n / 2]) / 2.);
  } else {
    mediana = v[(n) / 2];
  }
  return mediana;
}

void Print(const char *outputFilename, double *v, int n) {
  ofstream outputFile(outputFilename);
  if (!outputFile) {
    cout << "Error opening output file: exiting " << endl;
    exit(33);
  }
  for (unsigned int i = 0; i < n; i++) {
    outputFile << v[i] << endl;
  }
  outputFile.close();
}

void Print(const char *outputFilename, string text) {
  ofstream outputFile(outputFilename, ios::app);
  outputFile << text << endl;
  outputFile.close();
}

void Print(double *v, int n) {
  for (int k = 0; k < n; k++) {
    cout << v[k] << endl;
  }
}
