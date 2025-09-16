#include <iostream>

using namespace std;

double *readDataFromFile(const char *ifilename, int N);
double doMedia(double *v, int N);
double doVarianza(double *v, int N);
double doMediana(double *v, int N);
void Merge(double *&v, int left, int right);
void MergeSort(double *&v, int left, int right);
void Print(const char *ofilename, double *v, int n);
void Print(const char *ofilename, string s);
void Print(double *v, int N);
