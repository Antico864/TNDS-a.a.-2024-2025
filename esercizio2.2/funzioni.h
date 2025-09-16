#pragma once
#include "Vettore.h"
#include <iostream>

using namespace std;

Vettore Read(int N, const char *ifilename);
double CalcAvg(Vettore &v);
double CalcVar(Vettore &v);
double CalcMed(Vettore &v);
void Merge(Vettore &v, int left, int right);
void MergeSort(Vettore &v, int left, int right);
void Print(const char *ofilename, Vettore &v);
void Print(const char *ofilename, string s);
void Print(Vettore &v);
