#include <cmath>
#include <iostream>
#include <tgmath.h>

#include "Posizione.h"

using namespace std;

/////////////////////////////////
//
//  Inizializzazione struttura:
//
/////////////////////////////////

Posizione::Posizione() {
  m_x = 0;
  m_y = 0;
  m_z = 0;
}

Posizione::Posizione(double x, double y, double z) {
  m_x = x;
  m_y = y;
  m_z = z;
}

Posizione::~Posizione() { ; }

double Posizione::getR() const {
  double R = 0;
  R = sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));
  return R;
}

/////////////////////////////////
//
//  Metodi get:
//
/////////////////////////////////

//  Angolo con l'asse y:
double Posizione::getPhi() const {
  double Phi = 0;
  Phi = atan(m_x / m_y);
  return Phi;
}

//  Angolo con l'asse x:
double Posizione::getTheta() const {
  double Theta = 0;
  Theta = atan(m_y / m_x);
  return Theta;
}

//  Angolo con l'asse z:
double Posizione::getRho() const {
  double Rho = 0;
  Rho = acos(m_z / getR());
  return Rho;
}

double Posizione::getDist(const Posizione &a) const {
  Posizione diff = pos(a);
  double dist = diff.getR();
  return dist;
};

Posizione Posizione::pos(const Posizione &a) const {
  Posizione diff(getX() - a.getX(), getY() - a.getY(), getZ() - a.getZ());
  return diff;
};
