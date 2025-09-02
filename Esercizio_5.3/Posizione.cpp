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
};

Posizione::Posizione(double x, double y, double z) {
  m_x = x;
  m_y = y;
  m_z = z;
};

Posizione::~Posizione() { ; }

// Raggio rispetto all'origine degli assi:
double Posizione::getR() const {
  return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));
};

/////////////////////////////////
//
//  Metodi get:
//
/////////////////////////////////

//  Angolo con l'asse y:
double Posizione::getPhi() const { return atan2(m_x, m_y); };

//  Angolo con l'asse x:
double Posizione::getTheta() const { return atan2(m_y, m_x); };

//  Angolo con l'asse z:
double Posizione::getRho() const {
  Posizione p(m_x, m_y, m_z);
  double r = p.getR();
  if (r == 0) {
    cout << "Errore: divisione per zero nel calcolo di rho." << endl;
    exit(-1);
  }
  double cosRho = (m_z - p.getZ()) / r;
  if (cosRho > 1)
    cosRho = 1; // Correzione numerica
  if (cosRho < -1)
    cosRho = -1;
  return acos(cosRho);
};

/////////////////////////////////
//
//  Operatori:
//
/////////////////////////////////

Posizione &Posizione::operator=(const Posizione &p) {
  if (this != &p) { // Evita l'auto-assegnazione
    this->m_x = p.getX();
    this->m_y = p.getY();
    this->m_z = p.getZ();
  }
  return *this;
};

bool Posizione::operator!=(const Posizione &p) const {
  return m_x != p.getX() || m_y != p.getY() || m_z != p.getZ();
};

bool Posizione::operator==(const Posizione &p) const {
  return m_x == p.getX() && m_y == p.getY() && m_z == p.getZ();
};

Posizione Posizione::operator+(const Posizione &p) const {
  return Posizione(m_x + p.getX(), m_y + p.getY(), m_z + p.getZ());
};

Posizione Posizione::operator-(const Posizione &p) const {
  return Posizione(m_x - p.getX(), m_y - p.getY(), m_z - p.getZ());
};

/////////////////////////////////
//
//  Modifica di sdrif:
//
/////////////////////////////////

// Distanza da un punto:
double Posizione::getDist(const Posizione &a) const {
  Posizione diff = a - *this;
  return diff.getR();
};

// Posizione rispetto a un punto:
Posizione Posizione::pos(const Posizione &a) const {
  Posizione diff(getX() - a.getX(), getY() - a.getY(), getZ() - a.getZ());
  return diff;
};
