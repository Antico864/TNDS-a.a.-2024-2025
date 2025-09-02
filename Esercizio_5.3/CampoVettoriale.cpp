#include "CampoVettoriale.h"

using namespace std;

////////////////////////////////////////////////
//
//  Costruttori:
//
////////////////////////////////////////////////

CampoVettoriale::CampoVettoriale(const Posizione &p) : Posizione(p) {
  m_Ex = 0;
  m_Ey = 0;
  m_Ez = 0;
};

CampoVettoriale::CampoVettoriale(const Posizione &p, double Ex, double Ey,
                                 double Ez)
    : Posizione(p) {
  m_Ex = Ex;
  m_Ey = Ey;
  m_Ez = Ez;
};

CampoVettoriale::CampoVettoriale(double x, double y, double z, double Ex,
                                 double Ey, double Ez)
    : Posizione(x, y, z) {
  m_Ex = Ex;
  m_Ey = Ey;
  m_Ez = Ez;
};

////////////////////////////////////////////////
//
//  Operazioni:
//
////////////////////////////////////////////////

//  Questo e il modulo non modificano l'oggetto chiamante, possono essere
//  definiti come funzioni esterne.
CampoVettoriale CampoVettoriale::operator+(const CampoVettoriale &V) const {
  if (V.getX() != getX() || V.getY() != getY() || V.getZ() != getZ()) {
    throw invalid_argument(
        "La somma di campi vettoriali deve avvenire nello stesso punto.");
  }

  CampoVettoriale sum(Posizione(getX(), getY(), getZ()));
  sum.setEx(getEx() + V.getEx());
  sum.setEy(getEy() + V.getEy());
  sum.setEz(getEz() + V.getEz());
  return sum;
};

// Usa overloading di operator+...
CampoVettoriale &CampoVettoriale::operator+=(const CampoVettoriale &V) {
  return (*this) = (*this) + V;
};

double CampoVettoriale::Modulo() {
  double modulo = pow(getEx(), 2) + pow(getEy(), 2) + pow(getEz(), 2);
  return sqrt(modulo);
};
