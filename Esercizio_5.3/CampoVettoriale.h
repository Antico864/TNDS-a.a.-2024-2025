#pragma once

#include "Posizione.h"

using namespace std;

class CampoVettoriale : public Posizione {
public:
  // Costruttori
  CampoVettoriale(const Posizione &);
  CampoVettoriale(const Posizione &, double Ex, double Ey, double Ez);
  CampoVettoriale(double x, double y, double z, double Ex, double Ey,
                  double Ez);

  // Operazioni
  CampoVettoriale &operator+=(const CampoVettoriale &);
  CampoVettoriale operator+(const CampoVettoriale &) const;

  // Metodi
  double getEx() const { return m_Ex; };
  double getEy() const { return m_Ey; };
  double getEz() const { return m_Ez; };

  void setEx(double Ex) { m_Ex = Ex; };
  void setEy(double Ey) { m_Ey = Ey; };
  void setEz(double Ez) { m_Ez = Ez; };

  double Modulo();

private:
  double m_Ex, m_Ey, m_Ez;
};
