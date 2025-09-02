#pragma once
#include <cmath>
#include <iostream>

using namespace std;

class Posizione {

public:
  Posizione();
  Posizione(double x, double y, double z);

  // Distruttore. Se lo scrivo lo devo invocare.
  ~Posizione();

  double getX() const { return m_x; }; //
  double getY() const { return m_y; }; //
  double getZ() const { return m_z; }; //

  double getR() const;     //
  double getPhi() const;   //
  double getTheta() const; //
  double getRho() const;   //

  // Riporta la posizione di un oggetto posizione in coordinate cartesiane
  // rispetto alla posizione inserita:
  Posizione pos(const Posizione &) const;

  double getDist(const Posizione &) const;

private:
  double m_x, m_y, m_z;
};
