#pragma once
#include <cmath>
#include <iostream>

using namespace std;

class Posizione {

public:
  Posizione();
  Posizione(double x, double y, double z);

  ~Posizione();

  double getX() const { return m_x; };
  double getY() const { return m_y; };
  double getZ() const { return m_z; };

  double getR() const;
  double getPhi() const;
  double getTheta() const;
  double getRho() const;

  Posizione &operator=(const Posizione &p);
  bool operator!=(const Posizione &p) const;
  bool operator==(const Posizione &p) const;
  Posizione operator+(const Posizione &p) const;
  Posizione operator-(const Posizione &p) const;

  Posizione pos(const Posizione &) const;
  double getDist(const Posizione &) const;

private:
  double m_x, m_y, m_z;
};
