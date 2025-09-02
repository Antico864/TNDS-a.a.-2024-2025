#pragma once

#include "CampoVettoriale.h"
#include <cmath>
#include <iostream>

using namespace std;

class Particella {
public:
  // Costrutori
  Particella();
  Particella(double massa, double carica);

  // Distruttore
  virtual ~Particella();

  // Metodi:
  double getMassa() const { return m_massa; };
  double getCarica() const { return m_carica; };
  virtual void Print() const;

protected:
  double m_massa;
  double m_carica;
};

class Elettrone : public Particella {
public:
  // costruttore
  Elettrone();
  // distruttore
  virtual ~Elettrone();
  //
  virtual void Print() const;
};

class Protone : public Particella {
public:
  // costruttore
  Protone();
  // distruttore
  virtual ~Protone();
  //
  virtual void Print() const;
};
