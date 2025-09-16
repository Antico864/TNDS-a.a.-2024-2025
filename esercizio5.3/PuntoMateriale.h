#pragma once

#include "CampoVettoriale.h"
#include "Particella.h"
#include "Posizione.h"

using namespace std;

class PuntoMateriale : public Posizione, public Particella {
public:
  // Costruttori:
  PuntoMateriale(double m, double q, const Posizione &p);
  PuntoMateriale(double m, double q, double x, double y, double z);
  PuntoMateriale(const Particella &part, const Posizione &p);

  // Fisica:
  CampoVettoriale CampoElettrico(const Posizione &p) const;
  CampoVettoriale CampoGravitazionale(const Posizione &p) const;
};
