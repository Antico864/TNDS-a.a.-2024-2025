#include "PuntoMateriale.h"

using namespace std;

/////////////////////////////////
//
//  Costruttori:
//
/////////////////////////////////

PuntoMateriale::PuntoMateriale(double m, double q, const Posizione &p)
    : Posizione(p), Particella(m, q) {;};

PuntoMateriale::PuntoMateriale(double m, double q, double x, double y, double z)
    : Posizione(x, y, z), Particella(m, q) {;};

PuntoMateriale::PuntoMateriale(const Particella &part, const Posizione &p)
    : Posizione(p), Particella(part) {;};

/////////////////////////////////
//
//  Fisica:
//
/////////////////////////////////

CampoVettoriale PuntoMateriale::CampoElettrico(const Posizione &p) const {

  CampoVettoriale E(p);

  const double E_0 =
      getCarica() / (4 * (M_PI) * 8.8541878188e-12 * pow(getDist(p), 2));

  E.setEx(E_0 * pos(p).getX() / getDist(p));
  E.setEy(E_0 * pos(p).getY() / getDist(p));
  E.setEz(E_0 * pos(p).getZ() / getDist(p));

  return E;
};

CampoVettoriale PuntoMateriale::CampoGravitazionale(const Posizione &p) const {

  CampoVettoriale G(p);

  const double G_0 = 6.67430e-11 * getMassa() / pow(getDist(p), 2);

  G.setEx(G_0 * pos(p).getX());
  G.setEy(G_0 * pos(p).getY());
  G.setEz(G_0 * pos(p).getZ());

  return G;
};
