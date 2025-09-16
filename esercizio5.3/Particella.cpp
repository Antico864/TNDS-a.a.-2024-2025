#include "Particella.h"

using namespace std;

////////////////////////////
//
//  Classe madre Particella:
//
////////////////////////////

Particella::Particella() {
  m_massa = 0;
  m_carica = 0;
};

Particella::Particella(double m, double q) {
  m_massa = m;
  m_carica = q;
};

Particella::~Particella() {;};

void Particella::Print() const {
  cout << "Particella: m = " << getMassa() << "    q = " << getCarica() << endl;
};

////////////////////////////
//
//  Elettrone:
//
////////////////////////////

Elettrone::Elettrone() : Particella(9.1093826E-31, -1.60217653E-19) {};

Elettrone::~Elettrone() {;};

void Elettrone::Print() const {
  cout << "Elettrone: m = " << m_massa << "    q = " << m_carica << endl;
};

Protone::Protone() : Particella(9.1093826E-31, 1.60217653E-19) {};

Protone::~Protone() {};

void Protone::Print() const {
  cout << "Protone: m = " << m_massa << "\tq = " << m_carica << endl;
};
