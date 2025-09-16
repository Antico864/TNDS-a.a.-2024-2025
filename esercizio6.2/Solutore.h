#pragma once

#include <iostream>

#include "FunzioneBase.h"
#include "sign.h"

using namespace std;

class Solutore {

public:
  ///////////////////
  //  Costruttori
  ///////////////////

  Solutore() { ; };
  Solutore(double prec) { m_prec = prec; };
  Solutore(unsigned int n_max) { m_nmax = n_max; };

  virtual ~Solutore() { ; };

  ///////////////////
  //  CercaZeri:
  ///////////////////

  virtual double CercaZeriRef(double xmin, double xmax, const FunzioneBase &f,
                              double prec = 0.001, unsigned int nmax = 100) = 0;

  ///////////////////
  //  Metodi Get:
  ///////////////////

  unsigned int GetNMaxIt() { return m_nmax; };
  unsigned int GetNIt() { return m_niterations; };
  double GetPrec() { return m_prec; };

  ///////////////////
  //  Metodi Set:
  ///////////////////

  void SetNMaxIt(unsigned int nint) { m_nmax = nint; };
  void SetPrec(double epsilon) { m_prec = epsilon; };

protected:
  double m_a, m_b; //  Estremi intervallo di ricerca
  double m_prec;
  unsigned int m_nmax;        //  Massime
  unsigned int m_niterations; //  Effettuate
};

class Bisezione : public Solutore {

public:
  ///////////////////
  //  Costruttori:
  ///////////////////

  Bisezione() { ; };
  Bisezione(double prec) { m_prec = prec; };
  Bisezione(unsigned int n) { m_nmax = n; };
  virtual ~Bisezione() { ; };

  ///////////////////
  //  CercaZeri:
  ///////////////////

  double CercaZeriRef(double xmin, double xmax, const FunzioneBase &f,
                              double prec = 0.001,
                              unsigned int nmax = 100) override;
};
