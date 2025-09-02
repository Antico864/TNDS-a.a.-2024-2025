#pragma once

#include "RandomGen.hpp"
#include <iomanip>

using namespace std;

class PrismExperiment {

public:
  PrismExperiment(unsigned int seed);
  ~PrismExperiment() { ; };

  //  Metodi invocati in successione per ogni pseudomisura:
  void Execute();
  void Analyse();

  //  Eventuali metodi Get per accedere ai data members.
  //  Se ci serviranno li scriviamo.

  double Getth0inp() { return m_th0_input; };
  double Getth0outp() { return m_th0_meas; };

  double Getth1inp() { return m_th1_input; };
  double Getth1outp() { return m_th1_meas; };

  double Getth2inp() { return m_th2_input; };
  double Getth2outp() { return m_th2_meas; };

  double GetAinp() { return m_A_input; };
  double GetAoutp() { return m_A_meas; };

  double GetBinp() { return m_B_input; };
  double GetBoutp() { return m_B_meas; };

  double Getdm1inp() { return m_dm1_input; };
  double Getdm1outp() { return m_dm1_meas; };

  double Getdm2inp() { return m_dm2_input; };
  double Getdm2outp() { return m_dm2_meas; };

  double Getn1inp() { return m_n1_input; };
  double Getn1outp() { return m_n1_meas; };

  double Getn2inp() { return m_n2_input; };
  double Getn2outp() { return m_n2_meas; };

private:
  // generatore di numeri casuali

  RandomGen m_rgen;

  // parametri dell'apparato sperimentale:

  double m_lambda1, m_lambda2, m_alpha, m_sigmat;

  // valori delle quantita' misurabili :
  // input  : valori assunti come ipotesi nella simulazione
  // meas : valore dopo la simulazione di misura

  // A e B:
  double m_A_input, m_A_meas;
  double m_B_input, m_B_meas;

  // n1 e n2:
  double m_n1_input, m_n1_meas;
  double m_n2_input, m_n2_meas;

  // deltam1 e deltam2:
  double m_dm1_input, m_dm1_meas;
  double m_dm2_input, m_dm2_meas;

  // theta0, theta1 e theta2:
  double m_th0_input, m_th0_meas;
  double m_th1_input, m_th1_meas;
  double m_th2_input, m_th2_meas;

  //  In breve: utti i valori, "veri" e misurati,
  //  di tutte le grandezze coinvolte nel processo di simulazione.
};
