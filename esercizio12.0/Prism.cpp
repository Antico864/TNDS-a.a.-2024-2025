#include "Prism.hpp"

using namespace std;

PrismExperiment::PrismExperiment(unsigned int seed)
    : m_rgen(seed), m_lambda1(579.1E-9), m_lambda2(404.7E-9),
      m_alpha(60. * M_PI / 180.), m_sigmat(0.3E-3), m_A_input(2.7),
      m_B_input(60000E-18) {
  // Indici di rifrazione;
  m_n1_input = sqrt(m_A_input + m_B_input / (m_lambda1 * m_lambda1));
  m_n2_input = sqrt(m_A_input + m_B_input / (m_lambda2 * m_lambda2));

  // Imposto theta0 (completamente arbitrario);
  m_th0_input = M_PI / 2;

  // Calcolo le deviazioni minime
  // e theta1, theta2:
  m_dm1_input = 2. * asin(m_n1_input * sin(0.5 * m_alpha)) - m_alpha;
  m_th1_input = m_th0_input + m_dm1_input;
  m_dm2_input = 2. * asin(m_n2_input * sin(0.5 * m_alpha)) - m_alpha;
  m_th2_input = m_th0_input + m_dm2_input;
}

void PrismExperiment::Execute() {
  // Ottengo le pseudomisure degli angoli.
  // Ripercorro i passi che farebbe lo sperimentatore,
  // distribuendo le costanti in modo gaussiano.

  //  "Misuro" gli angoli perturbando i dati iniziali:
  m_th0_meas = m_rgen.GaussBM(m_th0_input, m_sigmat);
  m_th1_meas = m_rgen.GaussBM(m_th1_input, m_sigmat);
  m_th2_meas = m_rgen.GaussBM(m_th2_input, m_sigmat);
}

void PrismExperiment::Analyse() {
  // Determino le deviazioni minime:
  m_dm1_meas = m_th1_meas - m_th0_meas;
  m_dm2_meas = m_th2_meas - m_th0_meas;

  // Inverto per trovare gli indici di rifrazione:
  m_n1_meas = sin((m_dm1_meas + m_alpha) / 2) / sin(m_alpha / 2);
  m_n2_meas = sin((m_dm2_meas + m_alpha) / 2) / sin(m_alpha / 2);

  // Inverto la relazione di Cauchy per trovare A e B:
  m_A_meas = (pow(m_lambda2 * m_n2_meas, 2) - pow(m_lambda1 * m_n1_meas, 2)) /
             (m_lambda2 * m_lambda2 - m_lambda1 * m_lambda1);
  m_B_meas = (pow(m_n2_meas, 2) - pow(m_n1_meas, 2)) /
             (1 / (m_lambda2 * m_lambda2) - 1 / (m_lambda1 * m_lambda1));
}
