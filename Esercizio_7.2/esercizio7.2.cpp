#include <iomanip>
#include <iostream>

#include "Integral.h"
#include "funzioni.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 2) {
    cerr << "Use of program: " << argv[0] << " <precisione>" << endl;
    exit(-1);
  };
  double prec = atof(argv[1]);

  xsinx myfunction;

  Trapezoidi trapezio(0, M_PI / 2);

  double I = trapezio.Integrate(prec, myfunction);

  unsigned int nsteps = trapezio.GetN();
  double error = (4 / 3) * fabs(trapezio.GetIn1() - trapezio.GetIn());

  cout << "Number of steps = " << nsteps << "; " << endl;
  cout << "Integral = " << I << "; " << endl;
  cout << "Error = " << error << "; " << endl;

  return 0;
}