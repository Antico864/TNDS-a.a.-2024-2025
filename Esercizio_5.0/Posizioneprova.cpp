#include "Posizione.h"
#include <cmath>
#include <iostream>

using namespace std;

int main() {
  double x_1 = 0;
  double x_2 = 0;
  double x_3 = 0;

  cout << "Inserire x " << endl;
  cin >> x_1;

  cout << "Inserire y " << endl;
  cin >> x_2;

  cout << "Inserire z " << endl;
  cin >> x_3;

  Posizione vett(x_1, x_2, x_3);

  cout << "Coordinate cilindriche: " << endl;
  cout << "Raggio planare  #Rho = " << vett.getR() * cos(vett.getRho()) << endl;
  cout << "Angolo azimutale #Theta = " << vett.getTheta() << endl;
  cout << "z = " << vett.getZ() << "\n" << endl;

  cout << "Coordinate sferiche: " << endl;
  cout << "Angolo azimutale #Theta = " << vett.getTheta() << endl;
  cout << "Angolo polare #Rho = " << vett.getRho() << endl;
  cout << "Raggio R = " << vett.getR() << endl;

  return 0;
}
