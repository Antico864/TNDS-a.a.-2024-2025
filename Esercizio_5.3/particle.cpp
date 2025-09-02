#include <iostream>

#include "Particella.h"

using namespace std;

int main() {

  Particella p(1., -1.60217653e-19);
  Elettrone *e = new Elettrone;

  cout << "Particella p: m = " << p.getMassa() << "   q = " << p.getCarica()
       << endl;
  p.Print();

  cout << "Elettrone: m = " << e->getMassa() << "   q = " << e->getCarica()
       << endl;
  e->Print();

  return 0;
}
