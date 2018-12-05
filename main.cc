#include "TROOT.h"
#include "include/Guariguanchi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <assert.h>

using namespace std;

int main(int argc, char *argv[]) {
  
  cout << endl;
  cout << "**************************************************************************" << endl;
  cout << "*                W E L C O M E  to  GUARIGUANCHI                         *" << endl;
  cout << "*     Framework for numerical calculation of tracking performances       *" << endl;
  cout << "*                                                                        *" << endl;
  cout << "* Build info:                                                            *" << endl;
  cout << "*  Version : DEV on trunk                                                *" << endl;
  cout << "*  SVN url : http://sbgtrac.in2p3.fr/svn/CMOS/Guariguanchi/trunk/        *" << endl;
  cout << "*  SVN rev : 2.5M                                                        *" << endl;
  cout << "*  ROOT ver: 5.34/32                                                     *" << endl;
  cout << "*                                                                        *" << endl;
  cout << "* Authors:                                                               *" << endl;
  cout << "*  A. Perez Perez(1)                                                     *" << endl;
  cout << "* Institutions:                                                          *" << endl;
  cout << "*  (1) IPHC - Strasbourg                                                 *" << endl;
  cout << "*                                                                        *" << endl;
  cout << "**************************************************************************" << endl;
  cout << endl;
  
  const char* datacard = "DataCards/Datacard_file.txt";

  if(argc == 2) {
    datacard         = argv[1];
  }
  else {
    cout << endl;
    cout << "Wrong number of input parameters" << endl;
    cout << endl;
    assert(false);
  }

  cout << endl;
  cout << "=============================================" << endl;
  cout << "The function arguments are:"                   << endl;
  cout << "---------------------------------------------" << endl;
  cout << "datacard  = " << datacard                      << endl;
  cout << "=============================================" << endl;
  cout << endl;

  Guariguanchi t(datacard);
  t.DoAnalysis();

  return 0;  

}
