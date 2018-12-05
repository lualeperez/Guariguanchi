#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2D.h>
#include <TText.h>
#include <TSystem.h>
#include <TLine.h>
#include <TFile.h>
#include <TEllipse.h>
#include <TVector2.h>
#include <TVector3.h>
#include "include/GBFieldConstant.h"

//C++, C
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
#include <vector>
#include <cstring>

using namespace std;

//====================================================================
GBFieldConstant::GBFieldConstant(TString   aName,
				 TVector3  aBField,
				 GGlobalTools* aglobal)
                                 : GBField(aName,aglobal)
{
  
  GlobalBField = aBField;
  Type         = TString("Constant");
  
}
//====================================================================
GBFieldConstant::GBFieldConstant(TString   aName,
				 TVector3  aBFieldDirection,
				 float     aBFieldMagnitude,
				 GGlobalTools* aglobal)
                                 : GBField(aName,aglobal)
{
  
  GlobalBField = aBFieldMagnitude*(aBFieldDirection.Unit());
  Type         = TString("Constant");
  
}
//====================================================================
GBFieldConstant::GBFieldConstant(TString   aName,
				 float     aBFieldXComponent,
				 float     aBFieldYComponent,
				 float     aBFieldZComponent,
				 GGlobalTools* aglobal)
                                 : GBField(aName,aglobal)
{
  
  GlobalBField = TVector3(aBFieldXComponent,aBFieldYComponent,aBFieldZComponent);
  Type         = TString("Constant");
  
}
//====================================================================
GBFieldConstant::GBFieldConstant(const GBFieldConstant& other,TString aName)
                                 : GBField(aName,other.global)
{
  
  GlobalBField = other.GlobalBField;
  Type         = other.Type;
  
}
//====================================================================
GBFieldConstant::~GBFieldConstant() 
{
  
}
//====================================================================
GBField* GBFieldConstant::clone(TString aName) const
{
 
  return new GBFieldConstant(*this,aName);
  
}
//====================================================================
void  GBFieldConstant::Print()
{
  
  TVector3 unitVect = GlobalBField.Unit();
  double   Mag      = GlobalBField.Mag();
  
  cout << "  Begin Constant B Field" << endl;
  cout << "    Magnitude = "  << Mag/global->GetUnit("T") << " Tesla" << endl;
  cout << "    Direction = (" << unitVect.X() << "," << unitVect.Y() << "," << unitVect.Z() << ")" << endl;
  cout << "  End   Constant B Field" << endl;
  
}
//====================================================================


