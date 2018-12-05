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
#include "include/GEfficiencyModel.h"

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
GEfficiencyModel::GEfficiencyModel(TString   aName,
				   GGlobalTools* aglobal)
{
  
  Name    = aName;
  global  = aglobal;
  Type    = TString("");
  verbose = false;
  
}
//====================================================================
GEfficiencyModel::GEfficiencyModel(const GEfficiencyModel& other,TString aName)
{
  
  Name    = aName;
  global  = other.global;
  Type    = other.Type;
  verbose = other.verbose;
  
}
//====================================================================
GEfficiencyModel::~GEfficiencyModel() 
{
  
}
//====================================================================
GEfficiencyModel* GEfficiencyModel::clone(TString aName) const
{
 
  return new GEfficiencyModel(*this,aName);
  
}
//====================================================================


