/***************************************************************************//**
 * @brief:      
 * @Description: Base class for Detection Efficiency models
 *
 *
 * @createdby:  PEREZ PEREZ Alejandro <luis_alejandro.perez_perez@iphc.cnrs.fr> at 2017-10-01 14:07:38
 * @copyright:  (c)2017 IPHC - CNRS - Université de Strasbourg. All Rights Reserved.
 * 
 * @License: You are free to use this source files for your own development as long
 *           as it stays in a public research context. You are not allowed to use it
 *           for commercial purpose. You must put this header with laboratory and
 *           authors names in all development based on this library.
 *           When results obtained with this package (Guariguanchi) are communicated, 
 *           you should quote the package name.
 *
 * @lastchange: $Revision$
 *              $Author$
 *              $Date$
 *
 *******************************************************************************/

#ifndef GEfficiencyModel_h
#define GEfficiencyModel_h

#include "TMath.h"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TCutG.h>
#include "TStopwatch.h"
#include "include/GGlobalTools.h"

#include <algorithm>
#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

class GEfficiencyModel {
  
private:
  
public:
  
  TString   Name;  // Efficiency model name
  TString   Type;  // Efficiency model type: not other type defined thus far ...
  
  GGlobalTools* global;
  
  // Constructor
  GEfficiencyModel(TString   aName,
		   GGlobalTools* aglobal);
  
  GEfficiencyModel(const GEfficiencyModel& other, TString Name = TString(""));

  // Destructor
  virtual ~GEfficiencyModel();
  
  //Set of generic functions
  
  //Functions to access internal variables
  TString   GetName() const { return  Name; }
  TString   GetType() const { return  Type; }
  bool      GetVerbose() { return verbose; }
  
  //Functions to set internal variables
  void      SetName(TString aName)  { Name = aName; }
  void      SetType(TString aType)  { Type = aType; }
  void      SetVerbose(bool aVerbose) { verbose = aVerbose; }
  
  //Set of functions to be defined for the daughter classes
  virtual  GEfficiencyModel*  clone(TString aName) const;
  virtual  void      CheckInputs() {;}
  virtual  double    GetEfficiency(double DetEffic, TVector3 mom, TVector3 pos, TString aParticle) {return DetEffic;}
  virtual  void      Print() {;}

protected:

  bool verbose;
  
};

#endif //~ GEfficiencyModel_h

