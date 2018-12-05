/***************************************************************************//**
 * @brief:      
 * @Description: Base class for Resolution models
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

#ifndef GResolutionModel_h
#define GResolutionModel_h

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

class GResolutionModel {
  
private:
  
public:
  
  TString   Name;  // Resolution model name
  TString   Type;  // Resolution model type: TPC, ...
  
  GGlobalTools* global;
  
  // Constructor
  GResolutionModel(TString   aName,
		   GGlobalTools* aglobal);
  
  GResolutionModel(const GResolutionModel& other, TString Name = TString(""));

  // Destructor
  virtual ~GResolutionModel();
  
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
  virtual  GResolutionModel*  clone(TString aName) const;
  virtual  void      CheckInputs() {;}
  virtual  TVector2  GetResolution(TVector2 ResolutionUV, TVector3 mom, TVector3 pos) {return ResolutionUV;}
  virtual  void      Print() {;}

protected:

  bool verbose;
  
};

#endif //~ GResolutionModel_h

