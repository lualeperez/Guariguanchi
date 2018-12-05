/***************************************************************************//**
 * @brief:      
 * @Description: Base class for B-fields
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
 

#ifndef GBField_h
#define GBField_h

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

class GBField {
  
private:
  
public:
  
  TString   Name;                // B-field name
  TString   Type;                // B-field type: constant, ...
  
  GGlobalTools* global;
  
  // Constructor
  GBField(TString   aName,
	  GGlobalTools* aglobal);
  
  GBField(const GBField& other, TString Name = TString(""));

  // Destructor
  virtual ~GBField();
  
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
  virtual  GBField*  clone(TString aName) const;
  virtual  TVector3  GetBFieldValue(TVector3 PositionXYZ) {return TVector3(0.0,0.0,0.0);}
  virtual  void      Print() {;}

protected:

  bool verbose;
  
};

#endif //~ GBField_h

