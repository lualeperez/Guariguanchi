///////////////////////////////////////////////////////////////////////////
// This file handles the function for the readout of the input datacards //
///////////////////////////////////////////////////////////////////////////

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH2D.h>
#include <TText.h>
#include <TSystem.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TVector2.h>
#include <TVector3.h>
#include "include/Guariguanchi.h"

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
void  Guariguanchi::ReadDataCard(void)
{

  //Readout of the input datacard

  // create a file-reading object
  ifstream fin;
  fin.open(DataCard.Data()); // open a file
  if(!fin.good()) {
    cout << endl;
    cout << "File " << DataCard.Data() << " not found. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }
  
  int line_number = 0;
  //read each line of the file
  while(!fin.eof()) {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
    line_number++;

    // parse the line into blank-delimited tokens
    int n = 0; // a for-loop index
    
    // array to store memory addresses of the tokens in buf
    const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
    
    // parse the line
    token[0] = strtok(buf, DELIMITER); // first token
    
    // go to next line if the line is blank of it starts with comment "//"
    if(!token[0] || TString(token[0]).BeginsWith(COMMENT)) continue;
    
    for(n=1; n<MAX_TOKENS_PER_LINE; n++) {
      token[n] = strtok(0, DELIMITER); // subsequent tokens
      
      if(!token[n] || TString(token[n]).BeginsWith(COMMENT)) break; // no more tokens
    }
    
    if(n == 2 && TString(token[0]) == TString("ParticleType:")) particle = TString(token[1]);
    else if(n == 5 && TString(token[0]) == TString("ParticleOrigin:")) {
      ParticleOrigin  = TVector3(atof(token[1]),
	                         atof(token[2]),
				 atof(token[3]));
      if(!global->IsDistanceUnit(token[4])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      ParticleOrigin *= global->GetDistanceUnit(TString(token[4]));
    }
    else if(n == 5 && TString(token[0]) == TString("ReferencePoint:")) {
      ReferencePoint  = TVector3(atof(token[1]),
	                         atof(token[2]),
				 atof(token[3]));
      if(!global->IsDistanceUnit(token[4])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      ReferencePoint *= global->GetDistanceUnit(TString(token[4]));
    }
    else if(n == 2 && TString(token[0]) == TString("Verbose:")) {
      verbose = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("MCSeed:")) {
      MCSeed = atoi(token[1]);
    }
    else if(n == 3 && TString(token[0]) == TString("KappaElossFluctuation:")) {      
      if(!global->IsEnergyUnit(token[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << ", for KappaElossFluctuation." << endl;
      kappa  = atof(token[1])*global->GetEnergyUnit(TString(token[2]));
    }
    else if(n == 2 && TString(token[0]) == TString("FitPowerForImpactParam:")) {
      FitPowerForImpactParam = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("MonResolRepresentation:")) {
      if(TString(token[1]) == TString("sigma(1/Pt)") || TString(token[1]) == TString("sigma(Pt)/Pt") || TString(token[1]) == TString("sigma(Pt)")) MonResolRepresentation = TString(token[1]);
      else {
	cout << endl;
	cout << "WARNING:: the \"MonResolRepresentation\" parameter can only have the values: sigma(Pt)/Pt, sigma(1/Pt) and sigma(Pt) . Setting it to default value sigma(Pt)/Pt." << endl;
	cout << endl;
      }
    }
    else if(n == 2 && TString(token[0]) == TString("PrintGeometry:")) {
      DoPrintGeometry = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("PrintGeometryWeight:")) {
      DoPrintGeometryWeight = global->SetBoolFromString(TString(token[1]));
    }    
    else if(n == 2 && TString(token[0]) == TString("DoGeometryCheck:")) {
      DoGeometryCheck = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 3 && TString(token[0]) == TString("GeoCheckPrecision:")) {
      if(!global->IsDistanceUnit(token[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;      
      GeoCheckPrecision = atof(token[1])*global->GetDistanceUnit(TString(token[2]));
    }
    else if(n == 2 && TString(token[0]) == TString("SavePlots:")) {
      SavePlots = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("IncludeEloss:")) {
      IncludeEloss = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("PlotGeometry:")) {
      DoPlotGeometry = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("DoRZGeoRepresentation:")) {
      DoRZGeoRepresentation = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("PlotWorldVolume:")) {
      DoPlotWorldVolume = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("PlotStepBfieldVolume:")) {
      DoPlotStepBfieldVolume = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("PlotSomeTracks:")) {
      DoPlotSomeTracks = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("UseAllMomVals:")) {
      UseAllMomVals_GeoPrint = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("DoMaterialBugdetAnalysis:")) {
      DoMaterialBugdetAnalysis = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("DoTrkResolAnalysis:")) {
      DoTrkResolAnalysis = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("DoTelescopeAnalysis:")) {
      DoTelescopeAnalysis = global->SetBoolFromString(TString(token[1]));
    }
    else if(n == 2 && TString(token[0]) == TString("DoPseudoEfficVsMon:")) {
      DoPseudoEfficVsMon = global->SetBoolFromString(TString(token[1]));
    }
    else if(TString(token[0]) == TString("MomentumValues:") && momArr.size() == 0) {
      if(n <= 2) {
	cout << endl;
	cout << "The momentum values (\"MomentumValues:\") have to have at least 3 tokens, in the following way:" << endl;
	cout << " MomentumValues:   val1  val2 ... Units" << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      if(!global->IsMomentumUnit(token[n-1])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      double UnitValue = global->GetMomentumUnit(TString(token[n-1]));
      momArr.clear();
      for(int kt=1;kt<n-1;kt++) momArr.push_back(atof(token[kt])*UnitValue);
      momVariable = TString("p");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginMomentumScan") && momArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndMomentumScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < 0)     continue;
	  if(R[1] < 0)     continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  momArr.clear();
	  for(int i=0;i<Nbins+1;i++) momArr.push_back(R[0] + i*(R[1] - R[0])/Nbins);
	  momVariable = TString("p");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins")) Nbins = atoi(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("pMin")) {
	  if(!global->IsMomentumUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[0]  = atof(token1[1])*global->GetMomentumUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("pMax")) {
	  if(!global->IsMomentumUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[1]  = atof(token1[1])*global->GetMomentumUnit(TString(token1[2]));
	}
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of momentum scan parameters inside the BeginMomentumScan and EndMomentumScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginMomentumScan" << endl;
	cout << "  Nbins   N           // N > 0"    << endl;
	cout << "  pMin    val  units  // val >= 0" << endl;
	cout << "  pMax    val  units  // val >= 0" << endl;
	cout << "EndMomentumScan"   << endl;
	cout << "pMin has to be smaller than pMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("EnergyValues:") && momArr.size() == 0) {
      if(n <= 2) {
	cout << endl;
	cout << "The energy values (\"EnergyValues:\") have to have at least 3 tokens, in the following way:" << endl;
	cout << " EnergyValues:   val1  val2 ... Units" << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      if(!global->IsEnergyUnit(token[n-1])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      double UnitValue = global->GetEnergyUnit(TString(token[n-1]));
      momArr.clear();
      for(int kt=1;kt<n-1;kt++) momArr.push_back(atof(token[kt])*UnitValue);
      momVariable = TString("E");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginEnergyScan") && momArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndEnergyScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < 0)     continue;
	  if(R[1] < 0)     continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  momArr.clear();
	  for(int i=0;i<Nbins+1;i++) momArr.push_back(R[0] + i*(R[1] - R[0])/Nbins);
	  momVariable = TString("E");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins")) Nbins = atoi(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("EMin")) {
	  if(!global->IsEnergyUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[0]  = atof(token1[1])*global->GetEnergyUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("EMax")) {
	  if(!global->IsEnergyUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[1]  = atof(token1[1])*global->GetEnergyUnit(TString(token1[2]));
	}
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of momentum scan parameters inside the BeginMomentumScan and EndMomentumScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginEnergyScan" << endl;
	cout << "  Nbins   N           // N > 0"    << endl;
	cout << "  EMin    val  units  // val >= 0" << endl;
	cout << "  EMax    val  units  // val >= 0" << endl;
	cout << "EndEnergyScan"   << endl;
	cout << "pMin has to be smaller than pMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("EkinValues:") && momArr.size() == 0) {
      if(n <= 2) {
	cout << endl;
	cout << "The kinetic energy values (\"EkinValues:\") have to have at least 3 tokens, in the following way:" << endl;
	cout << " EkinValues:   val1  val2 ... Units" << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      if(!global->IsEnergyUnit(token[n-1])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      double UnitValue = global->GetEnergyUnit(TString(token[n-1]));
      momArr.clear();
      for(int kt=1;kt<n-1;kt++) momArr.push_back(atof(token[kt])*UnitValue);
      momVariable = TString("Ekin");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginEkinScan") && momArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndEkinScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < 0)     continue;
	  if(R[1] < 0)     continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  momArr.clear();
	  for(int i=0;i<Nbins+1;i++) momArr.push_back(R[0] + i*(R[1] - R[0])/Nbins);
	  momVariable = TString("Ekin");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins")) Nbins = atoi(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("EkinMin")) {
	  if(!global->IsEnergyUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[0]  = atof(token1[1])*global->GetEnergyUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("EkinMax")) {
	  if(!global->IsEnergyUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[1]  = atof(token1[1])*global->GetEnergyUnit(TString(token1[2]));
	}
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of momentum scan parameters inside the BeginMomentumScan and EndMomentumScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginEkinScan" << endl;
	cout << "  Nbins   N           // N > 0"    << endl;
	cout << "  EkinMin    val  units  // val >= 0" << endl;
	cout << "  EkinMax    val  units  // val >= 0" << endl;
	cout << "EndEkinScan"   << endl;
	cout << "pMin has to be smaller than pMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("EkinPerUValues:") && momArr.size() == 0) {
      if(n <= 2) {
	cout << endl;
	cout << "The kinetic energy per nucleon values (\"EkinPerUValues:\") have to have at least 3 tokens, in the following way:" << endl;
	cout << " EkinPerUValues:   val1  val2 ... Units" << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      if(!global->IsEnergyUnit(token[n-1])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      double UnitValue = global->GetEnergyUnit(TString(token[n-1]));
      momArr.clear();
      for(int kt=1;kt<n-1;kt++) momArr.push_back(atof(token[kt])*UnitValue);
      momVariable = TString("EkinPerU");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginEkinPerUScan") && momArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndEkinPerUScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < 0)     continue;
	  if(R[1] < 0)     continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  momArr.clear();
	  for(int i=0;i<Nbins+1;i++) momArr.push_back(R[0] + i*(R[1] - R[0])/Nbins);
	  momVariable = TString("EkinPerU");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins")) Nbins = atoi(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("EkinPerUMin")) {
	  if(!global->IsEnergyUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[0]  = atof(token1[1])*global->GetEnergyUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("EkinPerUMax")) {
	  if(!global->IsEnergyUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[1]  = atof(token1[1])*global->GetEnergyUnit(TString(token1[2]));
	}
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of momentum scan parameters inside the BeginMomentumScan and EndMomentumScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginEkinPerUScan" << endl;
	cout << "  Nbins   N           // N > 0"    << endl;
	cout << "  EkinPerUMin    val  units  // val >= 0" << endl;
	cout << "  EkinPerUMax    val  units  // val >= 0" << endl;
	cout << "EndEkinPerUScan"   << endl;
	cout << "pMin has to be smaller than pMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("PolarAngleValues:") && thetaArr.size() == 0) {
      if(n <= 2) {
	cout << endl;
	cout << "The polar angle values (\"PolarAngleValues:\") have to have at least 3 tokens, in the following way:" << endl;
	cout << " PolarAngleValues:   val1  val2 ... Units" << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      if(!global->IsAngleUnit(token[n-1])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      double UnitValue = global->GetAngleUnit(TString(token[n-1]));
      thetaArr.clear();
      for(int kt=1;kt<n-1;kt++) thetaArr.push_back(atof(token[kt])*UnitValue);
      polarVariable = TString("theta");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginPolarAngleScan") && thetaArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndPolarAngleScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < 0)     continue;
	  if(R[1] < 0)     continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  thetaArr.clear();
	  for(int i=0;i<Nbins+1;i++) thetaArr.push_back(R[0] + i*(R[1] - R[0])/Nbins);
	  polarVariable = TString("theta");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins"))    Nbins = atoi(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("thetaMin")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[0]  = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("thetaMax")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[1]  = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of polar angle scan parameters inside the BeginPolarAngleScan and EndPolarAngleScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginPolarAngleScan" << endl;
	cout << "  Nbins     N           // N > 0"    << endl;
	cout << "  thetaMin  val  units  // val >= 0" << endl;
	cout << "  thetaMax  val  units  // val >= 0" << endl;
	cout << "EndPolarAngleScan"   << endl;
	cout << "thetaMin has to be smaller than thetaMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("CosThetaValues:") && thetaArr.size() == 0) {
      if(n <= 1) {
	cout << endl;
	cout << "The cos(theta) angle values (\"CosThetaValues:\") have to have at least 2 tokens, in the following way:" << endl;
	cout << " CosThetaValues:   val1  val2 ..." << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      thetaArr.clear();
      for(int kt=1;kt<n;kt++) thetaArr.push_back(global->FromCosThetaToTheta(atof(token[kt])));
      polarVariable = TString("costheta");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginCosThetaScan") && thetaArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndCosThetaScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < -1 || R[0] > 1) continue;
	  if(R[1] < -1 || R[1] > 1) continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  thetaArr.clear();
	  for(int i=0;i<Nbins+1;i++) thetaArr.push_back(global->FromCosThetaToTheta(R[0] + i*(R[1] - R[0])/Nbins));
	  polarVariable = TString("costheta");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins"))       Nbins = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("costhetaMin")) R[0]  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("costhetaMax")) R[1]  = atof(token1[1]);
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of costheta scan parameters inside the BeginCosThetaScan and EndCosThetaScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginCosThetaScan" << endl;
	cout << "  Nbins     N           // N > 0"    << endl;
	cout << "  costhetaMin  val  // val inside (-1,1)" << endl;
	cout << "  costhetaMax  val  // val inside (-1,1)" << endl;
	cout << "EndCosThetaScan"   << endl;
	cout << "costhetaMin has to be smaller than costhetaMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("EtaValues:") && thetaArr.size() == 0) {
      if(n <= 1) {
	cout << endl;
	cout << "The cos(theta) angle values (\"EtaValues:\") have to have at least 2 tokens, in the following way:" << endl;
	cout << " EtaValues:   val1  val2 ..." << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      thetaArr.clear();
      for(int kt=1;kt<n;kt++) thetaArr.push_back(global->FromEtaToTheta(atof(token[kt])));
      polarVariable = TString("eta");
    }
    else if(n == 1 && TString(token[0]) == TString("BeginEtaScan") && thetaArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndEtaScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  thetaArr.clear();
	  for(int i=0;i<Nbins+1;i++) thetaArr.push_back(global->FromEtaToTheta(R[0] + i*(R[1] - R[0])/Nbins));
	  polarVariable = TString("eta");
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins"))  Nbins = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("etaMin")) R[0]  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("etaMax")) R[1]  = atof(token1[1]);
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of eta scan parameters inside the BeginEtaScan and EndEtaScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginEtaScan" << endl;
	cout << "  Nbins     N           // N > 0"    << endl;
	cout << "  etaMin  val  " << endl;
	cout << "  etaMax  val  " << endl;
	cout << "EndEtaScan"   << endl;
	cout << "etaMin has to be smaller than etaMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(TString(token[0]) == TString("AzimuthalAngleValues:") && phiArr.size() == 0) {
      if(n <= 2) {
	cout << endl;
	cout << "The azimuthal angle values (\"AzimuthalAngleValues:\") have to have at least 3 tokens, in the following way:" << endl;
	cout << " AzimuthalAngleValues:   val1  val2 ... Units" << endl;
	cout << "Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      
      if(!global->IsAngleUnit(token[n-1])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
      double UnitValue = global->GetAngleUnit(TString(token[n-1]));
      phiArr.clear();
      for(int kt=1;kt<n-1;kt++) phiArr.push_back(atof(token[kt])*UnitValue);
    }
    else if(n == 1 && TString(token[0]) == TString("BeginAzimuthalAngleScan") && phiArr.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 3;
      
      int Nbins = -1;
      double R[2];
      R[0] = Dummy_value;
      R[1] = Dummy_value;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndAzimuthalAngleScan")) {
	  if(Nbins <= 0)   continue;
	  if(R[0] < 0)     continue;
	  if(R[1] < 0)     continue;
	  if(R[0] >= R[1]) continue;
	  
	  IsEnd = true;
	  phiArr.clear();
	  for(int i=0;i<Nbins+1;i++) phiArr.push_back(R[0] + i*(R[1] - R[0])/Nbins);
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Nbins"))  Nbins = atoi(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("phiMin")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[0]  = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("phiMax")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  R[1]  = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of azimuthal angle scan parameters inside the BeginAzimuthalAngleScan and EndAzimuthalAngleScan block. Parameters have to be specified as follows:" << endl;
	cout << "BeginAzimuthalAngleScan" << endl;
	cout << "  Nbins     N           // N > 0"    << endl;
	cout << "  phiMin    val  units  // val >= 0" << endl;
	cout << "  phiMax    val  units  // val >= 0" << endl;
	cout << "EndAzimuthalAngleScan"   << endl;
	cout << "phiMin has to be smaller than phiMax" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginTrkResolAnalysisParams") && TrkResolAnalysisPars.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max = 20;
      
      //Initial values of parameters for track resolution analysis
      TrkResolAnalysisPars_t  aTrkResolAnalysisPars;
      aTrkResolAnalysisPars.NhitsMin                = Dummy_value;
      aTrkResolAnalysisPars.SameRange               = true;
      aTrkResolAnalysisPars.UseAllMomVals           = false;
      aTrkResolAnalysisPars.PlotMaterialBudget      = false;
      aTrkResolAnalysisPars.PlotDOCAatHighMom       = false;
      aTrkResolAnalysisPars.PlotDOCAvsMonFit        = false;
      aTrkResolAnalysisPars.PlotPhiAveraged         = false;
      aTrkResolAnalysisPars.PlotOnlyPhiAveraged     = false;
      aTrkResolAnalysisPars.PlotPerformancesVsTheta = false;
      aTrkResolAnalysisPars.UseLogYAxes             = false;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndTrkResolAnalysisParams")) {
	  IsEnd = true;
	  	  
	  if(aTrkResolAnalysisPars.NhitsMin == Dummy_value) {
	    cout << endl;
	    cout << "ERROR inside TrkResolAnalysisParams Block: NhitsMin not specified!!!"  << endl;
	    cout << endl;
	    IsEnd = false;
	  }
	  
	  if(IsEnd) TrkResolAnalysisPars.push_back(aTrkResolAnalysisPars);
	}
	else if(n1 == 2 && TString(token1[0]) == TString("NhitsMin")) aTrkResolAnalysisPars.NhitsMin = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("SameRange")) {
	  aTrkResolAnalysisPars.SameRange = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("UseAllMomVals")) {
	  aTrkResolAnalysisPars.UseAllMomVals = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotMaterialBudget")) {
	  aTrkResolAnalysisPars.PlotMaterialBudget = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotDOCAatHighMom")) {
	  aTrkResolAnalysisPars.PlotDOCAatHighMom = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotDOCAvsMonFit")) {
	  aTrkResolAnalysisPars.PlotDOCAvsMonFit = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotPhiAveraged")) {
	  aTrkResolAnalysisPars.PlotPhiAveraged = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotOnlyPhiAveraged")) {
	  aTrkResolAnalysisPars.PlotOnlyPhiAveraged = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotPerformancesVsTheta")) {
	  aTrkResolAnalysisPars.PlotPerformancesVsTheta = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("UseLogY")) {
	  aTrkResolAnalysisPars.UseLogYAxes = global->SetBoolFromString(TString(token1[1]));
	}
		
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of track parameters resolution calculation parameters inside the BeginTrkResolAnalysisParams and EndTrkResolAnalysisParams block. Parameters have to be specified as follows:" << endl;
	PrintTrkResolAnalysisParamsBlock();
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginEfficAnalysisParams") && EfficAnalysisPars.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max = 20;
      
      //Initial values of parameters for track resolution analysis
      EfficAnalysisPars_t  aEfficAnalysisPars;
      aEfficAnalysisPars.SameRange               = true;
      aEfficAnalysisPars.UseAllMomVals           = false;
      aEfficAnalysisPars.PlotPhiAveraged         = false;
      aEfficAnalysisPars.PlotOnlyPhiAveraged     = false;
      aEfficAnalysisPars.PlotPerformancesVsTheta = false;
      aEfficAnalysisPars.UseLogYAxes             = false;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndEfficAnalysisParams")) {
	  IsEnd = true;
	  
	  if(IsEnd) EfficAnalysisPars.push_back(aEfficAnalysisPars);
	}
	else if(n1 == 2 && TString(token1[0]) == TString("SameRange")) {
	  aEfficAnalysisPars.SameRange = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("UseAllMomVals")) {
	  aEfficAnalysisPars.UseAllMomVals = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotPhiAveraged")) {
	  aEfficAnalysisPars.PlotPhiAveraged = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotOnlyPhiAveraged")) {
	  aEfficAnalysisPars.PlotOnlyPhiAveraged = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("PlotPerformancesVsTheta")) {
	  aEfficAnalysisPars.PlotPerformancesVsTheta = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("UseLogY")) {
	  aEfficAnalysisPars.UseLogYAxes = global->SetBoolFromString(TString(token1[1]));
	}
		
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of track parameters resolution calculation parameters inside the BeginEfficAnalysisParams and EndEfficAnalysisParams block. Parameters have to be specified as follows:" << endl;
	PrintEfficAnalysisParamsBlock();
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginMatBudgetAnalysisParams") && MatBudgetAnalysisPars.size() == 0) {
      bool IsEnd = false;
      int counter = 0;
      int Max = 6;
      
      //Initial values of parameters for material budget analysis
      MatBudgetAnalysisPars_t aMatBudgetAnalysisPars;
      aMatBudgetAnalysisPars.mom_min =  1.0*global->GetMomentumUnit("GeV/c");
      aMatBudgetAnalysisPars.mom_max = 10.0*global->GetMomentumUnit("GeV/c");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndMatBudgetAnalysisParams")) {
	  IsEnd = true;

	  if(aMatBudgetAnalysisPars.mom_min >= aMatBudgetAnalysisPars.mom_max) {
	    cout << endl;
	    cout << "ERROR inside MatBudgetAnalysisPars Block: mom_min >= mom_max!!!"  << endl;
	    cout << endl;
	    IsEnd = false;
	  }
	  
	  if(IsEnd)  MatBudgetAnalysisPars.push_back(aMatBudgetAnalysisPars);
	}
	else if(n1 == 3 && TString(token1[0]) == TString("mom_min")) {
	  if(!global->IsMomentumUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  aMatBudgetAnalysisPars.mom_min = atof(token1[1])*global->GetMomentumUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("mom_max")) {
	  if(!global->IsMomentumUnit(token1[2])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	  aMatBudgetAnalysisPars.mom_max = atof(token1[1])*global->GetMomentumUnit(TString(token1[2]));
	}

	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of track parameters resolution calculation parameters inside the BeginMatBudgetAnalysisParams and EndMatBudgetAnalysisParams block. Parameters have to be specified as follows:" << endl;
	PrintMatBudgetAnalysisParamsBlock();
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginPseudoEfficParams")) {
      bool IsEnd = false;
      int counter = 0;
      int Max = 7;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndPseudoEfficParams")) IsEnd = true;
	else if(n1 == 2 && TString(token1[0]) == TString("NhitsMin"))  Nhits_min_sel              = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("NhitsSeed")) Nmin_Layers_track_seed_sel = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("NfakesMax")) Nfakes_max_sel             = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("ndf"))       ndf_sel                    = atoi(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("chi2Ondf"))  chi2Ondf_sel               = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("SameRange")) {
	  SameRange_sel = global->SetBoolFromString(TString(token1[1]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("SeedExternal")) {
	  SeedExternal_sel = global->SetBoolFromString(TString(token1[1]));
	}
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of Pseudo-efficiency calculation parameters inside the BeginPseudoEfficParams and EndPseudoEfficParams block. Parameters have to be specified as follows:" << endl;
	cout << "BeginPseudoEfficParams" << endl;
	cout << "  NhitsMin       N     // N   >  0"    << endl;
	cout << "  NhitsSeed      N     // N   >  0"    << endl;
	cout << "  NfakesMax      N     // N   >= 0"    << endl;
	cout << "  ndf            val   // val >= 0"    << endl;
	cout << "  chi2Ondf       val   // val >= 0"    << endl;
	cout << "  SameRange      bool  // either = true or false" << endl;
	cout << "  SeedExternal   bool  // either = true or false" << endl;
	cout << "EndPseudoEfficParams"   << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }

      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoRanges")) {
      bool IsEndGeoRanges = false;
      int counter = 0;
      int Max     = 100;
      
      int counter_georange = 0;
      
      while(!IsEndGeoRanges && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	counter++;
	
	if(n1 == 1 && TString(token1[0]) == TString("EndGeoRanges")) IsEndGeoRanges = true;
	else if(n1 == 1 && TString(token1[0]) == TString("BeginRange")) {
	  bool IsEndRange = false;
	  int counter2 = 0;
          int Max2     = 3;
	  
	  counter_georange++;
	  
	  ARange_t aRange;
	  aRange.Xmin = -1.0e+10;
	  aRange.Xmax = -1.0e+10;
	  aRange.Ymax = -1.0e+10;
	  aRange.Ymax = -1.0e+10;
	  aRange.Zmin = -1.0e+10;
	  aRange.Zmax = -1.0e+10;
	  
          while(!IsEndRange && counter2 <= Max2 && !fin.eof()) {
	    char buf2[MAX_CHARS_PER_LINE];
	    fin.getline(buf2, MAX_CHARS_PER_LINE);
	    line_number++;
	    int n2 = 0;
	    const char* token2[MAX_TOKENS_PER_LINE] = {};
	    
	    token2[0] = strtok(buf2, DELIMITER); // first token
	    if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
	    
	    for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
	      token2[n2] = strtok(0, DELIMITER);
	      if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
	    }
	    
	    if(n2 == 1 && TString(token2[0]) == TString("EndRange")) {
	      std::vector<TString> NonSetRangeParams;
	      NonSetRangeParams.clear();
	      bool GoodRange = true;
	      
	      if(aRange.Xmin == -1.0e+10) {
		NonSetRangeParams.push_back(TString("Xmin"));
		GoodRange = false;
	      }
	      if(aRange.Xmax == -1.0e+10) {
		NonSetRangeParams.push_back(TString("Xmax"));
		GoodRange = false;
	      }
	      if(aRange.Ymin == -1.0e+10) {
		NonSetRangeParams.push_back(TString("Ymin"));
		GoodRange = false;
	      }
	      if(aRange.Ymax == -1.0e+10) {
		NonSetRangeParams.push_back(TString("Ymax"));
		GoodRange = false;
	      }
	      if(aRange.Zmin == -1.0e+10) {
		NonSetRangeParams.push_back(TString("Zmin"));
		GoodRange = false;
	      }
	      if(aRange.Zmax == -1.0e+10) {
		NonSetRangeParams.push_back(TString("Zmax"));
		GoodRange = false;
	      }

	      if(aRange.Xmin >= aRange.Xmax) {
		cout << endl;
		cout << "For GeoRange " << counter_georange << ", Xmin value (" << aRange.Xmin/global->GetDistanceUnit("cm") << " cm) is higher than Xmax value (" << aRange.Xmax/global->GetDistanceUnit("cm") << " cm). " 
		     << "Check your inputs. Exiting now !!!"
		     << endl;
		cout << endl;
		assert(false);
	      }
	      else if(aRange.Ymin >= aRange.Ymax) {
		cout << endl;
		cout << "For GeoRange " << counter_georange << ", Ymin value (" << aRange.Ymin/global->GetDistanceUnit("cm") << " cm) is higher than Ymax value (" << aRange.Ymax/global->GetDistanceUnit("cm") << " cm). " 
		     << "Check your inputs. Exiting now !!!"
		     << endl;
		cout << endl;
		assert(false);
	      }
	      else if(aRange.Zmin >= aRange.Zmax) {
		cout << endl;
		cout << "For GeoRange " << counter_georange << ", Zmin value (" << aRange.Zmin/global->GetDistanceUnit("cm") << " cm) is higher than Zmax value (" << aRange.Zmax/global->GetDistanceUnit("cm") << " cm). " 
		     << "Check your inputs. Exiting now !!!"
		     << endl;
		cout << endl;
		assert(false);
	      }

	      if(GoodRange) GeoRanges.push_back(aRange);
	      else {
		cout << endl;
		cout << "WARNING:: Range parameters ";
		for(int irangepar=0;irangepar<int(NonSetRangeParams.size());irangepar++) {
		  cout << NonSetRangeParams[irangepar].Data();
		  if(irangepar < int(NonSetRangeParams.size()) - 1) cout << ", ";
		  else                                              cout << "  ";
		}
		cout << " for range " << counter_georange << " not set. Range not included in geometry range list!!!";
		cout << endl;
	      }
	      
	      IsEndRange = true;
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("XRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aRange.Xmin = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aRange.Xmax = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("YRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aRange.Ymin = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aRange.Ymax = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("ZRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aRange.Zmin = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aRange.Zmax = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    
	    counter2++;
	  }
	  
	  if(!IsEndRange) {
	    cout << endl;
	    cout << "Wrong specification of geometry range parameters inside the BeginRange and EndRange block. Parameters have to be specified as follows:" << endl;
	    cout << "BeginRange" << endl;
            cout << "  XRange    Xmin   Xmax  units" << endl;
	    cout << "  YRange    Ymin   Ymax  units" << endl;
	    cout << "  ZRange    Zmin   Zmax  units" << endl;
	    cout << "EndRange" << endl;
	    cout << "Check you inputs. Exiting now!!!" << endl;
	    cout << endl;
	    
	    assert(false);
	  }
	  continue;
	  
	}
      }
      
      if(!IsEndGeoRanges) {
	cout << endl;
	cout << "Wrong specification of geometry ranges parameters inside the BeginGeoRanges and EndGeoRanges block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoRanges" << endl;
	cout << "  BeginRange" << endl;
        cout << "    XRange    Xmin   Xmax  units" << endl;
	cout << "    YRange    Ymin   Ymax  units" << endl;
	cout << "    ZRange    Zmin   Zmax  units" << endl;
	cout << "  EndRange" << endl;
	cout << "  ... " << endl;
	cout << "  ... " << endl;
	cout << "  ... " << endl;
	cout << "  BeginRange" << endl;
        cout << "    XRange    Xmin   Xmax  units" << endl;
	cout << "    YRange    Ymin   Ymax  units" << endl;
	cout << "    ZRange    Zmin   Zmax  units" << endl;
	cout << "  EndRange" << endl;
	cout << "EndGeoRanges" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginVoxeling")) {
      bool IsEndVoxeling = false;
      int  counter = 0;
      int  Max     = 100;
      
      while(!IsEndVoxeling && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	counter++;
	
	if(n1 == 1 && TString(token1[0]) == TString("EndVoxeling")) IsEndVoxeling = true;
	else if(n1 == 1 && TString(token1[0]) == TString("BeginVoxel")) {
	  bool IsEndVoxel = false;
	  int counter2 = 0;
          int Max2     = 5;
	  
	  Voxel_t aVoxel;
	  aVoxel.Rx[0]     = -100.0*global->GetDistanceUnit("m");
	  aVoxel.Rx[1]     = +100.0*global->GetDistanceUnit("m");
	  aVoxel.Ry[0]     = -100.0*global->GetDistanceUnit("m");
	  aVoxel.Ry[1]     = +100.0*global->GetDistanceUnit("m");
	  aVoxel.Rz[0]     = -100.0*global->GetDistanceUnit("m");
	  aVoxel.Rz[1]     = +100.0*global->GetDistanceUnit("m");
	  aVoxel.Rtheta[0] =    0.0*global->GetAngleUnit("deg");
	  aVoxel.Rtheta[1] =  180.0*global->GetAngleUnit("deg");
	  aVoxel.Rphi[0]   =    0.0*global->GetAngleUnit("deg");
	  aVoxel.Rphi[1]   =  360.0*global->GetAngleUnit("deg");
	  aVoxel.Rr[0]     =    0.0;
	  aVoxel.Rr[1]     = +100.0*global->GetDistanceUnit("m");
	  
          while(!IsEndVoxel && counter2 <= Max2 && !fin.eof()) {
	    char buf2[MAX_CHARS_PER_LINE];
	    fin.getline(buf2, MAX_CHARS_PER_LINE);
	    line_number++;
	    int n2 = 0;
	    const char* token2[MAX_TOKENS_PER_LINE] = {};
	    
	    token2[0] = strtok(buf2, DELIMITER); // first token
	    if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
	    
	    for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
	      token2[n2] = strtok(0, DELIMITER);
	      if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
	    }
	    
	    if(n2 == 1 && TString(token2[0]) == TString("EndVoxel")) {
	      bool GoodVoxel = true;
	      
	      if(counter2 == 0) {
		cout << "No ranges specified for this voxel. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      } 
	      else if(aVoxel.Rx[0] >= aVoxel.Rx[1]) {
		cout << "Voxel xMin >= xMax. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      else if(aVoxel.Ry[0] >= aVoxel.Ry[1]) {
		cout << "Voxel yMin >= yMax. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      else if(aVoxel.Rz[0] >= aVoxel.Rz[1]) {
		cout << "Voxel zMin >= zMax. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      else if(aVoxel.Rr[0] < 0.0) {
		cout << "Voxel rMin < 0. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      else if(aVoxel.Rr[1] < 0.0) {
		cout << "Voxel rMax < 0. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      else if(aVoxel.Rr[0] >= aVoxel.Rr[1]) {
		cout << "Voxel rMin >= rMax. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      //else if(aVoxel.Rtheta[0] < 0.0 || aVoxel.Rtheta[0] > 180.0*global->GetAngleUnit("deg")) {
		//cout << "Voxel thetaMin not in (0,180) deg range. Voxel not included in voxel list" << endl;
		//GoodVoxel = false;
	      //}
	      //else if(aVoxel.Rtheta[1] < 0.0 || aVoxel.Rtheta[1] > 180.0*global->GetAngleUnit("deg")) {
		//cout << "Voxel thetaMax not in (0,180) deg range. Voxel not included in voxel list" << endl;
		//GoodVoxel = false;
	      //}
	      else if(aVoxel.Rtheta[0] >= aVoxel.Rtheta[1]) {
		cout << "Voxel thetaMin >= thetaMax. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }
	      //else if(aVoxel.Rphi[0] < 0.0 || aVoxel.Rphi[0] > 180.0*global->GetAngleUnit("deg")) {
		//cout << "Voxel phiMin not in (0,360) deg range. Voxel not included in voxel list" << endl;
		//GoodVoxel = false;
	      //}
	      //else if(aVoxel.Rphi[1] < 0.0 || aVoxel.Rphi[1] > 180.0*global->GetAngleUnit("deg")) {
		//cout << "Voxel phiMax not in (0,360) deg range. Voxel not included in voxel list" << endl;
		//GoodVoxel = false;
	      //}
	      else if(aVoxel.Rphi[0] >= aVoxel.Rphi[1]) {
		cout << "Voxel phiMin >= phiMax. Voxel not included in voxel list!!!" << endl;
		GoodVoxel = false;
	      }

	      if(GoodVoxel) VoxelList.push_back(aVoxel);
	      
	      IsEndVoxel = true;
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("xRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aVoxel.Rx[0] = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aVoxel.Rx[1] = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("yRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aVoxel.Ry[0] = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aVoxel.Ry[1] = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("zRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aVoxel.Rz[0] = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aVoxel.Rz[1] = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("rRange")) {
	      if(!global->IsDistanceUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aVoxel.Rr[0] = atof(token2[1])*global->GetDistanceUnit(TString(token2[3]));
	      aVoxel.Rr[1] = atof(token2[2])*global->GetDistanceUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("thetaRange")) {
	      if(!global->IsAngleUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aVoxel.Rtheta[0] = atof(token2[1])*global->GetAngleUnit(TString(token2[3]));
	      aVoxel.Rtheta[1] = atof(token2[2])*global->GetAngleUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("phiRange")) {
	      if(!global->IsAngleUnit(token2[3])) cout << "ERROR in DataCard " << DataCard.Data() << " line number " << line_number << endl;
	      
	      aVoxel.Rphi[0] = atof(token2[1])*global->GetAngleUnit(TString(token2[3]));
	      aVoxel.Rphi[1] = atof(token2[2])*global->GetAngleUnit(TString(token2[3]));
	    }
	    
	    counter2++;
	  }
	  
	  if(!IsEndVoxel) {
	    cout << endl;
	    cout << "Wrong specification of voxel parameters inside the BeginVoxel and EndVoxel block. Parameters have to be specified as follows:" << endl;
	    cout << "BeginVoxel" << endl;
	    cout << "  //At least one of the following parameters has to be specified" << endl;
	    cout << "  xRange        Xmin       Xmax       units // Xmin     < Xmax"     << endl;
	    cout << "  yRange        Ymin       Ymax       units // Ymin     < Ymax"     << endl;
	    cout << "  zRange        Zmin       Zmax       units // Zmin     < Zmax"     << endl;
	    cout << "  thetaRange    Thetamin   Thetamax   units // Thetamin < Thetamax" << endl;
	    cout << "  phiRange      Phimin     Phimax     units // Phimin   < Phimax"   << endl;
	    cout << "EndVoxel" << endl;
	    cout << "Check you inputs. Exiting now!!!" << endl;
	    cout << endl;
	    
	    assert(false);
	  }
	  continue;
	}
	
      }
      
      if(!IsEndVoxeling) {
	cout << endl;
	cout << "Wrong specification of Voxels parameters inside the BeginVoxeling and EndVoxeling block. Parameters have to be specified as follows:" << endl;
	cout << "BeginVoxeling" << endl;
	cout << "  BeginVoxel" << endl;
	cout << "    //At least one of the following parameters has to be specified" << endl;
	cout << "    xRange        Xmin       Xmax       units // Xmin     < Xmax"     << endl;
	cout << "    yRange        Ymin       Ymax       units // Ymin     < Ymax"     << endl;
	cout << "    zRange        Zmin       Zmax       units // Zmin     < Zmax"     << endl;
	cout << "    thetaRange    Thetamin   Thetamax   units // Thetamin < Thetamax" << endl;
	cout << "    phiRange      Phimin     Phimax     units // Phimin   < Phimax"   << endl;
	cout << "  EndVoxel" << endl;
	cout << "  ... " << endl;
	cout << "  ... " << endl;
	cout << "  ... " << endl;
	cout << "  BeginVoxel" << endl;
	cout << "    //At least one of the following parameters has to be specified" << endl;
	cout << "    xRange        Xmin       Xmax       units // Xmin     < Xmax"     << endl;
	cout << "    yRange        Ymin       Ymax       units // Ymin     < Ymax"     << endl;
	cout << "    zRange        Zmin       Zmax       units // Zmin     < Zmax"     << endl;
	cout << "    thetaRange    Thetamin   Thetamax   units // Thetamin < Thetamax" << endl;
	cout << "    phiRange      Phimin     Phimax     units // Phimin   < Phimax"   << endl;
	cout << "  EndVoxel" << endl;
	cout << "EndVoxeling" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 1 && TString(token[0]).BeginsWith("Begin") && TString(token[0]).Contains("Bfield") && GlobalBfield == NULL) {
      ReadBField(TString("Global B-field"),TString(token[0]),&fin,GlobalBfield,DataCard.Data(),line_number);
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeometries")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 100;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	if(n1 == 1 && TString(token1[0]) == TString("EndGeometries")) IsEnd = true;
	else if(n1 == 1) GeometryDataCardList.push_back(TString(token1[0]));
	counter++;
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of parameters inside the BeginGeometries and EndGeometries block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeometries" << endl;
        cout << "  Fullpath/GeometryDataCard_1.txt" << endl;
        cout << " ... " << endl;
	cout << "  Fullpath/GeometryDataCard_n.txt" << endl;
	cout << "EndGeometries" << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	
	assert(false);
      }
      
      continue;
    }
    else if(n == 2 && TString(token[0]) == TString("OutputFile:")) {
      TheOutputFile = TString(token[1]);
    }
  }

  cout << endl;
  if(verbose) cout << "verbose = true." << endl;
  else        cout << "verbose = false." << endl;
  cout << endl;

  TString  PartOriginUnit("cm");
  
  cout << endl;
  if(particle == Dummy_name) {
    particle = TString("pi+");
    cout << "particle type not specified. Setting it to default value:" << endl;
  }
  else if(!global->CheckParticleName(particle)) {
    particle = TString("pi+");
    cout << "Setting it to default value:" << endl;
  }
  cout << "particle        = " << particle.Data() << ", mass = " << global->GetParticleMass(particle)/global->GetMassUnit("GeV/c2") << " GeV/c^2" << endl;
  
  if(ParticleOrigin == Dummy_vector) {
    ParticleOrigin = TVector3(0.0,0.0,0.0);
    cout << "particle origin not specified. Setting it to default value:" << endl;
  }
  TVector3 ParticleOrigin_tmp = (1.0/global->GetDistanceUnit(PartOriginUnit))*ParticleOrigin;
  cout << "particle origin = (" << ParticleOrigin_tmp.X() << "," << ParticleOrigin_tmp.Y() << "," << ParticleOrigin_tmp.Z() << ") " << PartOriginUnit.Data() << endl;
  
  cout << endl;
  
  cout << endl;
  if(ReferencePoint == Dummy_vector) {
    ReferencePoint = TVector3(0.0,0.0,0.0);
    cout << "Reference point for tracking not specified. Setting it to default value:" << endl;
  }
  TVector3 ReferencePoint_tmp = (1.0/global->GetDistanceUnit(PartOriginUnit))*ReferencePoint;
  cout << "Reference point for tracking = (" << ReferencePoint_tmp.X() << "," << ReferencePoint_tmp.Y() << "," << ReferencePoint_tmp.Z() << ") " << PartOriginUnit.Data() << endl;
  cout << endl;

  cout << endl;
  if(FitPowerForImpactParam) cout << "FitPowerForImpactParam = true" << endl;
  else                       cout << "FitPowerForImpactParam = false" << endl;
  cout << endl;

  cout << endl;
  cout << "MonResolRepresentation = " << MonResolRepresentation.Data() << endl;
  cout << endl;
  
  bool SomeNegative = false;
  //MOMENTUM VALUES
  if(momArr.size() == 0) {
    cout << endl;
    cout << "WARNNING:  No polar momentum values specified. Uning single value of mom = 2.0 GeV/c!!!" << endl;
    cout << endl;
    momArr.push_back(2.0*global->GetMomentumUnit("GeV/c"));
  }
  SomeNegative = false;
  for(int k=0;k<int(momArr.size());k++) {
    if(momArr[k] < 0.0) {
      SomeNegative = true;
      break;
    }
  }
  if(SomeNegative) {
    cout << endl;
    cout << "Some of the momentum values specified are negative:" << endl;
    for(int k=0;k<int(momArr.size());k++) {
      if(momArr[k] < 0) cout << k+1 << " mom-value = " << momArr[k]/global->GetMomentumUnit("GeV/c") << " GeV/c" << endl;
    }
    cout << "Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  //Order momentum values from low to high
  for(int iii=2;iii<=int(momArr.size());iii++) {
    for(int jjj=0;jjj<=int(momArr.size())-iii;jjj++) {
      double val_jjj   = momArr[jjj];
      double val_jjjp1 = momArr[jjj+1];
      
      if(val_jjj > val_jjjp1) {
	double val_aux = momArr[jjj];
	momArr[jjj]    = momArr[jjj+1];
	momArr[jjj+1]  = val_aux;
      }
    }
  }
  cout << "Momentum values (GeV/c):       ";
  for(int k=0;k<int(momArr.size());k++) {
    cout << momArr[k]/global->GetMomentumUnit("GeV/c");
    if(k < int(momArr.size()) -1) cout << ",  ";
  }
  cout << endl;
  //cout << endl;
  
  //POLAR ANGLE VALUES
  if(thetaArr.size() == 0) {
    cout << endl;
    cout << "WARNNING:  No polar polar angle values specified. Uning single value of theta = 0.0 deg!!!" << endl;
    cout << endl;
    thetaArr.push_back(0.0*global->GetAngleUnit("deg"));
  }
  SomeNegative = false;
  for(int k=0;k<int(thetaArr.size());k++) {
    if(thetaArr[k] < 0.0) {
      SomeNegative = true;
      break;
    }
  }
  if(SomeNegative) {
    cout << endl;
    cout << "Some of the polar angle values specified are negative (allowed range from 0 -> 180 deg):" << endl;
    for(int k=0;k<int(thetaArr.size());k++) {
      if(thetaArr[k] < 0) cout << k+1 << " theta-value = " << thetaArr[k]/global->GetAngleUnit("deg") << " deg" << endl;
    }
    cout << "Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  //Order polar angle values from low to high
  for(int iii=2;iii<=int(thetaArr.size());iii++) {
    for(int jjj=0;jjj<=int(thetaArr.size())-iii;jjj++) {
      double val_jjj   = thetaArr[jjj];
      double val_jjjp1 = thetaArr[jjj+1];
      
      if(val_jjj > val_jjjp1) {
	double val_aux  = thetaArr[jjj];
	thetaArr[jjj]   = thetaArr[jjj+1];
	thetaArr[jjj+1] = val_aux;
      }
    }
  }
  
  if(polarVariable == TString("theta"))          cout << "Polar angle values (deg):      ";
  else  if(polarVariable == TString("costheta")) cout << "costheta values :              ";
  else  if(polarVariable == TString("eta"))      cout << "eta values :                   ";
  for(int k=0;k<int(thetaArr.size());k++) {
    int idx = k;
    if(polarVariable == TString("theta"))          idx = k;
    else  if(polarVariable == TString("costheta")) idx = thetaArr.size() - k - 1;
    else  if(polarVariable == TString("eta"))      idx = thetaArr.size() - k - 1;
    
    double value = thetaArr[idx]/global->GetAngleUnit("deg");
    if(polarVariable == TString("theta"))          value = thetaArr[idx]/global->GetAngleUnit("deg");
    else  if(polarVariable == TString("costheta")) value = global->FromThetaToCosTheta(thetaArr[idx]);
    else  if(polarVariable == TString("eta"))      value = global->FromThetaToEta(thetaArr[idx]);
    
    cout << value;
    if(k < int(thetaArr.size()) -1) cout << ",  ";
  }
  cout << endl;
  
  
  //AZIMUTHAL ANGLE VALUES
  if(phiArr.size() == 0) {
    cout << endl;
    cout << "WARNNING:  No azimuthal angle values specified. Uning single value of phi = 0.0 deg!!!" << endl;
    cout << endl;
    phiArr.push_back(0.0*global->GetAngleUnit("deg"));
  }
  SomeNegative = false;
  for(int k=0;k<int(phiArr.size());k++) {
    if(phiArr[k] < 0.0) {
      SomeNegative = true;
      break;
    }
  }
  if(SomeNegative) {
    cout << endl;
    cout << "Some of the polar angle values specified are negative (allowed range from 0 -> 360 deg):" << endl;
    for(int k=0;k<int(phiArr.size());k++) {
      if(phiArr[k] < 0) cout << k+1 << " phi-value = " << phiArr[k]/global->GetAngleUnit("deg") << " deg" << endl;
    }
    cout << "Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  //Order polar angle values from low to high
  for(int iii=2;iii<=int(phiArr.size());iii++) {
    for(int jjj=0;jjj<=int(phiArr.size())-iii;jjj++) {
      double val_jjj   = phiArr[jjj];
      double val_jjjp1 = phiArr[jjj+1];
      
      if(val_jjj > val_jjjp1) {
	double val_aux = phiArr[jjj];
	phiArr[jjj]    = phiArr[jjj+1];
	phiArr[jjj+1]  = val_aux;
      }
    }
  }
  cout << "Azimuthal angle values (deg):  ";
  for(int k=0;k<int(phiArr.size());k++) {
    cout << phiArr[k]/global->GetAngleUnit("deg");
    if(k < int(phiArr.size()) -1) cout << ",  ";
  }
  cout << endl;
  cout << endl;
  
  if(VoxelList.size() > 0) {
    cout << endl;
    cout << "Voxels parameters:" << endl;
    for(int ivoxel=0;ivoxel<int(VoxelList.size());ivoxel++) {
      TString VoxelDiskUnit("cm");
      TString VoxelAngleUnit("deg");
      cout << "  Begin Voxel " << ivoxel+1 << endl;
      cout << "    xRange     (" << VoxelList[ivoxel].Rx[0]/global->GetDistanceUnit(VoxelDiskUnit)   << "," << VoxelList[ivoxel].Rx[1]/global->GetDistanceUnit(VoxelDiskUnit)   << ") " << VoxelDiskUnit.Data() << endl;
      cout << "    yRange     (" << VoxelList[ivoxel].Ry[0]/global->GetDistanceUnit(VoxelDiskUnit)   << "," << VoxelList[ivoxel].Ry[1]/global->GetDistanceUnit(VoxelDiskUnit)   << ") " << VoxelDiskUnit.Data() << endl;
      cout << "    zRange     (" << VoxelList[ivoxel].Rz[0]/global->GetDistanceUnit(VoxelDiskUnit)   << "," << VoxelList[ivoxel].Rz[1]/global->GetDistanceUnit(VoxelDiskUnit)   << ") " << VoxelDiskUnit.Data() << endl;
      cout << "    thetaRange (" << VoxelList[ivoxel].Rtheta[0]/global->GetAngleUnit(VoxelAngleUnit) << "," << VoxelList[ivoxel].Rtheta[1]/global->GetAngleUnit(VoxelAngleUnit) << ") " << VoxelAngleUnit.Data() << endl;
      cout << "    phiRange   (" << VoxelList[ivoxel].Rphi[0]/global->GetAngleUnit(VoxelAngleUnit)   << "," << VoxelList[ivoxel].Rphi[1]/global->GetAngleUnit(VoxelAngleUnit)   << ") " << VoxelAngleUnit.Data() << endl;            
      cout << "  End   Voxel " << ivoxel+1 << endl;
      cout << endl;
    }
    cout << endl;
  }
  
  
  if(GeoRanges.size() > 0) {
    cout << endl;
    cout << "Visualization Geometry ranges:" << endl;
    for(int igeorange=0;igeorange<int(GeoRanges.size());igeorange++) {
      TString UnitRange("cm");
      cout << "  Begin Range " << igeorange+1 << endl;
      cout << "    XRange (" << GeoRanges[igeorange].Xmin/global->GetDistanceUnit(UnitRange) << "," << GeoRanges[igeorange].Xmax/global->GetDistanceUnit(UnitRange) << ") " << UnitRange.Data() << endl;
      cout << "    YRange (" << GeoRanges[igeorange].Ymin/global->GetDistanceUnit(UnitRange) << "," << GeoRanges[igeorange].Ymax/global->GetDistanceUnit(UnitRange) << ") " << UnitRange.Data() << endl;
      cout << "    ZRange (" << GeoRanges[igeorange].Zmin/global->GetDistanceUnit(UnitRange) << "," << GeoRanges[igeorange].Zmax/global->GetDistanceUnit(UnitRange) << ") " << UnitRange.Data() << endl;
      cout << "  End   Range " << igeorange+1 << endl;
      cout << endl;
    }
    cout << endl;
  }
  
  if(MCSeed != Dummy_value) {
    cout << endl;
    cout << "MCSeed set to " << MCSeed << endl;
    cout << endl;
    global->rand->SetSeed(MCSeed);
  }
  
  cout << endl;
  if(IncludeEloss) {
    cout << "IncludeEloss         = true"  << endl;
    if(kappa != Dummy_value) global->SetKappa(kappa);
    cout << "  kappa of Eloss fluctuations = " << global->GetKappa()/global->GetEnergyUnit("GeV") << " GeV" << endl;
  }
  else                      cout << "IncludeEloss         = false" << endl;
  cout << endl;
  
  cout << endl;
  if(SavePlots)             cout << "SavePlots            = true"  << endl;
  else                      cout << "SavePlots            = false" << endl;
  if(DoGeometryCheck) {
    cout << "GeometryCheck        = true   =>  GeoCheckPrecision = " << GeoCheckPrecision/global->GetDistanceUnit("mm") << " mm" << endl;
  }
  else                      cout << "GeometryCheck        = false" << endl;
  if(DoPrintGeometry)       cout << "PrintGeometry        = true"  << endl;
  else                      cout << "PrintGeometry        = false" << endl;  
  if(DoPrintGeometryWeight) cout << "PrintGeometryWeight  = true"  << endl;
  else                      cout << "PrintGeometryWeight  = false" << endl;
  
  if(DoPlotGeometry) {
    cout << "PlotGeometry         = true"  << endl;
    if(DoPlotWorldVolume) cout << "  PlotWorldVolume       = true"  << endl;
    else                  cout << "  PlotWorldVolume       = false" << endl;
    
    if(DoPlotStepBfieldVolume) cout << "  PlotStepBfieldVolume  = true"  << endl;
    else                       cout << "  PlotStepBfieldVolume  = false" << endl;
    
    if(DoRZGeoRepresentation) cout << "  DoRZGeoRepresentation = true"  << endl;
    else                      cout << "  DoRZGeoRepresentation = false" << endl;
        
    if(DoPlotSomeTracks) {
      cout << "  PlotSomeTracks       = true"  << endl;
      if(UseAllMomVals_GeoPrint) cout << "    Plotting tracks with all specified momenta." << endl;
      else                       cout << "    Plotting tracks a sub-set of specified momenta." << endl;
    }
    else  cout << "  PlotSomeTracks       = false" << endl;
    
  }
  else  cout << "PlotGeometry         = false" << endl;
  
  if(DoMaterialBugdetAnalysis) {
    cout << "DoMaterialBugdetAnalysis   = true"  << endl;
    
    if(MatBudgetAnalysisPars.size() == 0) {
      cout << endl;
      cout << "Requested material budget analysis, but parameters not set. Need to add to datacard a MatBudgetAnalysisParams block as follows:" << endl;
      PrintMatBudgetAnalysisParamsBlock();
      cout << endl;
      assert(false);
    }
    else PrintMatBudgetAnalysisParamsBlock();
  }
  else  cout << "DoMaterialBugdetAnalysis   = false" << endl;
  
  if(DoTelescopeAnalysis) {
    cout << "DoTelescopeAnalysis  = true"  << endl;
    
    if(TrkResolAnalysisPars.size() == 0) {
      cout << endl;
      cout << "Requested track resolution analysis, but parameters not set. Need to add to datacard a TrkResolAnalysisParams block as follows:" << endl;
      PrintTrkResolAnalysisParamsBlock();
      cout << endl;
      assert(false);
    }
    else PrintTrkResolAnalysisParamsBlock();
  }
  else {
    cout << "DoTelescopeAnalysis  = false"  << endl;
    
    if(DoTrkResolAnalysis)  {
      cout << "DoTrkResolAnalysis   = true"  << endl;
    
      if(TrkResolAnalysisPars.size() == 0) {
        cout << endl;
        cout << "Requested track resolution analysis, but parameters not set. Need to add to datacard a TrkResolAnalysisParams block as follows:" << endl;
        PrintTrkResolAnalysisParamsBlock();
        cout << endl;
        assert(false);
      }
      else PrintTrkResolAnalysisParamsBlock();
    }
    else cout << "DoTrkResolAnalysis   = false" << endl;
  }
  
  if(DoPseudoEfficVsMon) {
    cout << "DoPseudoEfficVsMon   = true"  << endl;
    
    if(EfficAnalysisPars.size() == 0) {
      cout << endl;
      cout << "Requested tracking efficiency analysis, but parameters not set. Need to add to datacard a EfficAnalysisParams block as follows:" << endl;
      PrintEfficAnalysisParamsBlock();
      cout << endl;
      assert(false);
    }
    else PrintEfficAnalysisParamsBlock();
  }
  else                     cout << "DoPseudoEfficVsMon   = false" << endl;
  cout << endl;
  
  //cout << endl;
  if(GlobalBfield == NULL) {
    GlobalBfield = new GBFieldConstant("Global B-field",TVector3(1.0,0.0,0.0),0.0,global);
    //cout << "Global B-field not set. Setting it to zero." << endl;
  }
  //cout << "Global B-field:" << endl;
  //GlobalBfield->Print();
  //cout << endl;
  
  cout << endl;
  if(GeometryDataCardList.size() == 0) {
    cout << endl;
    cout << " No geometry data-cards specified. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  else {
    cout << "Printing the geometry datacards list:" << endl;
    for(int k=0;k<int(GeometryDataCardList.size());k++) {
      cout << k+1 << "  " << GeometryDataCardList[k].Data() << endl;
    }
  }
  cout << endl;
  
  cout << endl;
  
  if(TheOutputFile == Dummy_name) {
    cout << "Output file name not specified. It has to be done as follows:" << endl;
    cout << "OutputFile: DIR/filename (with no extension)" << endl;
    cout << "Check your input. Exiting now!!!" << endl;
    assert(false);
  }
  else {
    TString Directory = global->GetOutputDirectory(TheOutputFile);
    TString command = TString("mkdir -p  ") + Directory;
    gSystem->Exec(command.Data());
    
    cout << "Output Directory name: " << Directory.Data()     << endl;
    cout << "Output file name:      " << TheOutputFile.Data() << endl;
  }
  cout << endl;

  return;
  
}

//====================================================================
void Guariguanchi::FillGeometryFromDataCard(const char* datacard, int index)
{
  
  TString GeoTmpName = TString("Geometry ") + long(index);
  GGeometry* aGeometry = new GGeometry(GeoTmpName.Data(),index,global);
  
  ReadGeometryDataCard(datacard,aGeometry);

  GeometryList.push_back(aGeometry);

  return;
  
}
//====================================================================
void Guariguanchi::ReadGeometryDataCard(const char* datacard, GGeometry* aGeometry)
{
  
  //Readout of the geometry datacard
  
  // create a file-reading object
  ifstream fin;
  fin.open(datacard); // open a file
  if(!fin.good()) {
    cout << endl;
    cout << "Geometry datacard file " << datacard << " not found. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }

  int line_number = 0;
  //read each line of the file
  while(!fin.eof()) {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
    line_number++;

    // parse the line into blank-delimited tokens
    int n = 0; // a for-loop index
    
    // array to store memory addresses of the tokens in buf
    const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
    
    // parse the line
    token[0] = strtok(buf, DELIMITER); // first token
    
    // go to next line if the line is blank of it starts with comment "//"
    if(!token[0] || TString(token[0]).BeginsWith(COMMENT)) continue;
    
    for(n=1; n<MAX_TOKENS_PER_LINE; n++) {
      token[n] = strtok(0, DELIMITER); // subsequent tokens
      
      if(!token[n] || TString(token[n]).BeginsWith(COMMENT)) break; // no more tokens
    }

    if(n >= 2 && TString(token[0]) == TString("GeometryName:")) {
      TString Name("");
      for(int i=1;i<n;i++) {
	Name += TString(token[i]);
	if(i < n-1) Name += TString(" ");
      }
      aGeometry->SetName(Name);
      if(verbose) cout << "Geometry name " << aGeometry->GetName().Data() << endl;
    }
    else if(n == 2 && TString(token[0]) == TString("InputFile:")) {
      ReadGeometryDataCard(token[1],aGeometry);
    }
    else if(n == 1 && TString(token[0]).BeginsWith("Begin") && TString(token[0]).Contains("Bfield") && aGeometry->GetBField() == NULL) {
      TString BfieldName = TString("geometry ") + aGeometry->GetName() + TString(" Bfield");
      GBField* aBfield = NULL;
      ReadBField(BfieldName,TString(token[0]),&fin,aBfield,datacard,line_number);
      if(aBfield != NULL) {
	aGeometry->SetBField(aBfield);
	delete aBfield;
      }
    }
    else if(n == 2 && TString(token[0]) == TString("GlobalBkgRateScaling:")) {
      float scaling = atof(token[1]);
      if(scaling < 0.0) {
	cout << endl;
	cout << "Error in datacard " << datacard << " at line number " << line_number << " : GlobalBkgRateScaling has to be >= 0. Check your inputs. Exiting now!!!" << endl; 
	cout << endl;
	assert(false);
      }
      aGeometry->SetBkgRateScaling(scaling);
    }
    else if(n == 1 && TString(token[0]).BeginsWith("Begin") && TString(token[0]).Contains("TrackFinderAlgo")) {
      ReadTrackFinderAlgo(TString(token[0]),&fin,aGeometry,datacard,line_number);
    }
    else if(n == 1 && TString(token[0]).BeginsWith("Begin") && TString(token[0]).Contains("WorldVolume") && aGeometry->GetWorldVolume() == NULL) {
      TString WorldVolName = TString("geometry ") + aGeometry->GetName() + TString(" World Volume");
      GGeoObject* aWorldolume = NULL;
      ReadWorldVolume(WorldVolName,TString(token[0]),&fin,aWorldolume,aGeometry->GetGeometryID(),datacard,line_number);
      if(aWorldolume != NULL) {
	aGeometry->SetWorldVolume(aWorldolume);
	delete aWorldolume;
      }
    }
    else if(n == 1 && TString(token[0]) == TString("BeginSystemsConfiguration")) {
      ReadSystemsConfiguration(&fin,aGeometry,datacard,line_number);
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGasDetector")) {
      ReadGasDetector(&fin,aGeometry,datacard,line_number);
    }
    else if(n == 1 && TString(token[0]) == TString("BeginBeamTestConfiguration")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 4;
    
      std::vector<TString>  TelescopeLayers;
      std::vector<TString>  DUTLayers;
      TelescopeLayers.clear();
      DUTLayers.clear();
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
            
        if(n1 == 1 && TString(token1[0]) == TString("EndBeamTestConfiguration")) {
	  bool GoodGeoObject = true;
	  
	  if(TelescopeLayers.size() == 0) {
	    cout << "ERROR in BeamTestConfiguration Block: Telescope system list not specified not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(DUTLayers.size() == 0) {
	    cout << "ERROR in BeamTestConfiguration Block: DUT system list not specified not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(GoodGeoObject) {
	    aGeometry->FillBeamTestConfigLayers(TelescopeLayers,DUTLayers);
	    TelescopeLayers.clear();
	    DUTLayers.clear();
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("TelescopeLayersList")) {
          for(int i=1;i<n1;i++) TelescopeLayers.push_back(TString(token1[i]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("DUTLayersList")) {
          for(int i=1;i<n1;i++) DUTLayers.push_back(TString(token1[i]));
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginBeamTestConfiguration and EndBeamTestConfiguration block. Parameters have to be specified as follows:" << endl;
	cout << "BeginBeamTestConfiguration" << endl;
	cout << "  TelescopeLayersList   system1  system2  system3  system4 ... " << endl;
	cout << "  DUTLayersList         system1  system2  system3  system4 ... " << endl;
	cout << "EndBeamTestConfiguration" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoPlane")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 27;
    
      TString   Name           = Dummy_name;
      TVector3  Position       = Dummy_vector;
      TVector3  RotAngles      = Dummy_vector;
      double    Thickness      = Dummy_value;
      TString   Material       = Dummy_name;
      double    XOX0           = Dummy_value;
      double    widthU         = Dummy_value;
      double    widthV         = Dummy_value;
      bool      IsSensitive    = false;
      double    Resolution     = Dummy_value;
      double    ResolutionU    = Dummy_value;
      double    ResolutionV    = Dummy_value;
      double    ROtime         = Dummy_value;
      double    DetEffic       = Dummy_value;
      double    InsensFracUneg = 0.0;
      double    InsensFracUpos = 0.0;
      double    InsensFracVneg = 0.0;
      double    InsensFracVpos = 0.0;
      double    BkgRate        = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoPlane")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoPlane Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoPlane Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Thickness == Dummy_value) {
	    cout << "ERROR in GeoPlane Block: Thickness not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoPlane Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(widthU == Dummy_value) {
	    cout << "ERROR in GeoPlane Block: widthU not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(widthV == Dummy_value) {
	    cout << "ERROR in GeoPlane Block: widthV not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoPlane Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoPlane Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoPlane Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoPlane Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoPlane Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoPlane Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	    
	    GGeoObject* aGeoObject = new GGeoPlane(Name,aGeometry->GetGeometryID(),0,
						   Position,Rot,Thickness,Material,
					           IsSensitive,
					           TVector2(widthU,widthV),
						   TVector2(InsensFracUneg,InsensFracUpos),
						   TVector2(InsensFracVneg,InsensFracVpos),
						   global,
					           ResolutionU,ResolutionV,DetEffic,ROtime,
					           BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness       = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("widthU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  widthU          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("widthV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  widthV          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoPlane and EndGeoPlane block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoPlane" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness       val units (Mandatory)" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  widthU          val units (Mandatory)" << endl;
	cout << "  widthV          val units (Mandatory)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoPlane" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoCylinder")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 27;
    
      TString   Name           = Dummy_name;
      TVector3  Position       = Dummy_vector;
      TVector3  RotAngles      = Dummy_vector;
      double    Thickness      = Dummy_value;
      TString   Material       = Dummy_name;
      double    XOX0           = Dummy_value;
      double    Length         = Dummy_value;
      double    Radius         = Dummy_value;
      bool      IsSensitive    = false;
      double    Resolution     = Dummy_value;
      double    ResolutionU    = Dummy_value;
      double    ResolutionV    = Dummy_value;
      double    ROtime         = Dummy_value;
      double    DetEffic       = Dummy_value;
      double    InsensFracVneg = 0.0;
      double    InsensFracVpos = 0.0;
      double    BkgRate        = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoCylinder")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoCylinder Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoCylinder Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Thickness == Dummy_value) {
	    cout << "ERROR in GeoCylinder Block: Thickness not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoCylinder Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Length == Dummy_value) {
	    cout << "ERROR in GeoCylinder Block: Length not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Radius == Dummy_value) {
	    cout << "ERROR in GeoCylinder Block: Radius not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoCylinder Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoCylinder Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoCylinder Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoCylinder Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoCylinder Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoCylinder Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	  
	    GGeoObject* aGeoObject = new GGeoCylinder(Name,aGeometry->GetGeometryID(),0,
						      Position,Rot,Thickness,Material,
					              IsSensitive,
					              Length,Radius,
						      TVector2(InsensFracVneg,InsensFracVpos),
						      global,
					              ResolutionU,ResolutionV,DetEffic,ROtime,
					              BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness       = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Length          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Radius")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Radius          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoCylinder and EndGeoCylinder block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoCylinder" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness       val units (Mandatory)" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  Length          val units (Mandatory)" << endl;
	cout << "  Radius          val units (Mandatory)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoCylinder" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoCylinderSection")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 27;
    
      TString   Name           = Dummy_name;
      TVector3  Position       = Dummy_vector;
      TVector3  RotAngles      = Dummy_vector;
      double    Thickness      = Dummy_value;
      TString   Material       = Dummy_name;
      double    XOX0           = Dummy_value;
      double    Length         = Dummy_value;
      double    Radius         = Dummy_value;
      double    DeltaPhi       = Dummy_value;
      bool      IsSensitive    = false;
      double    Resolution     = Dummy_value;
      double    ResolutionU    = Dummy_value;
      double    ResolutionV    = Dummy_value;
      double    ROtime         = Dummy_value;
      double    DetEffic       = Dummy_value;
      double    InsensFracUneg = 0.0;
      double    InsensFracUpos = 0.0;
      double    InsensFracVneg = 0.0;
      double    InsensFracVpos = 0.0;
      double    BkgRate        = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoCylinderSection")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoCylinderSection Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoCylinderSection Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Thickness == Dummy_value) {
	    cout << "ERROR in GeoCylinderSection Block: Thickness not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoCylinderSection Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Length == Dummy_value) {
	    cout << "ERROR in GeoCylinderSection Block: Length not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Radius == Dummy_value) {
	    cout << "ERROR in GeoCylinderSection Block: Radius not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoCylinderSection Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(DeltaPhi == Dummy_value) {
	    cout << "WARNING in GeoCylinderSection Block: DeltaPhi not specified. Setting it to 2*pi." << endl;
	    DeltaPhi = 2*TMath::Pi();
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoCylinderSection Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoCylinderSection Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoCylinderSection Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoCylinderSection Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoCylinderSection Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	  
	    GGeoObject* aGeoObject = new GGeoCylinderSection(Name,aGeometry->GetGeometryID(),0,
							     Position,Rot,Thickness,Material,
						             IsSensitive,
					                     Length,Radius,
						             DeltaPhi,
						             TVector2(InsensFracUneg,InsensFracUpos),
						             TVector2(InsensFracVneg,InsensFracVpos),
						             global,
					                     ResolutionU,ResolutionV,DetEffic,ROtime,
					                     BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness       = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Length          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Radius")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Radius          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("DeltaPhi")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  DeltaPhi        = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoCylinderSection and EndGeoCylinderSection block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoCylinderSection" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness       val units (Mandatory)" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  Length          val units (Mandatory)" << endl;
	cout << "  Radius          val units (Mandatory)" << endl;
	cout << "  DeltaPhi        val units (Optional, default value is 2*pi)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoCylinderSection" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoDisk")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 27;
    
      TString   Name           = Dummy_name;
      TVector3  Position       = Dummy_vector;
      TVector3  RotAngles      = Dummy_vector;
      double    Thickness      = Dummy_value;
      TString   Material       = Dummy_name;
      double    XOX0           = Dummy_value;
      double    Rin            = Dummy_value;
      double    Rout           = Dummy_value;
      bool      IsSensitive    = false;
      double    Resolution     = Dummy_value;
      double    ResolutionU    = Dummy_value;
      double    ResolutionV    = Dummy_value;
      double    ROtime         = Dummy_value;
      double    DetEffic       = Dummy_value;
      double    InsensFracVneg = 0.0;
      double    InsensFracVpos = 0.0;
      double    BkgRate        = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoDisk")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoDisk Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoDisk Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Thickness == Dummy_value) {
	    cout << "ERROR in GeoDisk Block: Thickness not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoDisk Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Rin == Dummy_value) {
	    cout << "ERROR in GeoDisk Block: Inner Radius not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Rout == Dummy_value) {
	    cout << "ERROR in GeoDisk Block: Outer Radius not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoDisk Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoDisk Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoDisk Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoDisk Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoDisk Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoDisk Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	  
	    GGeoObject* aGeoObject = new GGeoDisk(Name,aGeometry->GetGeometryID(),0,
						  Position,Rot,Thickness,Material,
					          IsSensitive,
					          Rin,Rout,
						  TVector2(InsensFracVneg,InsensFracVpos),
						  global,
					          ResolutionU,ResolutionV,DetEffic,ROtime,
					          BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness       = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("Rin")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Rin             = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Rout")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Rout            = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoDisk and EndGeoDisk block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoDisk" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness       val units (Mandatory)" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  Rin             val units (Mandatory)" << endl;
	cout << "  Rout            val units (Mandatory)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoDisk" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoDiskSection")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 28;
    
      TString   Name           = Dummy_name;
      TVector3  Position       = Dummy_vector;
      TVector3  RotAngles      = Dummy_vector;
      double    Thickness      = Dummy_value;
      TString   Material       = Dummy_name;
      double    XOX0           = Dummy_value;
      double    Rin            = Dummy_value;
      double    Rout           = Dummy_value;
      double    DeltaPhi       = Dummy_value;
      bool      IsSensitive    = false;
      double    Resolution     = Dummy_value;
      double    ResolutionU    = Dummy_value;
      double    ResolutionV    = Dummy_value;
      double    ROtime         = Dummy_value;
      double    DetEffic       = Dummy_value;
      double    InsensFracUneg = 0.0;
      double    InsensFracUpos = 0.0;
      double    InsensFracVneg = 0.0;
      double    InsensFracVpos = 0.0;
      double    BkgRate        = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoDiskSection")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoDiskSection Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoDiskSection Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Thickness == Dummy_value) {
	    cout << "ERROR in GeoDiskSection Block: Thickness not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoDiskSection Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Rin == Dummy_value) {
	    cout << "ERROR in GeoDiskSection Block: Inner Radius not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Rout == Dummy_value) {
	    cout << "ERROR in GeoDiskSection Block: Outer Radius not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoDiskSection Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  if(DeltaPhi == Dummy_value) {
	    cout << "WARNING in GeoDiskSection Block: DeltaPhi not specified. Setting it to 2*pi." << endl;
	    DeltaPhi = 2.0*TMath::Pi();
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoDiskSection Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoDiskSection Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoDiskSection Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoDiskSection Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoDiskSection Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	  
	    GGeoObject* aGeoObject = new GGeoDiskSection(Name,aGeometry->GetGeometryID(),0,
							 Position,Rot,Thickness,Material,
						         IsSensitive,
					                 Rin,Rout,
						         DeltaPhi,
						         TVector2(InsensFracUneg,InsensFracUpos),
						         TVector2(InsensFracVneg,InsensFracVpos),
						         global,
					                 ResolutionU,ResolutionV,DetEffic,ROtime,
					                 BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness       = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("Rin")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Rin             = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Rout")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Rout            = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("DeltaPhi")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  DeltaPhi        = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoDiskSection and EndGeoDiskSection block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoDiskSection" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness       val units (Mandatory)" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  Rin             val units (Mandatory)" << endl;
	cout << "  Rout            val units (Mandatory)" << endl;
	cout << "  DeltaPhi        va  units (Optional). Default value is 2*pi" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracUneg  val // insensitive angular-length fraction for low-edge  U (Optional, default is zero)" << endl;
	cout << "  InsensFracUpos  val // insensitive angular-length fraction for high-edge U (Optional, default is zero)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoDiskSection" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoCone")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 27;
    
      TString   Name             = Dummy_name;
      TVector3  Position         = Dummy_vector;
      TVector3  RotAngles        = Dummy_vector;
      double    Thickness1       = Dummy_value;
      double    Thickness2       = Dummy_value;
      double    XOX0             = Dummy_value;
      TString   Material         = Dummy_name;
      double    Length           = Dummy_value;
      double    Radius1          = Dummy_value;
      double    Radius2          = Dummy_value;
      bool      IsSensitive      = false;
      double    Resolution       = Dummy_value;
      double    ResolutionU      = Dummy_value;
      double    ResolutionV      = Dummy_value;
      double    ROtime           = Dummy_value;
      double    DetEffic         = Dummy_value;
      double    InsensFracVneg   = 0.0;
      double    InsensFracVpos   = 0.0;
      double    BkgRate          = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoCone")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoCone Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoCone Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoCone Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Length == Dummy_value) {
	    cout << "ERROR in GeoCone Block: Length not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Radius1 == Dummy_value) {
	    cout << "ERROR in GeoCone Block: Radius1 not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Radius2 == Dummy_value) {
	    cout << "ERROR in GeoCone Block: Radius2 not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoCone Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(Thickness1 == Dummy_value) {
	    cout << "ERROR in GeoCone Block: Thickness1 not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else {
	    if(Thickness2 == Dummy_value) {
	      //Use constant cross section along cone
	      if(Radius2 > 1.0*global->GetDistanceUnit("um")) Thickness2 = Thickness1*(Radius1/Radius2);
	      else                                            Thickness2 = Thickness1;
	    }
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoCone Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoCone Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoCone Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoCone Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoCone Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	  
	    GGeoObject* aGeoObject = new GGeoCone(Name,aGeometry->GetGeometryID(),0,
						  Position,Rot,Thickness1,Thickness2,Material,
					          IsSensitive,
					          Length,Radius1,Radius2,
						  TVector2(InsensFracVneg,InsensFracVpos),
						  global,
					          ResolutionU,ResolutionV,DetEffic,ROtime,
					          BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness1")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness1      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Thickness2")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness2      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Length          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Radius1")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Radius1         = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Radius2")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Radius2         = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoCone and EndGeoCone block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoCone" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness1      val units (Mandatory)" << endl;
	cout << "  Thickness2      val units (Optional. If not specified will assume constant cross-section along cone => Thickness2 = Thickness1*(Radius1/Radius2))" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  Length          val units (Mandatory)" << endl;
	cout << "  Radius1         val units (Mandatory)" << endl;
	cout << "  Radius2         val units (Mandatory)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoCone" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoConeSection")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 28;
    
      TString   Name             = Dummy_name;
      TVector3  Position         = Dummy_vector;
      TVector3  RotAngles        = Dummy_vector;
      double    Thickness1       = Dummy_value;
      double    Thickness2       = Dummy_value;
      double    XOX0             = Dummy_value;
      TString   Material         = Dummy_name;
      double    Length           = Dummy_value;
      double    Radius1          = Dummy_value;
      double    Radius2          = Dummy_value;
      double    DeltaPhi         = Dummy_value;
      bool      IsSensitive      = false;
      double    Resolution       = Dummy_value;
      double    ResolutionU      = Dummy_value;
      double    ResolutionV      = Dummy_value;
      double    ROtime           = Dummy_value;
      double    DetEffic         = Dummy_value;
      double    InsensFracUneg   = 0.0;
      double    InsensFracUpos   = 0.0;
      double    InsensFracVneg   = 0.0;
      double    InsensFracVpos   = 0.0;
      double    BkgRate          = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoConeSection")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoConeSection Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoConeSection Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoConeSection Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Length == Dummy_value) {
	    cout << "ERROR in GeoConeSection Block: Length not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Radius1 == Dummy_value) {
	    cout << "ERROR in GeoConeSection Block: Radius1 not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Radius2 == Dummy_value) {
	    cout << "ERROR in GeoConeSection Block: Radius2 not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoConeSection Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(DeltaPhi == Dummy_value) {
	    cout << "WARNING in GeoConeSection Block: DeltaPhi not specified. Setting it to 2*pi." << endl;
	    DeltaPhi = 2*TMath::Pi();
	  }
	  
	  if(Thickness1 == Dummy_value) {
	    cout << "ERROR in GeoConeSection Block: Thickness1 not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else {
	    if(Thickness2 == Dummy_value) {
	      //Use constant cross section along cone
	      if(Radius2 > 1.0*global->GetDistanceUnit("um")) Thickness2 = Thickness1*(Radius1/Radius2);
	      else                                            Thickness2 = Thickness1;
	    }
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoConeSection Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoConeSection Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoConeSection Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoConeSection Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoConeSection Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	  
	    GGeoObject* aGeoObject = new GGeoConeSection(Name,aGeometry->GetGeometryID(),0,
							 Position,Rot,Thickness1,Thickness2,Material,
						         IsSensitive,
					                 Length,Radius1,Radius2,
						         DeltaPhi,
						         TVector2(InsensFracUneg,InsensFracUpos),
							 TVector2(InsensFracVneg,InsensFracVpos),
						         global,
					                 ResolutionU,ResolutionV,DetEffic,ROtime,
					                 BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness1")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness1      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Thickness2")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness2      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Length          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Radius1")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Radius1         = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Radius2")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Radius2         = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("DeltaPhi")) {
	  if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  DeltaPhi        = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPlane parameters inside the BeginGeoConeSection and EndGeoConeSection block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoConeSection" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness1      val units (Mandatory)" << endl;
	cout << "  Thickness2      val units (Optional. If not specified will assume constant cross-section along cone => Thickness2 = Thickness1*(Radius1/Radius2))" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  Length          val units (Mandatory)" << endl;
	cout << "  Radius1         val units (Mandatory)" << endl;
	cout << "  Radius2         val units (Mandatory)" << endl;
	cout << "  DeltaPhi        va  units (Optional, default value is 2*pi)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoConeSection" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginGeoPetal")) {
      bool IsEnd  = false;
      int counter = 0;
      int Max     = 28;
    
      TString   Name           = Dummy_name;
      TVector3  Position       = Dummy_vector;
      TVector3  RotAngles      = Dummy_vector;
      double    Thickness      = Dummy_value;
      TString   Material       = Dummy_name;
      double    XOX0           = Dummy_value;
      double    bottomWidth    = Dummy_value;
      double    topWidth       = Dummy_value;
      double    Height         = Dummy_value;
      bool      IsSensitive    = false;
      double    Resolution     = Dummy_value;
      double    ResolutionU    = Dummy_value;
      double    ResolutionV    = Dummy_value;
      double    ROtime         = Dummy_value;
      double    DetEffic       = Dummy_value;
      double    InsensFracUneg = 0.0;
      double    InsensFracUpos = 0.0;
      double    InsensFracVneg = 0.0;
      double    InsensFracVpos = 0.0;
      double    BkgRate        = 0.0;
      TString   SystemName("");
      TString   LayerName("");
      TString   ResolutionModel("");
      TString   EfficiencyModel("");
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
      
        if(n1 == 1 && TString(token1[0]) == TString("EndGeoPetal")) {
	  bool GoodGeoObject = true;
	  
	  if(Name == Dummy_name) {
	    cout << "ERROR in GeoPetal Block: Name not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Position == Dummy_vector) {
	    cout << "ERROR in GeoPetal Block: Position not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Thickness == Dummy_value) {
	    cout << "ERROR in GeoPetal Block: Thickness not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  if(Material == Dummy_name) {
	    cout << "ERROR in GeoPetal Block: Material not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(bottomWidth == Dummy_value) {
	    cout << "ERROR in GeoPetal Block: bottomWidth not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(topWidth == Dummy_value) {
	    cout << "ERROR in GeoPetal Block: topWidth not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  else if(Height == Dummy_value) {
	    cout << "ERROR in GeoPetal Block: Height not specified!!!" << endl;
	    GoodGeoObject = false;
	  }
	  
	  if(RotAngles == Dummy_vector) {
	    cout << "WARNING in GeoPetal Block: RotAngles not specified. Setting all to zero." << endl;
	    RotAngles = TVector3(0,0,0);
	  }
	  
	  if(IsSensitive) {
	    if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	      if(Resolution == Dummy_value) {
		cout << "ERROR in GeoPetal Block: Resolution not specified!!!" << endl;
		GoodGeoObject = false;
	      }
	      else {
		ResolutionU = Resolution;
		ResolutionV = Resolution;
	      }
	    }
	    else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	      cout << "ERROR in GeoPetal Block: ResolutionU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	      cout << "ERROR in GeoPetal Block: ResolutionV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(ROtime == Dummy_value) {
	      cout << "ERROR in GeoPetal Block: ROtime not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(DetEffic == Dummy_value) {
	      cout << "ERROR in GeoPetal Block: DetEffic not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	  }
	  
	  if(GoodGeoObject) {
	    TMatrixD Rot;
	    global->GetGlobalRotationMatrix(RotAngles,Rot);
	    
	    GGeoObject* aGeoObject = new GGeoPetal(Name,aGeometry->GetGeometryID(),0,
						   Position,Rot,Thickness,Material,
					           IsSensitive,
					           bottomWidth,topWidth,Height,
						   TVector2(InsensFracUneg,InsensFracUpos),
						   TVector2(InsensFracVneg,InsensFracVpos),
						   global,
					           ResolutionU,ResolutionV,DetEffic,ROtime,
					           BkgRate);
	    aGeoObject->SetSystemName(SystemName);
	    aGeoObject->SetLayerName(LayerName);
	    aGeoObject->SetResolutionModel(ResolutionModel);
	    aGeoObject->SetEfficiencyModel(EfficiencyModel);
	    if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	    aGeometry->PushGeoElement(aGeoObject);
	    
	    IsEnd = true;
	  }
	  
        }
        else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  Name = TString("");
          for(int i=1;i<n1;i++) {
	    Name += TString(token1[i]);
	    if(i < n1-1) Name += TString(" ");
          }          
	}
        else if(n1 == 3 && TString(token1[0]) == TString("Thickness")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Thickness       = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("XOX0"))            XOX0            = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("Material"))        Material        = TString(token1[1]);
        else if(n1 == 3 && TString(token1[0]) == TString("bottomWidth")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  bottomWidth     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("topWidth")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  topWidth        = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Height")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Height          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("Resolution")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionU")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ResolutionV")) {
	  if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	}
	else if(n1 == 3 && TString(token1[0]) == TString("ROtime")) {
	  if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
	}
	else if(n1 == 2 && TString(token1[0]) == TString("Efficiency"))      DetEffic        = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token1[1]);
	else if(n1 == 2 && TString(token1[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token1[1]);
	else if(n1 == 3 && TString(token1[0]) == TString("BkgRate")) {
	  if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  SystemName = TString("");
          for(int i=1;i<n1;i++) {
	    SystemName += TString(token1[i]);
	    if(i < n1-1) SystemName += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayerName")) {
	  if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	  LayerName = TString(token1[1]);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("ResolutionModel")) {
	  ResolutionModel = TString("");
          for(int i=1;i<n1;i++) {
	    ResolutionModel += TString(token1[i]);
	    if(i < n1-1) ResolutionModel += TString(" ");
          }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("EfficiencyModel")) {
	  EfficiencyModel = TString("");
          for(int i=1;i<n1;i++) {
	    EfficiencyModel += TString(token1[i]);
	    if(i < n1-1) EfficiencyModel += TString(" ");
          }
	}
        else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	  if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetDistanceUnit(token1[4]);
	  Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	  if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	  double unit = global->GetAngleUnit(token1[4]);
	  RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
        }
        else if(n1 == 2 && TString(token1[0]) == TString("IsSensitive")) {
	  IsSensitive = global->SetBoolFromString(TString(token1[1])); 
	}
        counter++;
      }
      
      if(!IsEnd) {
        cout << endl;
        cout << "Wrong specification of GeoPetal parameters inside the BeginGeoPetal and EndGeoPetal block. Parameters have to be specified as follows:" << endl;
	cout << "BeginGeoPetal" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  Position        x       y       z    units (Mandatory)" << endl;
	cout << "  RotAngles    alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
	cout << "  Thickness       val units (Mandatory)" << endl;
	cout << "  Material     string (Mandatory)" << endl;
	cout << "  XOX0            val  (Optional)" << endl;
	cout << "  bottomWidth     val units (Mandatory)" << endl;
	cout << "  topWidth        val units (Mandatory)" << endl;
	cout << "  Height          val units (Mandatory)" << endl;
	cout << "  IsSensitive     bool (Mandatory)" << endl;
	cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	cout << "  SystemName      string (Optional)" << endl;
	cout << "  LayerName       string (Optional)" << endl;
	cout << "  ResolutionModel string (Optional)" << endl;
	cout << "  EfficiencyModel string (Optional)" << endl;
	cout << "EndGeoPetal" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
      
        assert(false);
      }
      continue;
    }
    else if(n == 1 && TString(token[0]).BeginsWith("Begin") && TString(token[0]).Contains("ResolutionModel")) {
      GResolutionModel* aResolModel = NULL;
      ReadResolutionModel(TString(token[0]),&fin,aResolModel,datacard,line_number);
      if(aResolModel != NULL) aGeometry->PushResolutionModelIntoGeometry(aResolModel);      
    }
    else if(n == 1 && TString(token[0]) == TString("BeginTrackCuts") && aGeometry->GetNTrackCuts() != 1) {
      bool IsEnd = false;
      int counter = 0;
      int Max     = 500;

      TrackCuts_t ATrackCut;
      ATrackCut.SystemName = Dummy_name;
      
      while(!IsEnd && counter <= Max && !fin.eof()) {
	char buf1[MAX_CHARS_PER_LINE];
        fin.getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
	
	token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
	}
	
	counter++;
	
	if(n1 == 1 && TString(token1[0]) == TString("EndTrackCuts")) {
	  
	  if(ATrackCut.SystemName == Dummy_name) {
	    cout << endl;
	    cout << "Track cut system name not specified. Check your inputs. Exiting now!!!" << endl;
	    cout << endl;
	    assert(false);
	  }
	  
	  IsEnd = true;
	  
	  if(ATrackCut.CutList.size() > 0) aGeometry->PushTrackCutIntoGeometry(ATrackCut);
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("SystemName")) {
	  ATrackCut.SystemName = TString("");
	  for(int i=1;i<n1;i++) {
	    ATrackCut.SystemName += TString(token1[i]);
	    if(i < n1-1) ATrackCut.SystemName += TString(" ");
          }
	}
        else if(n1 == 1 && TString(token1[0]) == TString("BeginCut")) {
          bool IsEnd2 = false;
	  int counter2 = 0;
          int Max2     = 5;
	  
	  Cut_t ACut;
	  ACut.pMin     =    0.0*global->GetMomentumUnit("GeV/c");
          ACut.pMax     = 1000.0*global->GetMomentumUnit("GeV/c");
	  ACut.phiMin   = -360.0*global->GetAngleUnit("deg");
	  ACut.phiMax   = +360.0*global->GetAngleUnit("deg");
	  ACut.thetaMin =    0.0*global->GetAngleUnit("deg");
	  ACut.thetaMax =  180.0*global->GetAngleUnit("deg");
	  ACut.NhitsMin = -999.0;
	  ACut.NhitsMax = -999.0;
	  
	  while(!IsEnd2 && counter2 <= Max2 && !fin.eof()) {
	    char buf2[MAX_CHARS_PER_LINE];
            fin.getline(buf2, MAX_CHARS_PER_LINE);
	    line_number++;
            int n2 = 0;
            const char* token2[MAX_TOKENS_PER_LINE] = {};
	
	    token2[0] = strtok(buf2, DELIMITER); // first token
            if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
    
            for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
              token2[n2] = strtok(0, DELIMITER);
	      if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
	    }
	    
	    counter++;
	    counter2++;
	
	    if(n2 == 1 && TString(token2[0]) == TString("EndCut")) {
	      IsEnd2 = true;
	      
	      if(ACut.pMin > ACut.pMax) {
		cout << endl;
		cout << "Track cut pRange has lower limit " << ACut.pMin/global->GetMomentumUnit("GeV/c") << " GeV/c higher than upper limit " 
		     << ACut.pMax/global->GetMomentumUnit("GeV/c") << " GeV/c. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
		
	      }
	      if(ACut.pMin < 0.0) {
		cout << endl;
		cout << "Track cut pRange lower limit " << ACut.pMin/global->GetMomentumUnit("GeV/c") << " GeV/c has negative value. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
	      if(ACut.pMax < 0.0) {
		cout << endl;
		cout << "Track cut pRange upper limit " << ACut.pMax/global->GetMomentumUnit("GeV/c") << " GeV/c has negative value. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
          
              if(ACut.phiMin > ACut.phiMax) {
		cout << endl;
		cout << "Track cut phiRange has lower limit " << ACut.phiMin/global->GetAngleUnit("deg") << " deg higher than upper limit " 
		     << ACut.phiMax/global->GetAngleUnit("deg") << " deg. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
	      if(ACut.thetaMin > ACut.thetaMax) {
		cout << endl;
		cout << "Track cut thetaRange has lower limit " << ACut.thetaMin/global->GetAngleUnit("deg") << " deg higher than upper limit " 
		     << ACut.thetaMax/global->GetAngleUnit("deg") << " deg. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
	      
	      if(ACut.NhitsMin > ACut.NhitsMax) {
		cout << endl;
		cout << "Track cut NhitsRange has lower limit " << ACut.NhitsMin << " higher than upper limit " << ACut.NhitsMax << ". Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
	      if(ACut.NhitsMin < 0.0) {
		cout << endl;
		cout << "Track cut NhitsRange lower limit " << ACut.NhitsMin << " has negative value. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
	      if(ACut.NhitsMax < 0.0) {
		cout << endl;
		cout << "Track cut NhitsRange upper limit " << ACut.NhitsMax << " has negative value. Check your inputs. Exiting now!!!" << endl;
		cout << endl;
		assert(false);
	      }
	  
	      ATrackCut.CutList.push_back(ACut);
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("pRange")) {
	      if(!global->IsMomentumUnit(token2[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	      ACut.pMin = atof(token2[1])*global->GetMomentumUnit(TString(token2[3]));
	      ACut.pMax = atof(token2[2])*global->GetMomentumUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("phiRange")) {
	      if(!global->IsAngleUnit(token2[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	      ACut.phiMin = atof(token2[1])*global->GetAngleUnit(TString(token2[3]));
	      ACut.phiMax = atof(token2[2])*global->GetAngleUnit(TString(token2[3]));
	    }
	    else if(n2 == 4 && TString(token2[0]) == TString("thetaRange")) {
	      if(!global->IsAngleUnit(token2[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	      ACut.thetaMin = atof(token2[1])*global->GetAngleUnit(TString(token2[3]));
	      ACut.thetaMax = atof(token2[2])*global->GetAngleUnit(TString(token2[3]));
	    }
	    else if(n2 == 3 && TString(token2[0]) == TString("NhitsRange")) {
	      ACut.NhitsMin = atof(token2[1]);
	      ACut.NhitsMax = atof(token2[2]);
	    }
	  }
	  
	  if(!IsEnd2) {
	    cout << endl;
	    cout << "Wrong specification of Cut parameters inside the BeginCut and EndCut block. Parameters have to be specified as follows:" << endl;
	    cout << "BeginCut" << endl;
	    cout << "  pRange       min  max  units" << endl;
	    cout << "  phiRange     min  max  units" << endl;
	    cout << "  thetaRange   min  max  units" << endl;
	    cout << "  NhitsRange   min  max"        << endl;
	    cout << "EndCut"   << endl;
	    cout << "Check you inputs. Exiting now!!!" << endl;
	    cout << endl;
	    assert(false);
	  }
	  
        }
	
      }
      
      if(!IsEnd) {
	cout << endl;
	cout << "Wrong specification of Track cut parameters inside the BeginTrackCuts and EndTrackCuts block. Parameters have to be specified as follows:" << endl;
	cout << "BeginTrackCuts" << endl;
	cout << "  SystemName  string" << endl;
	cout << "  BeginCut" << endl;
	cout << "    pRange       min  max  units" << endl;
	cout << "    phiRange     min  max  units" << endl;
	cout << "    thetaRange   min  max  units" << endl;
	cout << "    NhitsRange   min  max"        << endl;
	cout << "  EndCut"   << endl;
	cout << "  ...  " << endl;
	cout << "  BeginCut" << endl;
	cout << "    pRange       min  max  units" << endl;
	cout << "    phiRange     min  max  units" << endl;
	cout << "    thetaRange   min  max  units" << endl;
	cout << "    NhitsRange   min  max"        << endl;
	cout << "  EndCut"   << endl;
	cout << "EndTrackCuts"   << endl;
	cout << "Check you inputs. Exiting now!!!" << endl;
	cout << endl;
	assert(false);
      }
    }
    else if(n == 1 && TString(token[0]).BeginsWith("Begin") && TString(token[0]).Contains("Ladder")) {
      ReadLadderBlock(TString(token[0]),&fin,aGeometry,datacard,line_number);
    }

  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::ReadBField(TString FieldName, TString  BlockName,ifstream* fin, GBField* &aBField, const char* datacard, int &line_number)
{
  
  if(aBField != NULL) return;
  
  if(BlockName == TString("BeginConstantBfield")) {
    bool IsEnd  = false;
    int counter = 0;
    int Max     = 2;
    
    //Parameters of a constant B-field
    double    Magnitude = -1;
    TVector3  Direction(1.0,0.0,0.0);
      
    while(!IsEnd && counter <= Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      if(n1 == 1 && TString(token1[0]) == TString("EndConstantBfield")) {
	if(Magnitude < 0) continue;
	
	Direction = Direction.Unit();
	
	aBField = new GBFieldConstant(FieldName.Data(),Direction,Magnitude,global);
	
	IsEnd = true;
      }
      else if(n1 == 3 && TString(token1[0]) == TString("Magnitude"))  {
	if(!global->IsBfieldUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	Magnitude = atof(token1[1])*global->GetBfieldUnit(TString(token1[2]));
      }
      else if(n1 == 4 && TString(token1[0]) == TString("Direction"))  Direction = TVector3(atof(token1[1]),atof(token1[2]),atof(token1[3]));
      counter++;
    }
      
    if(!IsEnd) {
      cout << endl;
      cout << "Wrong specification of Magnetic Field parameters inside the BeginConstantBfield and EndConstantBfield block. Parameters have to be specified as follows:" << endl;
      cout << "BeginConstantBfield" << endl;
      cout << "  Magnitude    val  units //val >= 0" << endl;
      cout << "  Direction    x  y  z" << endl;
      cout << "EndConstantBfield"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }
  else if(BlockName == TString("BeginMultipleStepsBfield")) {
    bool IsEnd  = false;
    int counter = 0;
    int Max     = 10000;
    
    std::vector<GGeoObject*>  VolumeList;     //List of Volumes
    std::vector<TVector3>     Bfield_inList;  //List of inside B-fields
    TVector3     Bfield_out  = Dummy_vector;  //B-field outside volume
    
    VolumeList.clear();
    Bfield_inList.clear();
      
    while(!IsEnd && counter <= Max && !fin->eof()) {
      char buf[MAX_CHARS_PER_LINE];
      fin->getline(buf, MAX_CHARS_PER_LINE);
      line_number++;
      int n = 0;
      const char* token[MAX_TOKENS_PER_LINE] = {};
      
      token[0] = strtok(buf, DELIMITER); // first token
      if(!token[0] || TString(token[0]).BeginsWith(COMMENT)) continue;
    
      for(n=1; n<MAX_TOKENS_PER_LINE; n++) {
        token[n] = strtok(0, DELIMITER);
	if(!token[n] || TString(token[n]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      if(n == 1 && TString(token[0]) == TString("EndMultipleStepsBfield")) {
	bool GoodBfield = true;
	
	if(Bfield_inList.size() == 0) {
	  cout << endl;
	  cout << "ERROR inside MultipleStepsBfield block: both B-fields inside and outside volume not set!!!" << endl;
	  cout << endl;
	  GoodBfield = false;
	}
	else if(VolumeList.size() == 0) {
	  cout << endl;
	  cout << "ERROR inside MultipleStepsBfield block: volume not set!!!" << endl;
	  cout << endl;
	  GoodBfield = false;
	}
	
	if(Bfield_out == Dummy_vector) {
	  cout << endl;
	  cout << "WARNING inside MultipleStepsBfield block: B-field outside volume not set. Setting it to zero." << endl;
	  cout << endl;
	  Bfield_out = TVector3(0.0,0.0,0.0);
	}
	
	if(GoodBfield) {
	  aBField = new GBFieldMultipleSteps(FieldName.Data(),Bfield_inList,VolumeList,Bfield_out,global);
	  IsEnd = true;
	  for(int i=0;i<int(VolumeList.size());i++) delete VolumeList[i];
	}
	
      }
      else if(n == 1 && TString(token[0]) == TString("BeginInsideBfield"))  {
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 4;
    
	double    Magnitude = Dummy_value;
	TVector3  Direction = Dummy_vector;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
                    
          if(n1 == 1 && TString(token1[0]) == TString("EndInsideBfield")) {
	    bool Good = true;
	    
	    if(Magnitude == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in InsideBfield block: B-field magnitude not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Magnitude < 0.0) {
	      cout << endl;
	      cout << "  ERROR in InsideBfield block: B-field magnitude not negative value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Direction == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in InsideBfield block: B-field direction not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Direction.Mag() < 1.0e-6) {
	      cout << endl;
	      cout << "  ERROR in InsideBfield block: B-field direction vector set to zero!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(Good) {
	      Direction = Direction.Unit();
	      Bfield_inList.push_back(Magnitude*Direction);
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Magnitude"))  {
	    if(!global->IsBfieldUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Magnitude = atof(token1[1])*global->GetBfieldUnit(TString(token1[2]));
	  }
	  else if(n1 == 4 && TString(token1[0]) == TString("Direction"))  Direction = TVector3(atof(token1[1]),atof(token1[2]),atof(token1[3]));
	  
	  counter1++;
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of inside magnetic field parameters inside the BeginInsideBfield and EndInsideBfield block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginInsideBfield" << endl;
	  cout << "  Magnitude    val  units //val >= 0" << endl;
	  cout << "  Direction    x  y  z" << endl;
	  cout << "EndInsideBfield"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
	
      }
      else if(n == 1 && TString(token[0]) == TString("BeginOutsideBfield"))  {
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 4;
    
	double    Magnitude = Dummy_value;
	TVector3  Direction = Dummy_vector;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          if(n1 == 1 && TString(token1[0]) == TString("EndOutsideBfield")) {
	    bool Good = true;
	    
	    if(Magnitude == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in OutsideBfield block: B-field magnitude not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Magnitude < 0.0) {
	      cout << endl;
	      cout << "  ERROR in OutsideBfield block: B-field magnitude not negative value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Direction == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in OutsideBfield block: B-field direction not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Direction.Mag() < 1.0e-6) {
	      cout << endl;
	      cout << "  ERROR in OutsideBfield block: B-field direction vector set to zero!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(Good) {
	      Direction  = Direction.Unit();
	      Bfield_out = Magnitude*Direction;
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Magnitude"))  {
	    if(!global->IsBfieldUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Magnitude = atof(token1[1])*global->GetBfieldUnit(TString(token1[2]));
	  }
	  else if(n1 == 4 && TString(token1[0]) == TString("Direction"))  Direction = TVector3(atof(token1[1]),atof(token1[2]),atof(token1[3]));
	  
	  counter1++;
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of outside magnetic field parameters inside the BeginInsideBfield and EndInsideBfield block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginOutsideBfield" << endl;
	  cout << "  Magnitude    val  units //val >= 0" << endl;
	  cout << "  Direction    x  y  z" << endl;
	  cout << "EndOutsideBfield"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
      }
      else if(n == 1 && TString(token[0]) == TString("BeginBoxVolume"))  {
	//Box volume
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
    
	TVector3 Position  = Dummy_vector;
	TVector3 RotAngles = Dummy_vector;
	double   widthX    = Dummy_value;
	double   widthY    = Dummy_value;
	double   widthZ    = Dummy_value;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          if(n1 == 1 && TString(token1[0]) == TString("EndBoxVolume")) {
	    bool Good = true;
	    
	    if(Position == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: position not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(widthX == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: box widthX not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(widthX <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: box widthX set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(widthY == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: box widthY not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(widthY <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: box widthY set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(widthZ == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: box widthZ not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(widthZ <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield BoxVolume block: box widthZ set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(RotAngles == Dummy_vector) {
	      cout << endl;
	      cout << "  WARNING in StepBfield BoxVolume block: rotation angles not specified. Setting them to (0,0,0)."<< endl;
	      cout << endl;
	      RotAngles = TVector3(0,0,0);
	    }
	    
	    if(Good) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(RotAngles,Rot);
	      
	      TString Name     = FieldName + TString(" Volume ") + long(VolumeList.size()+1);
	      TString Material("Vacuum");
	      VolumeList.push_back(new GGeoPlane(Name,0,0,
						 Position,Rot,
					         widthZ,Material,
					         false,
					         TVector2(widthX,widthY),
					         TVector2(0,0),TVector2(0,0),
						 global));
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("widthX")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthX = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("widthY")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthY = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("widthZ")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthZ = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token1[4]);
	    Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	    if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(token1[4]);
	    RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
          
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of step B-field Box volume parameters inside the BeginBoxVolume and EndBoxVolume block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginBoxVolume" << endl;
	  cout << "  Position        x       y       z    units" << endl;
	  cout << "  RotAngles    alphaX  alphaY  alphaZ  units" << endl;
	  cout << "  widthX         val  units // val > 0" << endl;
	  cout << "  widthY         val  units // val > 0" << endl;
	  cout << "  widthZ         val  units // val > 0" << endl;
	  cout << "EndBoxVolume"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
      }
      else if(n == 1 && TString(token[0]) == TString("BeginCylinderVolume"))  {
	//Cylinder volume
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
    
	TVector3 Position  = Dummy_vector;
	TVector3 RotAngles = Dummy_vector;
	double   Radius    = Dummy_value;
	double   Length    = Dummy_value;
	
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          if(n1 == 1 && TString(token1[0]) == TString("EndCylinderVolume")) {
	    bool Good = true;
	    
	    if(Position == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in StepBfield CylinderVolume block: position not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield CylinderVolume block: cylinder Length not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield CylinderVolume block: cylinder Length set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield CylinderVolume block: cylinder Radius not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield CylinderVolume block: cylinder Radius set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(RotAngles == Dummy_vector) {
	      cout << endl;
	      cout << "  WARNING in StepBfield CylinderVolume block: rotation angles not specified. Setting them to (0,0,0)."<< endl;
	      cout << endl;
	      RotAngles = TVector3(0,0,0);
	    }
	    
	    if(Good) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(RotAngles,Rot);
	      
	      TString Name     = FieldName + TString(" Volume ") + long(VolumeList.size()+1);
	      TString Material("Vacuum");
	      VolumeList.push_back(new GGeoDisk(Name,0,0,
						Position,Rot,
					        Length,Material,
			                        false,
			                        0.0,Radius,
			                        TVector2(0,0),
				                global));
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Length = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Radius")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Radius = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token1[4]);
	    Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	    if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(token1[4]);
	    RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
          
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of step B-field Box volume parameters inside the BeginCylinderVolume and EndCylinderVolume block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginCylinderVolume" << endl;
	  cout << "  Position        x       y       z    units" << endl;
	  cout << "  RotAngles    alphaX  alphaY  alphaZ  units" << endl;
	  cout << "  Length         val  units // val > 0" << endl;
	  cout << "  Radius         val  units // val > 0" << endl;
	  cout << "EndCylinderVolume"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
      }
      else if(n == 1 && TString(token[0]) == TString("BeginDiskSectionVolume"))  {
	//Cylinder volume
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
    
	TVector3 Position  = Dummy_vector;
	TVector3 RotAngles = Dummy_vector;
	double   Rin       = Dummy_value;
	double   Rout      = Dummy_value;
	double   Length    = Dummy_value;
	double   DeltaPhi  = Dummy_value;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          if(n1 == 1 && TString(token1[0]) == TString("EndDiskSectionVolume")) {
	    bool Good = true;
	    
	    if(Position == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: position not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: disk Length not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: disk Length set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Rin == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: disk Rin not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Rin <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: disk Rin set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Rout == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: disk Rout not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Rout <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: disk Rout set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(DeltaPhi == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: cylinder DeltaPhi not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(DeltaPhi <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield DiskSectionVolume block: cylinder DeltaPhi set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(RotAngles == Dummy_vector) {
	      cout << endl;
	      cout << "  WARNING in StepBfield DiskSectionVolume block: rotation angles not specified. Setting them to (0,0,0)."<< endl;
	      cout << endl;
	      RotAngles = TVector3(0,0,0);
	    }
	    
	    if(Good) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(RotAngles,Rot);
	      
	      TString Name     = FieldName + TString(" Volume ") + long(VolumeList.size()+1);
	      TString Material("Vacuum");
	      VolumeList.push_back(new GGeoDiskSection(Name,0,0,
						       Position,Rot,
						       Length,Material,
			                               false,
			                               Rin,Rout,
						       DeltaPhi,
			                               TVector2(0,0),TVector2(0,0),
				                       global));
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Length = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Rin")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Rin = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Rout")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Rout = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("DeltaPhi")) {
	    if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    DeltaPhi = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token1[4]);
	    Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	    if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(token1[4]);
	    RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
          
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of step B-field Box volume parameters inside the BeginDiskSectionVolume and EndDiskSectionVolume block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginDiskSectionVolume" << endl;
	  cout << "  Position        x       y       z    units" << endl;
	  cout << "  RotAngles    alphaX  alphaY  alphaZ  units" << endl;
	  cout << "  Length         val  units // val > 0" << endl;
	  cout << "  Rin            val  units // val > 0" << endl;
          cout << "  Rout           val  unit  // val > 0" << endl;
	  cout << "  DeltaPhi       val  units // val > 0" << endl;
	  cout << "EndDiskSectionVolume"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
      }
      else if(n == 1 && TString(token[0]) == TString("BeginConeVolume"))  {
	//Cylinder volume
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
    
	TVector3 Position  = Dummy_vector;
	TVector3 RotAngles = Dummy_vector;
	double   Radius1   = Dummy_value;
	double   Radius2   = Dummy_value;
	double   Length    = Dummy_value;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          if(n1 == 1 && TString(token1[0]) == TString("EndConeVolume")) {
	    bool Good = true;
	    
	    if(Position == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: position not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: cone Length not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: cone Length set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius1 == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: cone Radius1 not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius1 <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: cone Radius1 set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius2 == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: cone Radius2 not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius2 <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeVolume block: cone Radius2 set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(RotAngles == Dummy_vector) {
	      cout << endl;
	      cout << "  WARNING in StepBfield ConeVolume block: rotation angles not specified. Setting them to (0,0,0)."<< endl;
	      cout << endl;
	      RotAngles = TVector3(0,0,0);
	    }
	    
	    if(Good) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(RotAngles,Rot);
	      
	      TString Name     = FieldName + TString(" Volume ") + long(VolumeList.size()+1);
	      TString Material("Vacuum");
	      VolumeList.push_back(new GGeoCone(Name,0,0,
						Position,Rot,
					        0.5*Radius1,0.5*Radius2,
			                        Material,
			                        false,
			                        Length,
			                        0.5*Radius1,0.5*Radius2,
			                        TVector2(0,0),
				                global));
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Length = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Radius1")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Radius1 = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Radius2")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Radius2 = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token1[4]);
	    Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	    if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(token1[4]);
	    RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
          
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of step B-field Box volume parameters inside the BeginConeVolume and EndConeVolume block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginConeVolume" << endl;
	  cout << "  Position        x       y       z    units" << endl;
	  cout << "  RotAngles    alphaX  alphaY  alphaZ  units" << endl;
	  cout << "  Length         val  units // val > 0" << endl;
	  cout << "  Radius1        val  units // val > 0" << endl;
	  cout << "  Radius2        val  units // val > 0" << endl;
	  cout << "EndConeVolume"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
      }
      else if(n == 1 && TString(token[0]) == TString("BeginConeSectionVolume"))  {
	//Cylinder volume
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
    
	TVector3 Position  = Dummy_vector;
	TVector3 RotAngles = Dummy_vector;
	double   Radius1   = Dummy_value;
	double   Radius2   = Dummy_value;
	double   Length    = Dummy_value;
	double   DeltaPhi  = Dummy_value;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          if(n1 == 1 && TString(token1[0]) == TString("EndConeSectionVolume")) {
	    bool Good = true;
	    
	    if(Position == Dummy_vector) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: position not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone Length not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Length <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone Length set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius1 == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone Radius1 not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius1 <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone Radius1 set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius2 == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone Radius2 not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(Radius2 <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone Radius2 set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(DeltaPhi == Dummy_value) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone DeltaPhi not set!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(DeltaPhi <= 0.0) {
	      cout << endl;
	      cout << "  ERROR in StepBfield ConeSectionVolume block: cone DeltaPhi set to <= 0 value!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(RotAngles == Dummy_vector) {
	      cout << endl;
	      cout << "  WARNING in StepBfield ConeSectionVolume block: rotation angles not specified. Setting them to (0,0,0)."<< endl;
	      cout << endl;
	      RotAngles = TVector3(0,0,0);
	    }
	    
	    if(Good) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(RotAngles,Rot);
	      
	      TString Name     = FieldName + TString(" Volume ") + long(VolumeList.size()+1);
	      TString Material("Vacuum");
	      VolumeList.push_back(new GGeoConeSection(Name,0,0,
						       Position,Rot,
					               0.5*Radius1,0.5*Radius2,
			                               Material,
			                               false,
			                               Length,
			                               0.5*Radius1,0.5*Radius2,
						       DeltaPhi,
			                               TVector2(0,0),TVector2(0,0),
				                       global));
	      
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Length = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Radius1")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Radius1 = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("Radius2")) {
	    if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Radius2 = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
	  }
	  else if(n1 == 3 && TString(token1[0]) == TString("DeltaPhi")) {
	    if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    DeltaPhi = atof(token1[1])*global->GetAngleUnit(TString(token1[2]));
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token1[4]);
	    Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
	  else if(n1 == 5 && TString(token1[0]) == TString("RotAngles")) {
	    if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(token1[4]);
	    RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
	  }
          
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of step B-field Box volume parameters inside the BeginConeSectionVolume and EndConeSectionVolume block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginConeSectionVolume" << endl;
	  cout << "  Position        x       y       z    units" << endl;
	  cout << "  RotAngles    alphaX  alphaY  alphaZ  units" << endl;
	  cout << "  Length         val  units // val > 0" << endl;
	  cout << "  Radius1        val  units // val > 0" << endl;
	  cout << "  Radius2        val  units // val > 0" << endl;
	  cout << "  DeltaPhi       val  units // val > 0" << endl;
	  cout << "EndConeSectionVolume"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
      }
      
      counter++;
    }
      
    if(!IsEnd) {
      cout << endl;
      cout << "Wrong specification of Step magnetic field parameters inside the BeginMultipleStepsBfield and EndMultipleStepsBfield block. Parameters have to be specified as follows:" << endl;
      cout << "BeginMultipleStepsBfield" << endl;
      cout << "  BeginInsideBfield" << endl;
      cout << "    Magnitude    val  units //val >= 0" << endl;
      cout << "    Direction    x  y  z" << endl;
      cout << "  EndInsideBfield"   << endl;
      cout << "  // Then the volume, which is either a Box " << endl;
      cout << "  BeginBoxVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    widthX         val  units // val > 0" << endl;
      cout << "    widthY         val  units // val > 0" << endl;
      cout << "    widthZ         val  units // val > 0" << endl;
      cout << "  EndBoxVolume"   << endl;
      cout << "  // or a Cylinder" << endl;
      cout << "  BeginCylinderVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Radius         val  units // val > 0" << endl;
      cout << "  EndCylinderVolume"   << endl;
      cout << "  // or a Cone" << endl;
      cout << "  BeginConeVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Radius1        val  units // val > 0" << endl;
      cout << "    Radius2        val  units // val > 0" << endl;
      cout << "  EndConeVolume"   << endl;
      cout << "  // or a Cylinder section" << endl;
      cout << "  BeginDiskSectionVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Rin            val  units // val > 0" << endl;
      cout << "    Rout           val  units // val > 0" << endl;
      cout << "    DeltaPhi       val  units // val > 0" << endl;
      cout << "  EndDiskSectionVolume"   << endl;
      cout << "  // or a Cone" << endl;
      cout << "  BeginConeSectionVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Radius1        val  units // val > 0" << endl;
      cout << "    Radius2        val  units // val > 0" << endl;
      cout << "    DeltaPhi       val  units // val > 0" << endl;
      cout << "  EndConeSectionVolume"   << endl;
      cout << " // ........" << endl;
      cout << "  BeginInsideBfield" << endl;
      cout << "    Magnitude    val  units //val >= 0" << endl;
      cout << "    Direction    x  y  z" << endl;
      cout << "  EndInsideBfield"   << endl;
      cout << "  // Then the volume, which is either a Box " << endl;
      cout << "  BeginBoxVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    widthX         val  units // val > 0" << endl;
      cout << "    widthY         val  units // val > 0" << endl;
      cout << "    widthZ         val  units // val > 0" << endl;
      cout << "  EndBoxVolume"   << endl;
      cout << "  // or a Cylinder" << endl;
      cout << "  BeginCylinderVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Radius         val  units // val > 0" << endl;
      cout << "  EndCylinderVolume"   << endl;
      cout << "  // or a Cone" << endl;
      cout << "  BeginConeVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Radius1        val  units // val > 0" << endl;
      cout << "    Radius2        val  units // val > 0" << endl;
      cout << "  EndConeVolume"   << endl;
      cout << "  // or a Cylinder section" << endl;
      cout << "  BeginDiskSectionVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Rin            val  units // val > 0" << endl;
      cout << "    Rout           val  units // val > 0" << endl;
      cout << "    DeltaPhi       val  units // val > 0" << endl;
      cout << "  EndDiskSectionVolume"   << endl;
      cout << "  // or a Cone" << endl;
      cout << "  BeginConeSectionVolume" << endl;
      cout << "    Position        x       y       z    units" << endl;
      cout << "    RotAngles    alphaX  alphaY  alphaZ  units" << endl;
      cout << "    Length         val  units // val > 0" << endl;
      cout << "    Radius1        val  units // val > 0" << endl;
      cout << "    Radius2        val  units // val > 0" << endl;
      cout << "    DeltaPhi       val  units // val > 0" << endl;
      cout << "  EndConeSectionVolume"   << endl;
      cout << "  BeginOutsideBfield" << endl;
      cout << "    Magnitude    val  units //val >= 0" << endl;
      cout << "    Direction    x  y  z" << endl;
      cout << "  EndOutsideBfield"   << endl;
      cout << "EndMultipleStepsBfield"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }
  
  
  return;
  
}
//====================================================================
void  Guariguanchi::ReadWorldVolume(TString WorldVolName,TString  BlockName,ifstream* fin,GGeoObject* &aWorldolume, int geoID, const char* datacard, int &line_number)
{
  
  if(aWorldolume != NULL) return;

  if(BlockName == TString("BeginBoxWorldVolume")) {
    bool IsEnd  = false;
    int counter = 0;
    int Max     = 4;
    
    double    widthX = -1.0;
    double    widthY = -1.0;
    double    widthZ = -1.0;
    TVector3  Position(0.0,0.0,0.0);
      
    while(!IsEnd && counter <= Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      if(n1 == 1 && TString(token1[0]) == TString("EndBoxWorldVolume")) {
	if(widthX < 0) continue;
	if(widthY < 0) continue;
	if(widthZ < 0) continue;
	
	TMatrixD Rot;
	global->GetGlobalRotationMatrix(TVector3(0,0,0),Rot);
	
	aWorldolume = new GGeoPlane(WorldVolName,geoID,0,
				    Position,Rot,widthZ,WorldVolMaterial,false,
				    TVector2(widthX,widthY),
				    TVector2(0,0),TVector2(0,0),
				    global);
	
	IsEnd = true;
      }
      else if(n1 == 3 && TString(token1[0]) == TString("widthX")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	widthX = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
      }
      else if(n1 == 3 && TString(token1[0]) == TString("widthY")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	widthY = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
      }
      else if(n1 == 3 && TString(token1[0]) == TString("widthZ")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	widthZ = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
      }
      else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit = global->GetDistanceUnit(token1[4]);
	Position  = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }      
      counter++;
    }
      
    if(!IsEnd) {
      cout << endl;
      cout << "Wrong specification of World Volume parameters inside the BeginBoxWorldVolume and EndBoxWorldVolume block. Parameters have to be specified as follows:" << endl;
      cout << "BeginBoxWorldVolume" << endl;
      cout << "  Position   x   y   z  units" << endl;
      cout << "  widthX     val  units  //val >= 0" << endl;
      cout << "  widthY     val  units  //val >= 0" << endl;
      cout << "  widthZ     val  units  //val >= 0" << endl;
      cout << "EndBoxWorldVolume"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }
  else if(BlockName == TString("BeginCylinderWorldVolume")) {
    bool IsEnd  = false;
    int counter = 0;
    int Max     = 3;
    
    double    Radius = -1.0;
    double    Length = -1.0;
    TVector3  Position(0.0,0.0,0.0);
      
    while(!IsEnd && counter <= Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      if(n1 == 1 && TString(token1[0]) == TString("EndCylinderWorldVolume")) {
	if(Radius < 0) continue;
	if(Length < 0) continue;
	
	TMatrixD Rot;
	global->GetGlobalRotationMatrix(TVector3(0,0,0),Rot);
	
	aWorldolume = new GGeoDisk(WorldVolName,geoID,0,
				   Position,Rot,Length,WorldVolMaterial,false,
			           0.0,Radius,
			           TVector2(0,0),
				   global);

	IsEnd = true;
      }
      else if(n1 == 3 && TString(token1[0]) == TString("Radius")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	Radius = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
      }
      else if(n1 == 3 && TString(token1[0]) == TString("Length")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	Length = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
      }
      else if(n1 == 5 && TString(token1[0]) == TString("Position")) {
	if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit = global->GetDistanceUnit(token1[4]);
	Position  = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }      
      counter++;
    }
      
    if(!IsEnd) {
      cout << endl;
      cout << "Wrong specification of World Volume parameters inside the BeginCylinderWorldVolume and EndCylinderWorldVolume block. Parameters have to be specified as follows:" << endl;
      cout << "BeginCylinderWorldVolume" << endl;
      cout << "  Position   x   y   z  units" << endl;
      cout << "  Radius     val  units  //val >= 0" << endl;
      cout << "  Length     val  units  //val >= 0" << endl;
      cout << "EndCylinderWorldVolume"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }

  return;
  
}
//====================================================================
void  Guariguanchi::ReadResolutionModel(TString  BlockName,ifstream* fin,GResolutionModel* &aResolModel, const char* datacard, int &line_number)
{
  
  if(aResolModel != NULL) return;

  if(BlockName == TString("BeginTPCResolutionModel")) {
    bool IsEnd  = false;
    int counter = 0;
    int Max     = 6;
    
    //Parameters of a constant B-field
    TString   Name = Dummy_name;
    double    A_U  = Dummy_value;
    double    B_U  = Dummy_value;
    double    C_U  = Dummy_value;
    double    A_V  = Dummy_value;
    double    B_V  = Dummy_value;
      
    while(!IsEnd && counter <= Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      if(n1 == 1 && TString(token1[0]) == TString("EndTPCResolutionModel")) {
	bool GoodResolutionModel = true;
	
	if(Name == Dummy_name) {
	  cout << "ERROR in TPC ResolutionModel Block: Name not specified!!!" << endl;
	  GoodResolutionModel = false;
	}
	else if(A_U == Dummy_value) {
	  cout << "ERROR in TPC ResolutionModel Block: A_U parameter not specified!!!" << endl;
	  GoodResolutionModel = false;
	}
	else if(B_U == Dummy_value) {
	  cout << "ERROR in TPC ResolutionModel Block: B_U parameter not specified!!!" << endl;
	  GoodResolutionModel = false;
	}
	else if(C_U == Dummy_value) {
	  cout << "ERROR in TPC ResolutionModel Block: C_U parameter not specified!!!" << endl;
	  GoodResolutionModel = false;
	}
	else if(A_V == Dummy_value) {
	  cout << "ERROR in TPC ResolutionModel Block: A_V parameter not specified!!!" << endl;
	  GoodResolutionModel = false;
	}
	else if(B_V == Dummy_value) {
	  cout << "ERROR in TPC ResolutionModel Block: A_B parameter not specified!!!" << endl;
	  GoodResolutionModel = false;
	}
	
	if(GoodResolutionModel) {
	  aResolModel = new GResolutionModelTPC(Name,
						A_U,B_U,C_U,
					        A_V,B_V,
					        global);
	  
	  IsEnd = true;
	}
	
      }
      else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	Name = TString("");
        for(int i=1;i<n1;i++) {
	  Name += TString(token1[i]);
	  if(i < n1-1) Name += TString(" ");
        }
      }
      else if(n1 == 3 && TString(token1[0]) == TString("A_U"))   A_U  = atof(token1[1])*global->GetUnit(TString(token1[2]));
      else if(n1 == 3 && TString(token1[0]) == TString("B_U"))   B_U  = atof(token1[1])*global->GetUnit(TString(token1[2]));
      else if(n1 == 3 && TString(token1[0]) == TString("C_U"))   C_U  = atof(token1[1])*global->GetUnit(TString(token1[2]));
      else if(n1 == 3 && TString(token1[0]) == TString("A_V"))   A_V  = atof(token1[1])*global->GetUnit(TString(token1[2]));
      else if(n1 == 3 && TString(token1[0]) == TString("B_V"))   B_V  = atof(token1[1])*global->GetUnit(TString(token1[2]));
      counter++;
    }
      
    if(!IsEnd) {
      cout << endl;
      cout << "Wrong specification of Resolution Model parameters inside the BeginTPCResolutionModel and EndTPCResolutionModel block. Parameters have to be specified as follows:" << endl;
      cout << "BeginTPCResolutionModel" << endl;
      cout << "  Name    string (Mandatory)" << endl;
      cout << "  A_U       val  units // val >= 0 (Mandatory)" << endl;
      cout << "  B_U       val  units // val >= 0 (Mandatory)" << endl;
      cout << "  C_U       val  units // val >= 0 (Mandatory)" << endl;
      cout << "  A_V       val  units // val >= 0 (Mandatory)" << endl;
      cout << "  B_V       val  units // val >= 0 (Mandatory)" << endl;
      cout << "EndTPCResolutionModel"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::ReadLadderBlock(TString  BlockName,ifstream* fin,GGeometry* aGeometry, const char* datacard, int &line_number)
{
  
  if(aGeometry == NULL) return;

  int  Ladder_Max = 10000;
  
  if(BlockName == TString("BeginLadderPlane")) {
    bool IsLadderEnd    = false;
    int  Ladder_counter = 0;
    
    //Parameters
    int       NPlanes         = 0;
    TString   LadderName      = Dummy_name;
    TVector3  LadderPosition  = Dummy_vector;
    TVector3  LadderRotAngles = Dummy_vector;
    
    while(!IsLadderEnd && Ladder_counter <= Ladder_Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      Ladder_counter++;
      
      if(n1 == 1 && TString(token1[0]) == TString("EndLadderPlane")) {
	if(LadderRotAngles == Dummy_vector) {
	  cout << "WARNING in LadderPlane Block: LadderRotAngles not specified. Setting all to zero." << endl;
	  LadderRotAngles = TVector3(0,0,0);
	}
	
	if(LadderName == Dummy_name) {
	  cout << "ERROR in LadderPlane Block: Name not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderPosition == Dummy_vector) {
	  cout << "ERROR in LadderPlane Block: LadderPosition parameter not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(NPlanes == 0) {
	  cout << "ERROR in LadderPlane Block: Nplanes is zero!!!" << endl;
	  IsLadderEnd = false;
	}
	else IsLadderEnd = true;
      }
      else if(n1 >= 2 && TString(token1[0]) == TString("LadderName")) {
	LadderName = TString("");
        for(int i=1;i<n1;i++) {
	  LadderName += TString(token1[i]);
	  if(i < n1-1) LadderName += TString(" ");
        }
      }
      else if(n1 == 5 && TString(token1[0]) == TString("LadderPosition")) {
	if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit    = global->GetDistanceUnit(token1[4]);
	LadderPosition = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 5 && TString(token1[0]) == TString("LadderRotAngles")) {
	if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit     = global->GetAngleUnit(token1[4]);
	LadderRotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 1 && TString(token1[0]) == TString("BeginPlane")) {
	bool IsPlaneEnd    = false;
        int  Plane_counter = 0;
        int  Plane_Max     = 24;
    
        TString   Name           = Dummy_name;
        TVector3  Position       = Dummy_vector;
        double    Thickness      = Dummy_value;
	double    XOX0           = Dummy_value;
        TString   Material       = Dummy_name;
        double    widthU         = Dummy_value;
        double    widthV         = Dummy_value;
        bool      IsSensitive    = false;
        double    Resolution     = Dummy_value;
        double    ResolutionU    = Dummy_value;
        double    ResolutionV    = Dummy_value;
        double    ROtime         = Dummy_value;
        double    DetEffic       = Dummy_value;
        double    InsensFracUneg = 0.0;
        double    InsensFracUpos = 0.0;
        double    InsensFracVneg = 0.0;
        double    InsensFracVpos = 0.0;
        double    BkgRate        = 0.0;
        TString   SystemName("");
	TString   LayerName("");
        TString   ResolutionModel("");
	TString   EfficiencyModel("");

        while(!IsPlaneEnd && Plane_counter <= Plane_Max && !fin->eof()) {
          char buf2[MAX_CHARS_PER_LINE];
          fin->getline(buf2, MAX_CHARS_PER_LINE);
	  line_number++;
          int n2 = 0;
          const char* token2[MAX_TOKENS_PER_LINE] = {};
      
          token2[0] = strtok(buf2, DELIMITER); // first token
          if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
    
          for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
            token2[n2] = strtok(0, DELIMITER);
	    if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          Ladder_counter++;
	  Plane_counter++;
          
          if(n2 == 1 && TString(token2[0]) == TString("EndPlane")) {
	    bool GoodGeoObject = true;
	    
	    if(LadderName == Dummy_name) {
	      cout << "ERROR in LadderPlane Block: Name not specified!!!" << endl;
	      assert(false);
	    }
	    else if(LadderPosition == Dummy_vector) {
	      cout << "ERROR in LadderPlane Block: LadderPosition parameter not specified!!!" << endl;
	      assert(false);
	    }

	    if(Name == Dummy_name) {
	      cout << "ERROR in Plane Block inside LadderPlane: Name not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Position == Dummy_vector) {
	      cout << "ERROR in Plane Block  inside LadderPlane: Position not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Thickness == Dummy_value) {
	      cout << "ERROR in Plane Block  inside LadderPlane: Thickness not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    if(Material == Dummy_name) {
	      cout << "ERROR in Plane Block  inside LadderPlane: Material not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(widthU == Dummy_value) {
	      cout << "ERROR in Plane Block  inside LadderPlane: widthU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(widthV == Dummy_value) {
	      cout << "ERROR in Plane Block  inside LadderPlane: widthV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    
	    if(LadderRotAngles == Dummy_vector) {
	      cout << "WARNING in LadderPlane Block: LadderRotAngles not specified. Setting all to zero." << endl;
	      LadderRotAngles = TVector3(0,0,0);
	    }
	  
	    if(IsSensitive) {
	      if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	        if(Resolution == Dummy_value) {
		  cout << "ERROR in Plane Block  inside LadderPlane: Resolution not specified!!!" << endl;
		  GoodGeoObject = false;
	        }
	        else {
		  ResolutionU = Resolution;
		  ResolutionV = Resolution;
	        }
	      }
	      else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	        cout << "ERROR in Plane Block  inside LadderPlane: ResolutionU not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	        cout << "ERROR in Plane Block  inside LadderPlane: ResolutionV not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ROtime == Dummy_value) {
	        cout << "ERROR in Plane Block  inside LadderPlane: ROtime not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(DetEffic == Dummy_value) {
	        cout << "ERROR in Plane Block  inside LadderPlane: DetEffic not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	    }
	  
	    if(GoodGeoObject) {
	      NPlanes++;
	      
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(LadderRotAngles,Rot);
	      TString FullName = LadderName + TString(" ") + Name;
	      
	      TMatrixD InvRot;
	      InvRot.ResizeTo(3,3);
	      InvRot = Rot;
	      InvRot.Invert();
	      TVector3  GlobalPosition = Position;
	      global->RotateVector(InvRot,GlobalPosition);
	      GlobalPosition += LadderPosition;
	      
	      GGeoObject* aGeoObject = new GGeoPlane(FullName,aGeometry->GetGeometryID(),0,
						     GlobalPosition,Rot,Thickness,Material,
					             IsSensitive,
					             TVector2(widthU,widthV),
						     TVector2(InsensFracUneg,InsensFracUpos),
						     TVector2(InsensFracVneg,InsensFracVpos),
						     global,
					             ResolutionU,ResolutionV,DetEffic,ROtime,
					             BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      aGeoObject->SetLadderType(TString("LadderPlane"));
	      if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	      aGeometry->PushGeoElement(aGeoObject);
	      
	      IsPlaneEnd = true;
	    }

          }
          else if(n2 >= 2 && TString(token2[0]) == TString("Name")) {
	    Name = TString("");
            for(int i=1;i<n2;i++) {
	      Name += TString(token2[i]);
	      if(i < n2-1) Name += TString(" ");
	    }
	  }
          else if(n2 == 3 && TString(token2[0]) == TString("Thickness")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Thickness       = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("XOX0"))            XOX0            = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("Material"))        Material        = TString(token2[1]);
          else if(n2 == 3 && TString(token2[0]) == TString("widthU")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthU          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("widthV")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthV          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("Resolution")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Resolution      = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionU")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionU     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionV")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionV     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ROtime")) {
	    if(!global->IsTimeUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ROtime          = atof(token2[1])*global->GetTimeUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("Efficiency"))      DetEffic        = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token2[1]);
	  else if(n2 == 3 && TString(token2[0]) == TString("BkgRate")) {
	    if(!global->IsRateDensityUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    BkgRate         = atof(token2[1])*global->GetRateDensityUnit(TString(token2[2]));
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("SystemName")) {
	    SystemName = TString("");
            for(int i=1;i<n2;i++) {
	      SystemName += TString(token2[i]);
	      if(i < n2-1) SystemName += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("LayerName")) {
	    if(n2 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	    LayerName = TString(token2[1]);
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("ResolutionModel")) {
	    ResolutionModel = TString("");
            for(int i=1;i<n2;i++) {
	      ResolutionModel += TString(token2[i]);
	      if(i < n2-1) ResolutionModel += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("EfficiencyModel")) {
	    EfficiencyModel = TString("");
            for(int i=1;i<n2;i++) {
	      EfficiencyModel += TString(token2[i]);
	      if(i < n2-1) EfficiencyModel += TString(" ");
            }
	  }
          else if(n2 == 5 && TString(token2[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token2[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token2[4]);
	    Position = TVector3(atof(token2[1])*unit,atof(token2[2])*unit,atof(token2[3])*unit);
          }
          else if(n2 == 2 && TString(token2[0]) == TString("IsSensitive")) {
	    IsSensitive = global->SetBoolFromString(TString(token2[1])); 
	  }
	}
      
        if(!IsPlaneEnd) {
          cout << endl;
          cout << "Wrong specification of GeoPlane parameters of LadderPlane object inside the BeginPlane and EndPlane block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginPlane" << endl;
	  cout << "  Name         string (Mandatory)" << endl;
	  cout << "  Position        x       y       z    units (Mandatory)" << endl;
	  cout << "  Thickness       val units (Mandatory)" << endl;
	  cout << "  Material     string (Mandatory)" << endl;
	  cout << "  XOX0            val  (Optional)" << endl;
	  cout << "  widthU          val units (Mandatory)" << endl;
	  cout << "  widthV          val units (Mandatory)" << endl;
	  cout << "  IsSensitive     bool (Mandatory)" << endl;
	  cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	  cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	  cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	  cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	  cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	  cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	  cout << "  SystemName      string (Optional)" << endl;
	  cout << "  LayerName       string (Optional)" << endl;
	  cout << "  ResolutionModel string (Optional)" << endl;
	  cout << "  EfficiencyModel string (Optional)" << endl;
	  cout << "EndPlane" << endl;
          cout << "Check you inputs. Exiting now!!!" << endl;
          cout << endl;
	  
          assert(false);
        } // End While Plane
	
      } // end if BeginPlane
        
    } // End While LadderPlane
      
    if(!IsLadderEnd) {
      cout << endl;
      if(Ladder_counter > Ladder_Max) cout << "IMPORTANT: Reached Maximum number of lines (" << Ladder_Max << ") to read a LadderPlane!!!" << endl;
      cout << "Wrong specification of LadderPlane parameters inside the BeginLadderPlane and EndLadderPlane block. Parameters have to be specified as follows:" << endl;
      cout << "BeginLadderPlane" << endl;
      cout << "  LadderName        string (Mandatory)"           << endl;
      cout << "  LadderPosition      x       y       z    units (Mandatory)" << endl;
      cout << "  LadderRotAngles  alphaX  alphaY  alphaZ  units (Optional)" << endl;
      cout << "  BeginPlane" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Position        x       y       z    units (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    XOX0            val  (Optional)" << endl;
      cout << "    widthU          val units (Mandatory)" << endl;
      cout << "    widthV          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
      cout << "    InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndPlane" << endl;
      cout << "  ... " << endl;
      cout << "  BeginPlane" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Position        x       y       z    units (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    XOX0            val  (Optional)" << endl;
      cout << "    widthU          val units (Mandatory)" << endl;
      cout << "    widthV          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
      cout << "    InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndPlane" << endl;
      cout << "EndLadderPlane"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
    
  }
  else if(BlockName == TString("BeginLadderCylinder")) {
    bool IsLadderEnd    = false;
    int  Ladder_counter = 0;
        
    //Parameters
    int       NCylinders      = 0;
    TString   LadderName      = Dummy_name;
    TVector3  LadderPosition  = Dummy_vector;
    TVector3  LadderRotAngles = Dummy_vector;
    double    LadderRadius    = Dummy_value;

    while(!IsLadderEnd && Ladder_counter <= Ladder_Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      Ladder_counter++;
      
      if(n1 == 1 && TString(token1[0]) == TString("EndLadderCylinder")) {
	if(LadderRotAngles == Dummy_vector) {
	  cout << "WARNING in LadderCylinder Block: LadderRotAngles not specified. Setting all to zero." << endl;
	  LadderRotAngles = TVector3(0,0,0);
	}
	
	if(LadderName == Dummy_name) {
	  cout << "ERROR in LadderCylinder Block: Name not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderPosition == Dummy_vector) {
	  cout << "ERROR in LadderCylinder Block: LadderPosition parameter not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderRadius == Dummy_value) {
	  cout << "ERROR in LadderCylinder Block: LadderRadius parameter not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(NCylinders == 0) {
	  cout << "ERROR in LadderCylinder Block: NCylinders is zero!!!" << endl;
	  IsLadderEnd = false;
	}
	else IsLadderEnd = true;
      }
      else if(n1 >= 2 && TString(token1[0]) == TString("LadderName")) {
	LadderName = TString("");
        for(int i=1;i<n1;i++) {
	  LadderName += TString(token1[i]);
	  if(i < n1-1) LadderName += TString(" ");
        }
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderRadius")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderRadius = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 5 && TString(token1[0]) == TString("LadderPosition")) {
	if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit    = global->GetDistanceUnit(token1[4]);
	LadderPosition = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 5 && TString(token1[0]) == TString("LadderRotAngles")) {
	if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit     = global->GetAngleUnit(token1[4]);
	LadderRotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 1 && TString(token1[0]) == TString("BeginCylinder")) {
	bool IsCylinderEnd    = false;
        int  Cylinder_counter = 0;
        int  Cylinder_Max     = 24;
    
        TString   Name           = Dummy_name;
        double    Thickness      = Dummy_value;
        TString   Material       = Dummy_name;
	double    XOX0           = Dummy_value;
        double    Radius         = Dummy_value;
        double    Length         = Dummy_value;
        bool      IsSensitive    = false;
        double    Resolution     = Dummy_value;
        double    ResolutionU    = Dummy_value;
        double    ResolutionV    = Dummy_value;
        double    ROtime         = Dummy_value;
        double    DetEffic       = Dummy_value;
        double    InsensFracVneg = 0.0;
        double    InsensFracVpos = 0.0;
        double    BkgRate        = 0.0;
        TString   SystemName("");
	TString   LayerName("");
        TString   ResolutionModel("");
	TString   EfficiencyModel("");

        while(!IsCylinderEnd && Cylinder_counter <= Cylinder_Max && !fin->eof()) {
          char buf2[MAX_CHARS_PER_LINE];
          fin->getline(buf2, MAX_CHARS_PER_LINE);
	  line_number++;
          int n2 = 0;
          const char* token2[MAX_TOKENS_PER_LINE] = {};
      
          token2[0] = strtok(buf2, DELIMITER); // first token
          if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
    
          for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
            token2[n2] = strtok(0, DELIMITER);
	    if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          Ladder_counter++;
	  Cylinder_counter++;
          
          if(n2 == 1 && TString(token2[0]) == TString("EndCylinder")) {
	    bool GoodGeoObject = true;

	    if(LadderName == Dummy_name) {
	      cout << "ERROR in LadderCylinder Block: Name not specified!!!" << endl;
	      assert(false);
	    }
	    else if(LadderPosition == Dummy_vector) {
	      cout << "ERROR in LadderCylinder Block: LadderPosition parameter not specified!!!" << endl;
	      assert(false);
	    }
	    else if(LadderRadius == Dummy_value) {
	      cout << "ERROR in LadderCylinder Block: LadderRadius parameter not specified!!!" << endl;
	      assert(false);
	    }
	    
	    if(Name == Dummy_name) {
	      cout << "ERROR in Cylinder Block inside LadderCylinder: Name not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Thickness == Dummy_value) {
	      cout << "ERROR in Cylinder Block  inside LadderCylinder: Thickness not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    if(Material == Dummy_name) {
	      cout << "ERROR in Cylinder Block  inside LadderCylinder: Material not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Radius == Dummy_value) {
	      cout << "ERROR in Cylinder Block  inside LadderCylinder: Radius not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Length == Dummy_value) {
	      cout << "ERROR in Cylinder Block  inside LadderCylinder: Length not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    
	    if(LadderRotAngles == Dummy_vector) {
	      cout << "WARNING in LadderCylinder Block: LadderRotAngles not specified. Setting all to zero." << endl;
	      LadderRotAngles = TVector3(0,0,0);
	    }
	  
	    if(IsSensitive) {
	      if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	        if(Resolution == Dummy_value) {
		  cout << "ERROR in Cylinder Block  inside LadderCylinder: Resolution not specified!!!" << endl;
		  GoodGeoObject = false;
	        }
	        else {
		  ResolutionU = Resolution;
		  ResolutionV = Resolution;
	        }
	      }
	      else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	        cout << "ERROR in Cylinder Block  inside LadderCylinder: ResolutionU not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	        cout << "ERROR in Cylinder Block  inside LadderCylinder: ResolutionV not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ROtime == Dummy_value) {
	        cout << "ERROR in Cylinder Block  inside LadderCylinder: ROtime not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(DetEffic == Dummy_value) {
	        cout << "ERROR in Cylinder Block  inside LadderCylinder: DetEffic not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	    }
	  
	    if(GoodGeoObject) {
	      NCylinders++;
	      
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(LadderRotAngles,Rot);
	      TString FullName = LadderName + TString(" ") + Name;
	      
	      TVector3  GlobalPosition = LadderPosition;
	      double    GlobalRadius   = LadderRadius + Radius;
	      
	      GGeoObject* aGeoObject = new GGeoCylinder(FullName,aGeometry->GetGeometryID(),0,
						        GlobalPosition,Rot,Thickness,Material,
					                IsSensitive,
					                Length,GlobalRadius,
						        TVector2(InsensFracVneg,InsensFracVpos),
						        global,
					                ResolutionU,ResolutionV,DetEffic,ROtime,
					                BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      aGeoObject->SetLadderType(TString("LadderCylinder"));
	      if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	      aGeometry->PushGeoElement(aGeoObject);
	      
	      IsCylinderEnd = true;
	    }

          }
          else if(n2 >= 2 && TString(token2[0]) == TString("Name")) {
	    Name = TString("");
            for(int i=1;i<n2;i++) {
	      Name += TString(token2[i]);
	      if(i < n2-1) Name += TString(" ");
	    }
	  }
          else if(n2 == 3 && TString(token2[0]) == TString("Thickness")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Thickness       = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("XOX0"))            XOX0            = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("Material"))        Material        = TString(token2[1]);
          else if(n2 == 3 && TString(token2[0]) == TString("Radius")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Radius          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("Length")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Length          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("Resolution")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Resolution      = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionU")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionU     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionV")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionV     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ROtime")) {
	    if(!global->IsTimeUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ROtime          = atof(token2[1])*global->GetTimeUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("Efficiency"))      DetEffic        = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token2[1]);
	  else if(n2 == 3 && TString(token2[0]) == TString("BkgRate")) {
	    if(!global->IsRateDensityUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    BkgRate         = atof(token2[1])*global->GetRateDensityUnit(TString(token2[2]));
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("SystemName")) {
	    SystemName = TString("");
            for(int i=1;i<n2;i++) {
	      SystemName += TString(token2[i]);
	      if(i < n2-1) SystemName += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("LayerName")) {
	    if(n2 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	    LayerName = TString(token2[1]);
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("ResolutionModel")) {
	    ResolutionModel = TString("");
            for(int i=1;i<n2;i++) {
	      ResolutionModel += TString(token2[i]);
	      if(i < n2-1) ResolutionModel += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("EfficiencyModel")) {
	    EfficiencyModel = TString("");
            for(int i=1;i<n2;i++) {
	      EfficiencyModel += TString(token2[i]);
	      if(i < n2-1) EfficiencyModel += TString(" ");
            }
	  }
          else if(n2 == 2 && TString(token2[0]) == TString("IsSensitive")) {
	    IsSensitive = global->SetBoolFromString(TString(token2[1])); 
	  }
	}
      
        if(!IsCylinderEnd) {
          cout << endl;
          cout << "Wrong specification of GeoCylinder parameters of LadderCylinder object inside the BeginGeoCylinder and EndGeoCylinder block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginCylinder" << endl;
	  cout << "  Name         string (Mandatory)" << endl;
	  cout << "  Thickness       val units (Mandatory)" << endl;
	  cout << "  Material     string (Mandatory)" << endl;
	  cout << "  XOX0            val  (Optional)" << endl;
	  cout << "  Radius          val units (Mandatory)" << endl;
	  cout << "  Length          val units (Mandatory)" << endl;
	  cout << "  IsSensitive     bool (Mandatory)" << endl;
	  cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	  cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	  cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	  cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	  cout << "  SystemName      string (Optional)" << endl;
	  cout << "  LayerName       string (Optional)" << endl;
	  cout << "  ResolutionModel string (Optional)" << endl;
	  cout << "  EfficiencyModel string (Optional)" << endl;
	  cout << "EndCylinder" << endl;
          cout << "Check you inputs. Exiting now!!!" << endl;
          cout << endl;
	  
          assert(false);
        } // End While Cylinder
	
      } // end if BeginCylinder
        
    } // End While LadderCylinder
      
    if(!IsLadderEnd) {
      cout << endl;
      if(Ladder_counter > Ladder_Max) cout << "IMPORTANT: Reached Maximum number of lines (" << Ladder_Max << ") to read a LadderPlane!!!" << endl;
      cout << "Wrong specification of LadderCylinder parameters inside the BeginLadderCylinder and EndLadderCylinder block. Parameters have to be specified as follows:" << endl;
      cout << "BeginLadderCylinder" << endl;
      cout << "  LadderName        string (Mandatory)"           << endl;
      cout << "  LadderPosition      x       y       z    units (Mandatory)" << endl;
      cout << "  LadderRotAngles  alphaX  alphaY  alphaZ  units (Optional)" << endl;
      cout << "  LadderRadius        val   units  (Mandatory)"  << endl;
      cout << "  BeginCylinder" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    XOX0            val  (Optional)" << endl;
      cout << "    Radius          val units (Mandatory)" << endl;
      cout << "    Length          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndCylinder" << endl;
      cout << "  ...  " << endl;
      cout << "  BeginCylinder" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    XOX0            val  (Optional)" << endl;
      cout << "    Radius          val units (Mandatory)" << endl;
      cout << "    Length          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndCylinder" << endl;
      cout << "EndLadderCylinder"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }

  }
  else if(BlockName == TString("BeginMosaicLadderPlane")) {
    bool IsLadderEnd    = false;
    int  Ladder_counter = 0;
    
    //Parameters
    TString   LadderName      = Dummy_name;
    TString   LadderMosaic    = Dummy_name;
    TVector3  LadderPosition  = Dummy_vector;
    TVector3  LadderRotAngles = Dummy_vector;
    double    LadderRadius    = Dummy_value;
    //Variables only used in case of Spiral or Alternating geometries
    double  LadderClearanceH = 1.0*global->GetDistanceUnit("um");
    double  LadderClearanceV = 1.0*global->GetDistanceUnit("um");
    double  LadderRmin       = Dummy_value;
    double  LadderRmax       = Dummy_value;
    double  LadderThickness  = Dummy_value;
    double  LadderOverlap    = Dummy_value;
    double  LadderWidth      = Dummy_value;
    double  LadderAlpha0     = 0.0;
    double  LadderAlpha      = 0.0;
    double  LadderShift      = Dummy_value;
    int     LadderN          = 0;
    bool    LadderShiftFix   = false;
    TString LadderVarPar     = TString("Overlap");
    
    //List with template GGeoPlane Objects used as basis to built mosaic
    std::vector<GGeoObject*>  LadderPlanesList;
    LadderPlanesList.clear();
    
    while(!IsLadderEnd && Ladder_counter <= Ladder_Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      Ladder_counter++;

      if(n1 == 1 && TString(token1[0]) == TString("EndMosaicLadderPlane")) {
	if(LadderRotAngles == Dummy_vector) {
	  cout << "WARNING in MosaicLadderPlane Block: LadderRotAngles not specified. Setting all to zero." << endl;
	  LadderRotAngles = TVector3(0,0,0);
	}
	
	if(LadderName == Dummy_name) {
	  cout << "ERROR in MosaicLadderPlane Block: Name not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderMosaic == Dummy_name) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderMosaic type not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderMosaic != TString("Spiral") && LadderMosaic != TString("Alternating")) {
	  cout << "ERROR in MosaicLadderPlane Block: unknown LadderMosaic type " << LadderMosaic.Data() << ". Options are Spiral and Alternating!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderPosition == Dummy_vector) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderPosition parameter not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderClearanceH < 0.0) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderClearanceH is negative!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderClearanceV < 0.0) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderClearanceV is negative!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderAlpha0 < 0.0 || LadderAlpha0 > 360.0*global->GetAngleUnit("deg")) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderAlpha0 is outside (0,360) deg range!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderVarPar != TString("Overlap") && LadderVarPar != TString("Radius")) {
	  cout << endl;
          cout << "WARNING: LadderVarPar can only have values: Overlap, Width or Radius. Setting it to default value = Overlap !!!" << endl;
          cout << endl;
	  LadderVarPar = TString("Overlap");
        }
        else if(LadderVarPar != TString("Overlap") && LadderOverlap == Dummy_value) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderOverlap parameter not being varied and not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderVarPar != TString("Radius") && LadderRadius == Dummy_value) {
	  cout << "ERROR in MosaicLadderPlane Block: LadderRadius parameter not being varied and not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderN == 0) {
	  cout << "ERROR in MosaicLadderPlane Block: zero LadderN, which has to be > 0!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderMosaic == TString("Alternating") && LadderN%2) {
	  cout << "WARNLadderShiftING in MosaicLadderPlane Block: for LadderMosaic = Alternating LadderN has to be an even number. Setting it to from " << LadderN << " to " << LadderN+1<< "." << endl;
	  LadderN++;
	}
	else IsLadderEnd = true;
	
	if(LadderShiftFix == true && LadderShift == Dummy_value) {
	  cout << "ERROR in MosaicLadderPlane Block: When LadderShiftFix is true LadderShift has to be specified!!!" << endl;
	  IsLadderEnd = false;
	}
	
	if(IsLadderEnd == true) {
	  
	  TMatrixD LadderRot;
          global->GetGlobalRotationMatrix(LadderRotAngles,LadderRot);
          TMatrixD LadderInvRot;
          LadderInvRot.ResizeTo(3,3);
          LadderInvRot = LadderRot;
          LadderInvRot.Invert();
          TVector3 UVector(1.0,0.0,0.0);
          TVector3 VVector(0.0,1.0,0.0);
          TVector3 WVector(0.0,0.0,1.0);
          global->RotateVector(LadderInvRot,UVector);
          global->RotateVector(LadderInvRot,VVector);
          global->RotateVector(LadderInvRot,WVector);
	  
	  //Now calculate the ladder minimum radius, thickness and width
	  LadderRmin       = +1.0e+20;
	  LadderRmax       = -1.0e+20;
	  LadderThickness  = 0.0;
	  LadderWidth      = -1.0e+20;
	  for(int iplane=0;iplane<int(LadderPlanesList.size());iplane++) {
	    GGeoPlane* aPlane = dynamic_cast<GGeoPlane*>(LadderPlanesList[iplane]);
	    TVector3   PlanePosition = aPlane->GetMainSurface()->GetPosition();

	    double r  = LadderRadius;
	    r        += PlanePosition.Z();
	    double t  = 0.5*aPlane->GetThickness();
	    double wU = 0.5*aPlane->GetUVWidth().X();
	      
	    if(wU <= 0) {
	      cout << endl;
	      cout << "For the Spiral and Alternating LadderType need to specify a Uwidth for each plane of the ladder. For ladder with name " << LadderName.Data()
		   << " its plane number " << iplane+1 << " with name " << aPlane->GetName().Data() << " has Uwidth = " << wU/global->GetDistanceUnit("cm") << " cm. "
		   << "Check your inputs. Exiting now!!!" << endl;
	      cout << endl;
	      assert(false);
	    }
	    
	    double wU1 = TMath::Max(TMath::Abs(PlanePosition.X() - wU),TMath::Abs(PlanePosition.X() + wU));
	    
	    if(LadderRmin  > r - t) LadderRmin  = r - t;
	    if(LadderRmax  < r + t) LadderRmax  = r + t;
	    if(LadderWidth < wU1)   LadderWidth = wU1;
	  }
	  LadderThickness  = LadderRmax - LadderRmin;
	  LadderThickness += LadderClearanceV;
	  LadderWidth     *= 2;
	  LadderWidth     += LadderClearanceH;
	  double delta_R   = LadderRadius - LadderRmin;

	  //Now calculate the values of the parameters to dispose ladders in position
	  global->GetMosaicGeoParams(LadderMosaic,LadderRmin,LadderThickness,LadderOverlap,LadderN,LadderWidth,LadderShift,LadderAlpha,LadderVarPar,LadderShiftFix);
	  
	  if(verbose || true) {
	    cout << endl;
	    cout << "Layer is " << LadderMosaic.Data() << " with name " << LadderName.Data()     << " : " << endl;
	    cout << " - Ladder radius      = " << LadderRadius/global->GetDistanceUnit("cm")     << " cm"  << endl;
	    cout << " - Ladder clearance H = " << LadderClearanceH/global->GetDistanceUnit("um") << " um"  << endl;
	    cout << " - Ladder clearance V = " << LadderClearanceV/global->GetDistanceUnit("um") << " um"  << endl;
	    cout << " - Ladder alpha0      = " << LadderAlpha0/global->GetAngleUnit("deg")       << " deg" << endl;
	    cout << " - Ladder alpha       = " << LadderAlpha/global->GetAngleUnit("deg")        << " deg" << endl;
	    cout << " - Ladder overlap     = " << LadderOverlap*100                              << " %"   << endl;
	    cout << " - Ladder Var Par     = " << LadderVarPar.Data()                            << endl;
	    cout << " - Ladder N           = " << LadderN                                        << endl;
	    cout << " - Ladder Thickness   = " << LadderThickness/global->GetDistanceUnit("um")  << " um"  << endl;
	    cout << " - Ladder Width       = " << LadderWidth/global->GetDistanceUnit("mm")      << " mm"  << endl;
	    cout << " - Ladder Rmin        = " << LadderRmin/global->GetDistanceUnit("cm")       << " cm"  << endl;
	    cout << " - Ladder Rmax        = " << LadderRmax/global->GetDistanceUnit("cm")       << " cm"  << endl;
	    cout << " - Ladder deltaR      = " << delta_R/global->GetDistanceUnit("um")          << " um"  << endl;
	    cout << " - Ladder shift       = " << LadderShift                                    << endl;
	    cout << endl;
	  }
	  
	  std::vector<double>  angle_rot;
	  std::vector<TString> axis_rot;
	  angle_rot.clear();
	  axis_rot.clear();
	  axis_rot.push_back(TString("X")); angle_rot.push_back(-90.0*global->GetAngleUnit("deg"));
	  axis_rot.push_back(TString("Y")); angle_rot.push_back(  0.0*global->GetAngleUnit("deg"));
	  
	  for(int ialpha=0;ialpha<LadderN;ialpha++) { //Begin of loop over ladders for this radius
	    double alpha_tmp = LadderAlpha0 + LadderAlpha*ialpha;
	    angle_rot[1] = +alpha_tmp;
	    double cos_alpha_tmp = TMath::Cos(alpha_tmp);
            double sin_alpha_tmp = TMath::Sin(alpha_tmp);

	    TMatrixD Rot;
	    Rot.ResizeTo(3,3);
	    global->GetGlobalRotationMatrix_FromList(axis_rot,angle_rot,Rot);
	    TVector3 RotAngles = global->GetRotationAngles(Rot*LadderRot);

	    //double Ladder_R_tmp;
	    TVector3 MyLadderPosition;
	    if(LadderMosaic == TString("Spiral")) {
	      double r,w,fx,delta;
	      r     = LadderRmin;
	      w     = LadderWidth;
	      fx    = LadderShift;
	      delta = w*fx;
		
	      double x0 = 0.5*w - delta;
	      double y0 = r;
		
	      double x =  x0*cos_alpha_tmp + y0*sin_alpha_tmp;
	      double y = -x0*sin_alpha_tmp + y0*cos_alpha_tmp;
	      double z =  0.0;
	      
	      //Ladder_R_tmp = r;
	      MyLadderPosition = x*UVector + y*VVector + z*WVector + LadderPosition;
	    }
	    else if(LadderMosaic == TString("Alternating")) {
	      double r,fr;
	      r     = LadderRmin;
	      fr    = LadderShift;
	      if(ialpha%2) r *= fr;
		
	      double x0 = 0.0;
	      double y0 = r;
		
	      double x =  x0*cos_alpha_tmp + y0*sin_alpha_tmp;
	      double y = -x0*sin_alpha_tmp + y0*cos_alpha_tmp;
	      double z =  0.0;
		
	      //Ladder_R_tmp = r;
	      MyLadderPosition = x*UVector + y*VVector + z*WVector + LadderPosition;
	    }
	    
	    char name[300];
	    sprintf(name,"%.1f",alpha_tmp/global->GetAngleUnit("deg"));
	    
	    for(int iplane=0;iplane<int(LadderPlanesList.size());iplane++) {
	      GGeoPlane* aPlane = dynamic_cast<GGeoPlane*>(LadderPlanesList[iplane]);
	      
	      TString   Name           = LadderName + TString(", Plane ") + aPlane->GetName() + TString(" at alpha = ") + TString(name) + TString(" deg");
	      TMatrixD aRot;
	      global->GetGlobalRotationMatrix(RotAngles,aRot);
	      TMatrixD aInvRot;
	      aInvRot.ResizeTo(3,3);
	      aInvRot = aRot;
	      aInvRot.Invert();
	      TVector3  Position       = aPlane->GetMainSurface()->GetPosition();
	      global->RotateVector(aInvRot,Position);
	      Position += MyLadderPosition;
	      
	      double    Thickness       = aPlane->GetThickness();
	      TString   Material        = aPlane->GetMaterial();
	      double    XOX0            = aPlane->GetXOX0();
	      double    widthU          = aPlane->GetUVWidth().X();
	      double    widthV          = aPlane->GetUVWidth().Y();
	      bool      IsSensitive     = aPlane->GetIsSensitive();
	      double    ResolutionU     = aPlane->GetResolutionU();
	      double    ResolutionV     = aPlane->GetResolutionV();
	      double    ROtime          = aPlane->GetROtime();
	      double    DetEffic        = aPlane->GetDetEfficiency();
	      double    InsensFracUneg  = aPlane->GetUInsensitive().X();
	      double    InsensFracUpos  = aPlane->GetUInsensitive().Y();
	      double    InsensFracVneg  = aPlane->GetVInsensitive().X();
	      double    InsensFracVpos  = aPlane->GetVInsensitive().Y();
	      double    BkgRate         = aPlane->GetBkgRate();
	      TString   SystemName      = aPlane->GetSystemName();
	      TString   LayerName       = aPlane->GetLayerName();
	      TString   ResolutionModel = aPlane->GetResolutionModel();
	      TString   EfficiencyModel = aPlane->GetEfficiencyModel();
	      
	      GGeoObject* aGeoObject = new GGeoPlane(Name,aGeometry->GetGeometryID(),0,
						     Position,aRot,Thickness,Material,
					             IsSensitive,
					             TVector2(widthU,widthV),
						     TVector2(InsensFracUneg,InsensFracUpos),
						     TVector2(InsensFracVneg,InsensFracVpos),
						     global,
					             ResolutionU,ResolutionV,DetEffic,ROtime,
					             BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      aGeoObject->SetLadderType(aPlane->GetLadderType());
	      aGeoObject->SetXOX0(XOX0);
	      aGeometry->PushGeoElement(aGeoObject);
	    }
	    
	  } //End of loop over ladders for this radius
	  
	  //Freeing memory of the template GGeoPlane Objects
	  for(int kkk=0;kkk<int(LadderPlanesList.size());kkk++) delete LadderPlanesList[kkk];
	  LadderPlanesList.clear();
	}
	
      }
      else if(n1 >= 2 && TString(token1[0]) == TString("LadderName")) {
	LadderName = TString("");
        for(int i=1;i<n1;i++) {
	  LadderName += TString(token1[i]);
	  if(i < n1-1) LadderName += TString(" ");
        }
      }
      else if(n1 >= 2 && TString(token1[0]) == TString("LadderMosaic")) {
	LadderMosaic = TString("");
        for(int i=1;i<n1;i++) {
	  LadderMosaic += TString(token1[i]);
	  if(i < n1-1) LadderMosaic += TString(" ");
        }
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderRadius")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderRadius     = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderClearanceH")){
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderClearanceH = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderClearanceV")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderClearanceV = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 2 && TString(token1[0]) == TString("LadderOverlap"))     LadderOverlap    = atof(token1[1]);
      else if(n1 == 3 && TString(token1[0]) == TString("LadderRadius")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderRadius     = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderAlpha0")) {
	if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderAlpha0     = atof(token1[1])*global->GetAngleUnit(token1[2]);
      }
      else if(n1 == 2 && TString(token1[0]) == TString("LadderN"))           LadderN          = atoi(token1[1]);
      else if(n1 == 2 && TString(token1[0]) == TString("LadderShift"))       LadderShift      = atof(token1[1]);
      else if(n1 == 2 && TString(token1[0]) == TString("LadderVarPar"))      LadderVarPar     = TString(token1[1]);
      else if(n1 == 5 && TString(token1[0]) == TString("LadderPosition")) {
	if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit    = global->GetDistanceUnit(token1[4]);
	LadderPosition = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 5 && TString(token1[0]) == TString("LadderRotAngles")) {
	if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit     = global->GetAngleUnit(token1[4]);
	LadderRotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 2 && TString(token1[0]) == TString("LadderShiftFix")) {
	LadderShiftFix = global->SetBoolFromString(TString(token1[1])); 
      }
      else if(n1 == 1 && TString(token1[0]) == TString("BeginPlane")) {
	bool IsPlaneEnd    = false;
        int  Plane_counter = 0;
        int  Plane_Max     = 24;
    
	TString   Name           = Dummy_name;
        TVector3  Position       = Dummy_vector;
        double    Thickness      = Dummy_value;
	double    XOX0           = Dummy_value;
        TString   Material       = Dummy_name;
        double    widthU         = Dummy_value;
        double    widthV         = Dummy_value;
        bool      IsSensitive    = false;
        double    Resolution     = Dummy_value;
        double    ResolutionU    = Dummy_value;
        double    ResolutionV    = Dummy_value;
        double    ROtime         = Dummy_value;
        double    DetEffic       = Dummy_value;
        double    InsensFracUneg = 0.0;
        double    InsensFracUpos = 0.0;
        double    InsensFracVneg = 0.0;
        double    InsensFracVpos = 0.0;
        double    BkgRate        = 0.0;
        TString   SystemName("");
	TString   LayerName("");
        TString   ResolutionModel("");
	TString   EfficiencyModel("");

        while(!IsPlaneEnd && Plane_counter <= Plane_Max && !fin->eof()) {
          char buf2[MAX_CHARS_PER_LINE];
          fin->getline(buf2, MAX_CHARS_PER_LINE);
	  line_number++;
          int n2 = 0;
          const char* token2[MAX_TOKENS_PER_LINE] = {};
      
          token2[0] = strtok(buf2, DELIMITER); // first token
          if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
    
          for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
            token2[n2] = strtok(0, DELIMITER);
	    if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          Ladder_counter++;
	  Plane_counter++;
          
          if(n2 == 1 && TString(token2[0]) == TString("EndPlane")) {
	    bool GoodGeoObject = true;

	    if(Name == Dummy_name) {
	      cout << "ERROR in Plane Block inside MosaicLadderPlane: Name not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Position == Dummy_vector) {
	      cout << "ERROR in Plane Block inside MosaicLadderPlane: Position not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Thickness == Dummy_value) {
	      cout << "ERROR in Plane Block  inside MosaicLadderPlane: Thickness not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    if(Material == Dummy_name) {
	      cout << "ERROR in Plane Block  inside MosaicLadderPlane: Material not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(widthU == Dummy_value) {
	      cout << "ERROR in Plane Block  inside MosaicLadderPlane: widthU not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(widthV == Dummy_value) {
	      cout << "ERROR in Plane Block  inside MosaicLadderPlane: widthV not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    
	    if(LadderRotAngles == Dummy_vector) {
	      cout << "WARNING in Plane Block  inside MosaicLadderPlane: RotAngles not specified. Setting all to zero." << endl;
	      LadderRotAngles = TVector3(0,0,0);
	    }
	  
	    if(IsSensitive) {
	      if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	        if(Resolution == Dummy_value) {
		  cout << "ERROR in Plane Block  inside MosaicLadderPlane: Resolution not specified!!!" << endl;
		  GoodGeoObject = false;
	        }
	        else {
		  ResolutionU = Resolution;
		  ResolutionV = Resolution;
	        }
	      }
	      else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	        cout << "ERROR in Plane Block  inside MosaicLadderPlane: ResolutionU not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	        cout << "ERROR in Plane Block  inside MosaicLadderPlane: ResolutionV not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ROtime == Dummy_value) {
	        cout << "ERROR in Plane Block  inside MosaicLadderPlane: ROtime not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(DetEffic == Dummy_value) {
	        cout << "ERROR in Plane Block  inside MosaicLadderPlane: DetEffic not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	    }
	  
	    if(GoodGeoObject) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(LadderRotAngles,Rot);
	      
	      TVector3  GlobalPosition = Position;
	      
	      GGeoObject* aGeoObject = new GGeoPlane(Name,aGeometry->GetGeometryID(),0,
						     GlobalPosition,Rot,Thickness,Material,
					             IsSensitive,
					             TVector2(widthU,widthV),
						     TVector2(InsensFracUneg,InsensFracUpos),
						     TVector2(InsensFracVneg,InsensFracVpos),
						     global,
					             ResolutionU,ResolutionV,DetEffic,ROtime,
					             BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      aGeoObject->SetLadderType(LadderMosaic);
	      if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	      LadderPlanesList.push_back(aGeoObject);
	      
	      IsPlaneEnd = true;
	    }

          }
          else if(n2 >= 2 && TString(token2[0]) == TString("Name")) {
	    Name = TString("");
            for(int i=1;i<n2;i++) {
	      Name += TString(token2[i]);
	      if(i < n2-1) Name += TString(" ");
	    }
	  }
          else if(n2 == 3 && TString(token2[0]) == TString("Thickness")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Thickness       = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("XOX0"))            XOX0            = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("Material"))        Material        = TString(token2[1]);
          else if(n2 == 3 && TString(token2[0]) == TString("widthU")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthU          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("widthV")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    widthV          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("Resolution")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Resolution      = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionU")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionU     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionV")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionV     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ROtime")) {
	    if(!global->IsTimeUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ROtime          = atof(token2[1])*global->GetTimeUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("Efficiency"))      DetEffic        = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token2[1]);
	  else if(n2 == 3 && TString(token2[0]) == TString("BkgRate")) {
	    if(!global->IsRateDensityUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    BkgRate         = atof(token2[1])*global->GetRateDensityUnit(TString(token2[2]));
	  }
	  else if(n2 == 5 && TString(token2[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token2[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token2[4]);
	    Position = TVector3(atof(token2[1])*unit,atof(token2[2])*unit,atof(token2[3])*unit);
          }
	  else if(n2 >= 2 && TString(token2[0]) == TString("SystemName")) {
	    SystemName = TString("");
            for(int i=1;i<n2;i++) {
	      SystemName += TString(token2[i]);
	      if(i < n2-1) SystemName += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("LayerName")) {
	    if(n2 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	    LayerName = TString(token2[1]);
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("ResolutionModel")) {
	    ResolutionModel = TString("");
            for(int i=1;i<n2;i++) {
	      ResolutionModel += TString(token2[i]);
	      if(i < n2-1) ResolutionModel += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("EfficiencyModel")) {
	    EfficiencyModel = TString("");
            for(int i=1;i<n2;i++) {
	      EfficiencyModel += TString(token2[i]);
	      if(i < n2-1) EfficiencyModel += TString(" ");
            }
	  }
          else if(n2 == 2 && TString(token2[0]) == TString("IsSensitive")) {
	    IsSensitive = global->SetBoolFromString(TString(token2[1])); 
	  }
	}
      
        if(!IsPlaneEnd) {
          cout << endl;
          cout << "Wrong specification of GeoPlane parameters of MosaicLadderPlane object inside the BeginPlane and EndPlane block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginPlane" << endl;
	  cout << "  Name         string (Mandatory)" << endl;
	  cout << "  Position        x       y       z    units (Mandatory)" << endl;
	  cout << "  Thickness       val units (Mandatory)" << endl;
	  cout << "  Material     string (Mandatory)" << endl;
	  cout << "  XOX0            val   (Optional)" << endl;
	  cout << "  widthU          val units (Mandatory)" << endl;
	  cout << "  widthV          val units (Mandatory)" << endl;
	  cout << "  IsSensitive     bool (Mandatory)" << endl;
	  cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	  cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	  cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	  cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	  cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	  cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	  cout << "  SystemName      string (Optional)" << endl;
	  cout << "  LayerName       string (Optional)" << endl;
	  cout << "  ResolutionModel string (Optional)" << endl;
	  cout << "  EfficiencyModel string (Optional)" << endl;
	  cout << "EndPlane" << endl;
          cout << "Check you inputs. Exiting now!!!" << endl;
          cout << endl;
	  
          assert(false);
	  
        } // End While Plane

      } // end if BeginPlane

    } // End While MosaicLadderPlane
    
    if(!IsLadderEnd) {
      cout << endl;
      if(Ladder_counter > Ladder_Max) cout << "IMPORTANT: Reached Maximum number of lines (" << Ladder_Max << ") to read a LadderPlane!!!" << endl;
      cout << "Wrong specification of MosaicLadderPlane parameters inside the BeginMosaicLadderPlane and EndMosaicLadderPlane block. Parameters have to be specified as follows:" << endl;
      cout << "BeginMosaicLadderPlane" << endl;
      cout << "  LadderName        string (Mandatory)"           << endl;
      cout << "  LadderMosaic      string (Mandatory, either Spiral or Alternating)" << endl;
      cout << "  LadderPosition      x       y       z    units (Mandatory)" << endl;
      cout << "  LadderRotAngles  alphaX  alphaY  alphaZ  units (Optional)" << endl;
      cout << "  LadderVarPar     string  (Optional, can take values Overlap and Radius, default is Overlap)" << endl;
      cout << "  LadderRadius       val   units  (Mandatory if LadderVarPar = Overlap)" << endl;
      cout << "  LadderOverlap      val   units  (Mandatory if LadderVarPar = Radius)" << endl;
      cout << "  LadderClearanceH   val   units  (Optional)" << endl;
      cout << "  LadderClearanceV   val   units  (Optional)" << endl;
      cout << "  LadderAlpha0       val   units  (Optional, default value is zero)" << endl;
      cout << "  LadderN            N            (Mandatory and > 0)" << endl;
      cout << "  LadderShift        value        (Optional)" << endl;
      cout << "  LadderShiftFix     bool         (Optional)" << endl;
      cout << "  BeginPlane" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Position        x       y       z    units (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    XOX0            val   (Optional)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    widthU          val units (Mandatory)" << endl;
      cout << "    widthV          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
      cout << "    InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndPlane" << endl;
      cout << "  ...  " << endl;
      cout << "  BeginPlane" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Position        x       y       z    units (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    XOX0            val   (Optional)" << endl;
      cout << "    widthU          val units (Mandatory)" << endl;
      cout << "    widthV          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
      cout << "    InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndPlane" << endl;
      cout << "EndMosaicLadderPlane"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }    
  }
  else if(BlockName == TString("BeginMosaicLadderPetal")) {
    bool IsLadderEnd    = false;
    int  Ladder_counter = 0;
    
    //Parameters
    TString   LadderName       = Dummy_name;
    TVector3  LadderPosition   = Dummy_vector;
    TVector3  LadderRotAngles  = Dummy_vector;
    double    LadderRadius     = Dummy_value;
    double    LadderAlpha0     = 0.0;
    double    LadderAlpha      = 0.0;
    int       LadderN          = 0;
    double    LadderShift      = Dummy_value;
    
    //List with template GGeoPlane Objects used as basis to built mosaic
    std::vector<GGeoObject*>  LadderPetalsList;
    LadderPetalsList.clear();
    
    while(!IsLadderEnd && Ladder_counter <= Ladder_Max && !fin->eof()) {
      char buf1[MAX_CHARS_PER_LINE];
      fin->getline(buf1, MAX_CHARS_PER_LINE);
      line_number++;
      int n1 = 0;
      const char* token1[MAX_TOKENS_PER_LINE] = {};
      
      token1[0] = strtok(buf1, DELIMITER); // first token
      if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
      for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
        token1[n1] = strtok(0, DELIMITER);
	if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      Ladder_counter++;

      if(n1 == 1 && TString(token1[0]) == TString("EndMosaicLadderPetal")) {
	if(LadderRotAngles == Dummy_vector) {
	  cout << "WARNING in MosaicLadderPetal Block: LadderRotAngles not specified. Setting all to zero." << endl;
	  LadderRotAngles = TVector3(0,0,0);
	}
	
	if(LadderName == Dummy_name) {
	  cout << "ERROR in MosaicLadderPetal Block: Name not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderPosition == Dummy_vector) {
	  cout << "ERROR in MosaicLadderPetal Block: LadderPosition parameter not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderAlpha0 < 0.0 || LadderAlpha0 > 360.0*global->GetAngleUnit("deg")) {
	  cout << "ERROR in MosaicLadderPetal Block: LadderAlpha0 is outside (0,360) deg range!!!" << endl;
	  IsLadderEnd = false;
	}
	else if(LadderN == 0) {
	  cout << "ERROR in MosaicLadderPetal Block: zero LadderN, which has to be > 0!!!" << endl;
	  IsLadderEnd = false;
	}
	if(LadderShift == Dummy_value) {
	  cout << "ERROR in MosaicLadderPetal Block: LadderShift not specified!!!" << endl;
	  IsLadderEnd = false;
	}
	else IsLadderEnd = true;
	
        if(IsLadderEnd == true) {
	  LadderAlpha          = 360.0*global->GetAngleUnit("deg")/LadderN;
	  double LadderAlpha02 = LadderAlpha0 + 0.5*LadderAlpha;
	  
	  TMatrixD LadderRot;
          global->GetGlobalRotationMatrix(LadderRotAngles,LadderRot);
          TMatrixD LadderInvRot;
          LadderInvRot.ResizeTo(3,3);
          LadderInvRot = LadderRot;
          LadderInvRot.Invert();
          TVector3 UVector(1.0,0.0,0.0);
          TVector3 VVector(0.0,1.0,0.0);
          TVector3 WVector(0.0,0.0,1.0);
          global->RotateVector(LadderInvRot,UVector);
          global->RotateVector(LadderInvRot,VVector);
          global->RotateVector(LadderInvRot,WVector);
	  
	  for(int ialpha=0;ialpha<LadderN;ialpha++) { //Begin of loop over ladders for this radius
	    char name[300];
	    double alpha_tmp;
	    double cos_alpha_tmp;
            double sin_alpha_tmp;
	    TMatrixD Rot;
	    Rot.ResizeTo(3,3);
	    TMatrixD Rot0;
	    global->GetRotationMatrix(-90*global->GetAngleUnit("deg"),TString("Z"),Rot0);
	    Rot0.ResizeTo(3,3);
	    
	    TVector3 MyPetalPosition;
	    TVector3 RotAngles;
	    
	    //Positive shift Petals
	    alpha_tmp = LadderAlpha0 + LadderAlpha*ialpha;
	    cos_alpha_tmp = TMath::Cos(alpha_tmp);
            sin_alpha_tmp = TMath::Sin(alpha_tmp);
	    global->GetRotationMatrix(alpha_tmp,TString("Z"),Rot);
	    Rot = Rot0*Rot;
	    RotAngles       = global->GetRotationAngles(Rot*LadderRot);
	    MyPetalPosition = LadderPosition + LadderRadius*cos_alpha_tmp*UVector + LadderRadius*sin_alpha_tmp*VVector + (+0.5*LadderShift)*WVector;
	    
	    sprintf(name,"%.1f",alpha_tmp/global->GetAngleUnit("deg"));
	    for(int iplane=0;iplane<int(LadderPetalsList.size());iplane++) {
	      GGeoPetal* aPetal = dynamic_cast<GGeoPetal*>(LadderPetalsList[iplane]);
	      
	      TString   Name           = LadderName + TString(", Disk ") + aPetal->GetName() + TString(" positive shift at alpha = ") + TString(name) + TString(" deg");
	      TMatrixD aRot;
	      global->GetGlobalRotationMatrix(RotAngles,aRot);
	      TMatrixD aInvRot;
	      aInvRot.ResizeTo(3,3);
	      aInvRot = aRot;
	      aInvRot.Invert();
	      TVector3  Position = aPetal->GetMainSurface()->GetPosition();
	      global->RotateVector(aInvRot,Position);
	      Position += MyPetalPosition;
	      
	      double    Thickness       = aPetal->GetThickness();
	      TString   Material        = aPetal->GetMaterial();
	      double    XOX0            = aPetal->GetXOX0();
	      double    bottomWidth     = aPetal->GetPetalBaseWidth();
	      double    topWidth        = aPetal->GetPetalTopWidth();
	      double    Height          = aPetal->GetPetalHeight();
	      bool      IsSensitive     = aPetal->GetIsSensitive();
	      double    ResolutionU     = aPetal->GetResolutionU();
	      double    ResolutionV     = aPetal->GetResolutionV();
	      double    ROtime          = aPetal->GetROtime();
	      double    DetEffic        = aPetal->GetDetEfficiency();
	      double    InsensFracUneg  = aPetal->GetUInsensitive().X();
	      double    InsensFracUpos  = aPetal->GetUInsensitive().Y();
	      double    InsensFracVneg  = aPetal->GetVInsensitive().X();
	      double    InsensFracVpos  = aPetal->GetVInsensitive().Y();
	      double    BkgRate         = aPetal->GetBkgRate();
	      TString   SystemName      = aPetal->GetSystemName();
	      TString   LayerName       = aPetal->GetLayerName();
	      TString   ResolutionModel = aPetal->GetResolutionModel();
	      TString   EfficiencyModel = aPetal->GetEfficiencyModel();
	      
	      GGeoObject* aGeoObject = new GGeoPetal(Name,aGeometry->GetGeometryID(),0,
						     Position,aRot,Thickness,Material,
					             IsSensitive,
					             bottomWidth,topWidth,Height,
						     TVector2(InsensFracUneg,InsensFracUpos),
						     TVector2(InsensFracVneg,InsensFracVpos),
						     global,
					             ResolutionU,ResolutionV,DetEffic,ROtime,
					             BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      aGeoObject->SetXOX0(XOX0);
	      aGeometry->PushGeoElement(aGeoObject);
	    }
	    
	    //Negative shift
	    alpha_tmp = LadderAlpha02 + LadderAlpha*ialpha;
	    cos_alpha_tmp = TMath::Cos(alpha_tmp);
            sin_alpha_tmp = TMath::Sin(alpha_tmp);
	    global->GetRotationMatrix(alpha_tmp,TString("Z"),Rot);
	    Rot = Rot0*Rot;
	    RotAngles       = global->GetRotationAngles(Rot*LadderRot);
	    MyPetalPosition = LadderPosition + LadderRadius*cos_alpha_tmp*UVector + LadderRadius*sin_alpha_tmp*VVector + (-0.5*LadderShift)*WVector;
	    
	    sprintf(name,"%.1f",alpha_tmp/global->GetAngleUnit("deg"));
	    for(int iplane=0;iplane<int(LadderPetalsList.size());iplane++) {
	      GGeoPetal* aPetal = dynamic_cast<GGeoPetal*>(LadderPetalsList[iplane]);
	      
	      TString   Name           = LadderName + TString(", Disk ") + aPetal->GetName() + TString(" negative shift at alpha = ") + TString(name) + TString(" deg");
	      TMatrixD aRot;
	      global->GetGlobalRotationMatrix(RotAngles,aRot);
	      TMatrixD aInvRot;
	      aInvRot.ResizeTo(3,3);
	      aInvRot = aRot;
	      aInvRot.Invert();
	      TVector3  Position = aPetal->GetMainSurface()->GetPosition();
	      global->RotateVector(aInvRot,Position);
	      Position += MyPetalPosition;
	      
	      double    Thickness       = aPetal->GetThickness();
	      TString   Material        = aPetal->GetMaterial();
	      double    XOX0            = aPetal->GetXOX0();
	      double    bottomWidth     = aPetal->GetPetalBaseWidth();
	      double    topWidth        = aPetal->GetPetalTopWidth();
	      double    Height          = aPetal->GetPetalHeight();
	      bool      IsSensitive     = aPetal->GetIsSensitive();
	      double    ResolutionU     = aPetal->GetResolutionU();
	      double    ResolutionV     = aPetal->GetResolutionV();
	      double    ROtime          = aPetal->GetROtime();
	      double    DetEffic        = aPetal->GetDetEfficiency();
	      double    InsensFracUneg  = aPetal->GetUInsensitive().X();
	      double    InsensFracUpos  = aPetal->GetUInsensitive().Y();
	      double    InsensFracVneg  = aPetal->GetVInsensitive().X();
	      double    InsensFracVpos  = aPetal->GetVInsensitive().Y();
	      double    BkgRate         = aPetal->GetBkgRate();
	      TString   SystemName      = aPetal->GetSystemName();
	      TString   LayerName       = aPetal->GetLayerName();
	      TString   ResolutionModel = aPetal->GetResolutionModel();
	      TString   EfficiencyModel = aPetal->GetEfficiencyModel();
	      
	      GGeoObject* aGeoObject = new GGeoPetal(Name,aGeometry->GetGeometryID(),0,
						     Position,aRot,Thickness,Material,
					             IsSensitive,
					             bottomWidth,topWidth,Height,
						     TVector2(InsensFracUneg,InsensFracUpos),
						     TVector2(InsensFracVneg,InsensFracVpos),
						     global,
					             ResolutionU,ResolutionV,DetEffic,ROtime,
					             BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      aGeoObject->SetXOX0(XOX0);
	      aGeometry->PushGeoElement(aGeoObject);
	    }
	    
	  } //End of loop over ladders for this radius

	  //Freeing memory of the template GGeoPlane Objects
	  for(int kkk=0;kkk<int(LadderPetalsList.size());kkk++) delete LadderPetalsList[kkk];
	  LadderPetalsList.clear();
	}
	
      }
      else if(n1 >= 2 && TString(token1[0]) == TString("LadderName")) {
	LadderName = TString("");
        for(int i=1;i<n1;i++) {
	  LadderName += TString(token1[i]);
	  if(i < n1-1) LadderName += TString(" ");
        }
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderRadius")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderRadius     = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderShift")) {
	if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderShift      = atof(token1[1])*global->GetDistanceUnit(token1[2]);
      }
      else if(n1 == 3 && TString(token1[0]) == TString("LadderAlpha0")) {
	if(!global->IsAngleUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	LadderAlpha0     = atof(token1[1])*global->GetAngleUnit(token1[2]);
      }
      else if(n1 == 2 && TString(token1[0]) == TString("LadderN"))           LadderN          = atoi(token1[1]);
      else if(n1 == 5 && TString(token1[0]) == TString("LadderPosition")) {
	if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit    = global->GetDistanceUnit(token1[4]);
	LadderPosition = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 5 && TString(token1[0]) == TString("LadderRotAngles")) {
	if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit     = global->GetAngleUnit(token1[4]);
	LadderRotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
      }
      else if(n1 == 1 && TString(token1[0]) == TString("BeginPetal")) {
	bool IsPetalEnd    = false;
        int  Petal_counter = 0;
        int  Petal_Max     = 25;
    
	TString   Name           = Dummy_name;
        TVector3  Position       = Dummy_vector;
        double    Thickness      = Dummy_value;
	double    XOX0           = Dummy_value;
        TString   Material       = Dummy_name;
        double    bottomWidth    = Dummy_value;
        double    topWidth       = Dummy_value;
        double    Height         = Dummy_value;
        bool      IsSensitive    = false;
        double    Resolution     = Dummy_value;
        double    ResolutionU    = Dummy_value;
        double    ResolutionV    = Dummy_value;
        double    ROtime         = Dummy_value;
        double    DetEffic       = Dummy_value;
        double    InsensFracUneg = 0.0;
        double    InsensFracUpos = 0.0;
        double    InsensFracVneg = 0.0;
        double    InsensFracVpos = 0.0;
        double    BkgRate        = 0.0;
        TString   SystemName("");
	TString   LayerName("");
        TString   ResolutionModel("");
	TString   EfficiencyModel("");

        while(!IsPetalEnd && Petal_counter <= Petal_Max && !fin->eof()) {
          char buf2[MAX_CHARS_PER_LINE];
          fin->getline(buf2, MAX_CHARS_PER_LINE);
	  line_number++;
          int n2 = 0;
          const char* token2[MAX_TOKENS_PER_LINE] = {};
      
          token2[0] = strtok(buf2, DELIMITER); // first token
          if(!token2[0] || TString(token2[0]).BeginsWith(COMMENT)) continue;
    
          for(n2=1; n2<MAX_TOKENS_PER_LINE; n2++) {
            token2[n2] = strtok(0, DELIMITER);
	    if(!token2[n2] || TString(token2[n2]).BeginsWith(COMMENT)) break; // no more tokens
          }
          
          Ladder_counter++;
	  Petal_counter++;
          
          if(n2 == 1 && TString(token2[0]) == TString("EndPetal")) {
	    bool GoodGeoObject = true;

	    if(Name == Dummy_name) {
	      cout << "ERROR in Petal Block inside MosaicLadderPetal: Name not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Position == Dummy_vector) {
	      cout << "ERROR in Petal Block inside MosaicLadderPetal: Position not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Thickness == Dummy_value) {
	      cout << "ERROR in Petal Block  inside MosaicLadderPetal: Thickness not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    if(Material == Dummy_name) {
	      cout << "ERROR in Petal Block  inside MosaicLadderPetal: Material not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(bottomWidth == Dummy_value) {
	      cout << "ERROR in Petal Block  inside MosaicLadderPetal: bottomWidth not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(topWidth == Dummy_value) {
	      cout << "ERROR in Petal Block  inside MosaicLadderPetal: topWidth not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    else if(Height == Dummy_value) {
	      cout << "ERROR in Petal Block  inside MosaicLadderPetal: Height not specified!!!" << endl;
	      GoodGeoObject = false;
	    }
	    
	    if(LadderRotAngles == Dummy_vector) {
	      cout << "WARNING in Petal Block  inside MosaicLadderPetal: RotAngles not specified. Setting all to zero." << endl;
	      LadderRotAngles = TVector3(0,0,0);
	    }
	  
	    if(IsSensitive) {
	      if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	        if(Resolution == Dummy_value) {
		  cout << "ERROR in Petal Block  inside MosaicLadderPetal: Resolution not specified!!!" << endl;
		  GoodGeoObject = false;
	        }
	        else {
		  ResolutionU = Resolution;
		  ResolutionV = Resolution;
	        }
	      }
	      else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	        cout << "ERROR in Petal Block  inside MosaicLadderPetal: ResolutionU not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	        cout << "ERROR in Petal Block  inside MosaicLadderPetal: ResolutionV not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(ROtime == Dummy_value) {
	        cout << "ERROR in Petal Block  inside MosaicLadderPetal: ROtime not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	      else if(DetEffic == Dummy_value) {
	        cout << "ERROR in Petal Block  inside MosaicLadderPetal: DetEffic not specified!!!" << endl;
	        GoodGeoObject = false;
	      }
	    }
	  
	    if(GoodGeoObject) {
	      TMatrixD Rot;
	      global->GetGlobalRotationMatrix(LadderRotAngles,Rot);
	      
	      TVector3  GlobalPosition = Position;
	      
	      GGeoObject* aGeoObject = new GGeoPetal(Name,aGeometry->GetGeometryID(),0,
						     GlobalPosition,Rot,Thickness,Material,
					             IsSensitive,
					             bottomWidth,topWidth,Height,
						     TVector2(InsensFracUneg,InsensFracUpos),
						     TVector2(InsensFracVneg,InsensFracVpos),
						     global,
					             ResolutionU,ResolutionV,DetEffic,ROtime,
					             BkgRate);
	      aGeoObject->SetSystemName(SystemName);
	      aGeoObject->SetLayerName(LayerName);
	      aGeoObject->SetResolutionModel(ResolutionModel);
	      aGeoObject->SetEfficiencyModel(EfficiencyModel);
	      if(XOX0 != Dummy_value) aGeoObject->SetXOX0(XOX0);
	      LadderPetalsList.push_back(aGeoObject);
	      
	      IsPetalEnd = true;
	    }

          }
          else if(n2 >= 2 && TString(token2[0]) == TString("Name")) {
	    Name = TString("");
            for(int i=1;i<n2;i++) {
	      Name += TString(token2[i]);
	      if(i < n2-1) Name += TString(" ");
	    }
	  }
          else if(n2 == 3 && TString(token2[0]) == TString("Thickness")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Thickness       = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("XOX0"))            XOX0            = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("Material"))        Material        = TString(token2[1]);
          else if(n2 == 3 && TString(token2[0]) == TString("bottomWidth")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    bottomWidth     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("topWidth")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    topWidth        = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("Height")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Height          = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("Resolution")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    Resolution      = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionU")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionU     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ResolutionV")) {
	    if(!global->IsDistanceUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ResolutionV     = atof(token2[1])*global->GetDistanceUnit(TString(token2[2]));
	  }
	  else if(n2 == 3 && TString(token2[0]) == TString("ROtime")) {
	    if(!global->IsTimeUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    ROtime          = atof(token2[1])*global->GetTimeUnit(TString(token2[2]));
	  }
	  else if(n2 == 2 && TString(token2[0]) == TString("Efficiency"))      DetEffic        = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracUneg"))  InsensFracUneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracUpos"))  InsensFracUpos  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVneg"))  InsensFracVneg  = atof(token2[1]);
	  else if(n2 == 2 && TString(token2[0]) == TString("InsensFracVpos"))  InsensFracVpos  = atof(token2[1]);
	  else if(n2 == 3 && TString(token2[0]) == TString("BkgRate")) {
	    if(!global->IsRateDensityUnit(token2[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    BkgRate         = atof(token2[1])*global->GetRateDensityUnit(TString(token2[2]));
	  }
	  else if(n2 == 5 && TString(token2[0]) == TString("Position")) {
	    if(!global->IsDistanceUnit(token2[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(token2[4]);
	    Position = TVector3(atof(token2[1])*unit,atof(token2[2])*unit,atof(token2[3])*unit);
          }
	  else if(n2 >= 2 && TString(token2[0]) == TString("SystemName")) {
	    SystemName = TString("");
            for(int i=1;i<n2;i++) {
	      SystemName += TString(token2[i]);
	      if(i < n2-1) SystemName += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("LayerName")) {
	    if(n2 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
	    LayerName = TString(token2[1]);
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("ResolutionModel")) {
	    ResolutionModel = TString("");
            for(int i=1;i<n2;i++) {
	      ResolutionModel += TString(token2[i]);
	      if(i < n2-1) ResolutionModel += TString(" ");
            }
	  }
	  else if(n2 >= 2 && TString(token2[0]) == TString("EfficiencyModel")) {
	    EfficiencyModel = TString("");
            for(int i=1;i<n2;i++) {
	      EfficiencyModel += TString(token2[i]);
	      if(i < n2-1) EfficiencyModel += TString(" ");
            }
	  }
          else if(n2 == 2 && TString(token2[0]) == TString("IsSensitive")) {
	    IsSensitive = global->SetBoolFromString(TString(token2[1])); 
	  }
	}
      
        if(!IsPetalEnd) {
          cout << endl;
	  if(Ladder_counter > Ladder_Max) cout << "IMPORTANT: Reached Maximum number of lines (" << Ladder_Max << ") to read a LadderPlane!!!" << endl;
          cout << "Wrong specification of GeoPetal parameters of MosaicLadderPetal object inside the BeginPetal and EndPetal block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginPetal" << endl;
	  cout << "  Name         string (Mandatory)" << endl;
	  cout << "  Position        x       y       z    units (Mandatory)" << endl;
	  cout << "  Thickness       val units (Mandatory)" << endl;
	  cout << "  Material     string (Mandatory)" << endl;
	  cout << "  XOX0            val   (Optional)" << endl;
	  cout << "  bottomWidth     val units (Mandatory)" << endl;
	  cout << "  topWidth        val units (Mandatory)" << endl;
	  cout << "  Height          val units (Mandatory)" << endl;
	  cout << "  IsSensitive     bool (Mandatory)" << endl;
	  cout << "  Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  ROtime          val units (Mandatory if IsSensitive = true)" << endl;
	  cout << "  Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
	  cout << "  InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
	  cout << "  InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
	  cout << "  InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
	  cout << "  InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
	  cout << "  BkgRate         val units (Optional, default is zero)" << endl;
	  cout << "  SystemName      string (Optional)" << endl;
	  cout << "  LayerName       string (Optional)" << endl;
	  cout << "  ResolutionModel string (Optional)" << endl;
	  cout << "  EfficiencyModel string (Optional)" << endl;
	  cout << "EndPetal" << endl;
          cout << "Check you inputs. Exiting now!!!" << endl;
          cout << endl;
	  
          assert(false);
	  
        } // End While Petal

      } // end if BeginPetal

    } // End While MosaicLadderPetal
    
    if(!IsLadderEnd) {
      cout << endl;
      cout << "Wrong specification of MosaicLadderPetal parameters inside the BeginMosaicLadderPetal and EndMosaicLadderPetal block. Parameters have to be specified as follows:" << endl;
      cout << "BeginMosaicLadderPetal" << endl;
      cout << "  LadderName        string (Mandatory)"           << endl;
      cout << "  LadderPosition      x       y       z    units (Mandatory)" << endl;
      cout << "  LadderRotAngles  alphaX  alphaY  alphaZ  units (Optional)" << endl;
      cout << "  LadderRadius       val   units  (Mandatory)" << endl;
      cout << "  LadderAlpha0       val   units  (Optional, default value is zero)" << endl;
      cout << "  LadderShift        val   units  (Mandatory)" << endl;
      cout << "  LadderN            N  (Mandatory and > 0)" << endl;
      cout << "  BeginPetal" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Position        x       y       z    units (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    XOX0            val   (Optional)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    bottomWidth     val units (Mandatory)" << endl;
      cout << "    topWidth        val units (Mandatory)" << endl;
      cout << "    Height          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
      cout << "    InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndPetal" << endl;
      cout << "  ...  " << endl;
      cout << "  BeginPetal" << endl;
      cout << "    Name         string (Mandatory)" << endl;
      cout << "    Position        x       y       z    units (Mandatory)" << endl;
      cout << "    Thickness       val units (Mandatory)" << endl;
      cout << "    Material     string (Mandatory)" << endl;
      cout << "    XOX0            val   (Optional)" << endl;
      cout << "    bottomWidth     val units (Mandatory)" << endl;
      cout << "    topWidth        val units (Mandatory)" << endl;
      cout << "    Height          val units (Mandatory)" << endl;
      cout << "    IsSensitive     bool (Mandatory)" << endl;
      cout << "    Resolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    ROtime          val units (Mandatory if IsSensitive = true)" << endl;
      cout << "    Efficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
      cout << "    InsensFracUneg  val // insensitive length fraction for low-edge  U (Optional, default is zero)" << endl;
      cout << "    InsensFracUpos  val // insensitive length fraction for high-edge U (Optional, default is zero)" << endl;
      cout << "    InsensFracVneg  val // insensitive length fraction for low-edge  V (Optional, default is zero)" << endl;
      cout << "    InsensFracVpos  val // insensitive length fraction for high-edge V (Optional, default is zero)" << endl;
      cout << "    BkgRate         val units (Optional, default is zero)" << endl;
      cout << "    SystemName      string (Optional)" << endl;
      cout << "    LayerName       string (Optional)" << endl;
      cout << "    ResolutionModel string (Optional)" << endl;
      cout << "    EfficiencyModel string (Optional)" << endl;
      cout << "  EndPetal" << endl;
      cout << "EndMosaicLadderPetal"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }    
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::ReadSystemsConfiguration(ifstream* fin,GGeometry* aGeometry,const char* datacard, int &line_number)
{
  
  bool IsSystemsEnd    = false;
  int  Systems_counter = 0;
  int  Systems_Max     = 1000;
    
  while(!IsSystemsEnd && Systems_counter <= Systems_Max && !fin->eof()) {
    char buf[MAX_CHARS_PER_LINE];
    fin->getline(buf, MAX_CHARS_PER_LINE);
    line_number++;
    int n = 0;
    const char* token[MAX_TOKENS_PER_LINE] = {};
      
    token[0] = strtok(buf, DELIMITER); // first token
    if(!token[0] || TString(token[0]).BeginsWith(COMMENT)) continue;
    
    for(n=1; n<MAX_TOKENS_PER_LINE; n++) {
      token[n] = strtok(0, DELIMITER);
      if(!token[n] || TString(token[n]).BeginsWith(COMMENT)) break; // no more tokens
    }
      
    Systems_counter++;
      
    if(n == 1 && TString(token[0]) == TString("EndSystemsConfiguration")) {
      IsSystemsEnd = true;
    }
    else if(n == 1 && TString(token[0]) == TString("BeginSystem")) {
      bool IsEnd    = false;
      int  counter = 0;
      int  Max     = 5;
    
      System_t aSystem;
      aSystem.Name = Dummy_name;
      aSystem.LayersList.clear();

      while(!IsEnd && counter <= Max && !fin->eof()) {
        char buf1[MAX_CHARS_PER_LINE];
        fin->getline(buf1, MAX_CHARS_PER_LINE);
	line_number++;
        int n1 = 0;
        const char* token1[MAX_TOKENS_PER_LINE] = {};
      
        token1[0] = strtok(buf1, DELIMITER); // first token
        if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
        for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
          token1[n1] = strtok(0, DELIMITER);
	  if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
        }
          
        Systems_counter++;
	counter++;
          
        if(n1 == 1 && TString(token1[0]) == TString("EndSystem")) {
	  bool GoodSystem = true;
	    
	  if(aSystem.Name == Dummy_name) {
	    cout << "ERROR in System Block: Name not specified!!!" << endl;
	    GoodSystem = false;
	  }
	  else if(aSystem.LayersList.size() == 0) {
	    cout << "ERROR in System Block: Layer list not specified!!!" << endl;
	    GoodSystem = false;
	  }
	  
	  if(GoodSystem) {
	    aGeometry->PushASystemIntoGeometry(aSystem);
	    IsEnd = true;
	  }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("Name")) {
	  aSystem.Name = TString("");
	  for(int i=1;i<n1;i++) {
	    aSystem.Name += TString(token1[i]);
	    if(i < n1-1) aSystem.Name += TString(" ");
	  }
	}
	else if(n1 >= 2 && TString(token1[0]) == TString("LayersList")) {
	  for(int i=1;i<n1;i++) aSystem.LayersList.push_back(TString(token1[i]));
	}
      } //end of System while
      
      if(!IsEnd) {
	cout << endl;
        cout << "Wrong specification of System parameters inside the BeginSystem and EndSystem block. Parameters have to be specified as follows:" << endl;
	cout << "BeginSystem" << endl;
	cout << "  Name         string (Mandatory)" << endl;
	cout << "  LayersList   string1  string2  string3  (Mandatory)" << endl;
	cout << "EndSystem" << endl;
        cout << "Check you inputs. Exiting now!!!" << endl;
        cout << endl;
	  
        assert(false);
      }
	
    } // end if BeginSystem
    
  } //end of Systems configuration while
  
  if(!IsSystemsEnd) {
    cout << endl;
    if(Systems_counter > Systems_Max) cout << "IMPORTANT: Reached Maximum number of lines (" << Systems_Max << ") to read a LadderPlane!!!" << endl;
    cout << "Wrong specification of system configurtion parameters inside the BeginSystemsConfiguration and EndSystemsConfiguration block. Parameters have to be specified as follows:" << endl;
    cout << "BeginSystemsConfiguration" << endl;
    cout << "  BeginSystem" << endl;
    cout << "    Name         string (Mandatory)" << endl;
    cout << "    LayersList   string1  string2  string3  (Mandatory)" << endl;
    cout << "  EndSystem" << endl;
    cout << "  ... " << endl;
    cout << "  BeginSystem" << endl;
    cout << "    Name         string (Mandatory)" << endl;
    cout << "    LayersList   string1  string2  string3  (Mandatory)" << endl;
    cout << "  EndSystem" << endl;
    cout << "EndSystemsConfiguration" << endl;
    cout << "Check you inputs. Exiting now!!!" << endl;
    cout << endl;
      
    assert(false);
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::ReadGasDetector(ifstream* fin,GGeometry* aGeometry, const char* datacard, int &line_number)
{
  
  bool IsEnd  = false;
  int counter = 0;
  int Max     = 25;
  
  TString   Name           = Dummy_name;
  TVector3  Position       = Dummy_vector;
  TVector3  RotAngles      = Dummy_vector;
  TString   Material       = Dummy_name;
  double    Rin            = Dummy_value;
  double    Rout           = Dummy_value;
  double    Length         = Dummy_value;
  int       NLayers        = Dummy_value;
  double    Resolution     = Dummy_value;
  double    ResolutionU    = Dummy_value;
  double    ResolutionV    = Dummy_value;
  double    ROtime         = Dummy_value;
  double    DetEffic       = Dummy_value;
  double    BkgRate        = 0.0;
  TString   SystemName("");
  TString   LayerName("");
  TString   ResolutionModel("");
  TString   EfficiencyModel("");
  
  while(!IsEnd && counter <= Max && !fin->eof()) {
    char buf1[MAX_CHARS_PER_LINE];
    fin->getline(buf1, MAX_CHARS_PER_LINE);
    line_number++;
    int n1 = 0;
    const char* token1[MAX_TOKENS_PER_LINE] = {};
    
    token1[0] = strtok(buf1, DELIMITER); // first token
    if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
    for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
      token1[n1] = strtok(0, DELIMITER);
      if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
    }
    
    if(n1 == 1 && TString(token1[0]) == TString("EndGasDetector")) {
      bool GoodGeoObject = true;
      
      if(Name == Dummy_name) {
	cout << "ERROR in GasDetector Block: Name not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(Position == Dummy_vector) {
	cout << "ERROR in GasDetector Block: Position not specified!!!" << endl;
	GoodGeoObject = false;
      }
      if(Material == Dummy_name) {
	cout << "ERROR in GasDetector Block: Material not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(Rin == Dummy_value) {
	cout << "ERROR in GasDetector Block: Rin not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(Rout == Dummy_value) {
	cout << "ERROR in GasDetector Block: Rout not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(Rout <= Rin) {
	cout << "ERROR in GasDetector Block: Rout is smaller than Rin!!!" << endl;
	GoodGeoObject = false;
      }
      else if(Length == Dummy_value) {
	cout << "ERROR in GasDetector Block: Length not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(NLayers == int(Dummy_value)) {
	cout << "ERROR in GasDetector Block: NLayers not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(NLayers <= 0) {
	cout << "ERROR in GasDetector Block: NLayers is <= 0 !!!" << endl;
	GoodGeoObject = false;
      }
      
      if(RotAngles == Dummy_vector) {
	cout << "WARNING in GasDetector Block: RotAngles not specified. Setting all to zero." << endl;
	RotAngles = TVector3(0,0,0);
      }
      
      if(ResolutionU == Dummy_value && ResolutionV == Dummy_value) {
	if(Resolution == Dummy_value) {
	  cout << "ERROR in GasDetector Block: Resolution not specified!!!" << endl;
	  GoodGeoObject = false;
	}
	else {
	  ResolutionU = Resolution;
	  ResolutionV = Resolution;
	}
      }
      else if(ResolutionU == Dummy_value && ResolutionV != Dummy_value) {
	cout << "ERROR in GasDetector Block: ResolutionU not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(ResolutionU != Dummy_value && ResolutionV == Dummy_value) {
	cout << "ERROR in GasDetector Block: ResolutionV not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(ROtime == Dummy_value) {
	cout << "ERROR in GasDetector Block: ROtime not specified!!!" << endl;
	GoodGeoObject = false;
      }
      else if(DetEffic == Dummy_value) {
	cout << "ERROR in GasDetector Block: DetEffic not specified!!!" << endl;
	GoodGeoObject = false;
      }
      
      if(GoodGeoObject) {
	TMatrixD Rot;
	global->GetGlobalRotationMatrix(RotAngles,Rot);
	
	double Thickness = (Rout - Rin)/NLayers;
	
	for(int ilayer=0;ilayer<NLayers;ilayer++) {
	  double R_tmp     = Rin + Thickness*(0.5 + ilayer);
	  TString Name_tmp = Name + TString(" ") + long(ilayer+1);
	  
	  GGeoObject* aGeoObject = new GGeoCylinder(Name_tmp,
						    aGeometry->GetGeometryID(),0,
						    Position,Rot,Thickness,Material,
					            true,
					            Length,R_tmp,
					            TVector2(0,0),
						    global,
					            ResolutionU,ResolutionV,DetEffic,ROtime,
					            BkgRate);
	  aGeoObject->SetSystemName(SystemName);
	  aGeoObject->SetLayerName(LayerName);
	  aGeoObject->SetResolutionModel(ResolutionModel);
	  aGeoObject->SetEfficiencyModel(EfficiencyModel);
	  aGeometry->PushGeoElement(aGeoObject);
	}
	
	IsEnd = true;
      }
      
    }
    else if(n1 >= 2 && TString(token1[0]) == TString("GasDetName")) {
      Name = TString("");
      for(int i=1;i<n1;i++) {
	Name += TString(token1[i]);
	if(i < n1-1) Name += TString(" ");
      } 
    }
    else if(n1 == 2 && TString(token1[0]) == TString("GasDetMaterial"))        Material        = TString(token1[1]);
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetRin")) {
      if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      Rin          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
    }
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetRout")) {
      if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      Rout          = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
    }
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetLength")) {
      if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      Length         = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
    }
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetResolution")) {
      if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      Resolution      = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
    }
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetResolutionU")) {
      if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      ResolutionU     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
    }
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetResolutionV")) {
      if(!global->IsDistanceUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      ResolutionV     = atof(token1[1])*global->GetDistanceUnit(TString(token1[2]));
    }
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetROtime")) {
      if(!global->IsTimeUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      ROtime          = atof(token1[1])*global->GetTimeUnit(TString(token1[2]));
    }
    else if(n1 == 2 && TString(token1[0]) == TString("GasDetEfficiency"))      DetEffic        = atof(token1[1]);
    else if(n1 == 2 && TString(token1[0]) == TString("GasDetNLayers"))         NLayers         = atoi(token1[1]);
    else if(n1 == 3 && TString(token1[0]) == TString("GasDetBkgRate")) {
      if(!global->IsRateDensityUnit(token1[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      BkgRate         = atof(token1[1])*global->GetRateDensityUnit(TString(token1[2]));
    }
    else if(n1 >= 2 && TString(token1[0]) == TString("GasDetSystemName")) {
      SystemName = TString("");
      for(int i=1;i<n1;i++) {
	SystemName += TString(token1[i]);
	if(i < n1-1) SystemName += TString(" ");
      }
    }
    else if(n1 >= 2 && TString(token1[0]) == TString("GasDetLayerName")) {
      if(n1 > 2) cout << "WARNING in specification of LayerName, file name " << datacard << ", line-number " << line_number << ". More than one string specified. Taking the first one." << endl;
      LayerName = TString(token1[1]);
    }
    else if(n1 >= 2 && TString(token1[0]) == TString("GasDetResolutionModel")) {
      ResolutionModel = TString("");
      for(int i=1;i<n1;i++) {
	ResolutionModel += TString(token1[i]);
	if(i < n1-1) ResolutionModel += TString(" ");
      }
    }
    else if(n1 >= 2 && TString(token1[0]) == TString("GasDetEfficiencyModel")) {
      EfficiencyModel = TString("");
      for(int i=1;i<n1;i++) {
	EfficiencyModel += TString(token1[i]);
	if(i < n1-1) EfficiencyModel += TString(" ");
      }
    }
    else if(n1 == 5 && TString(token1[0]) == TString("GasDetPosition")) {
      if(!global->IsDistanceUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      double unit = global->GetDistanceUnit(token1[4]);
      Position = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
    }
    else if(n1 == 5 && TString(token1[0]) == TString("GasDetRotAngles")) {
      if(!global->IsAngleUnit(token1[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
      double unit = global->GetAngleUnit(token1[4]);
      RotAngles = TVector3(atof(token1[1])*unit,atof(token1[2])*unit,atof(token1[3])*unit);
    }
    
    counter++;
  }
  
  if(!IsEnd) {
    cout << endl;
    cout << "Wrong specification of GeoPlane parameters inside the BeginGasDetector and EndGasDetector block. Parameters have to be specified as follows:" << endl;
    cout << "BeginGasDetector" << endl;
    cout << "  GasDetName            string (Mandatory)" << endl;
    cout << "  GasDetPosition           x       y       z    units (Mandatory)" << endl;
    cout << "  GasDetRotAngles       alphaX  alphaY  alphaZ  units (Optional, default is zero vector)" << endl;
    cout << "  GasDetMaterial        string (Mandatory)" << endl;
    cout << "  GasDetRin             val units (Mandatory)" << endl;
    cout << "  GasDetRout            val units (Mandatory)" << endl;
    cout << "  GasDetLength          val units (Mandatory)" << endl;
    cout << "  GasDetResolution      val units // If ResolutionU/V not specified will use same resolution for U & V directions (Mandatory if IsSensitive = true)" << endl;
    cout << "  GasDetResolutionU     val units (Mandatory if IsSensitive = true)" << endl;
    cout << "  GasDetResolutionV     val units (Mandatory if IsSensitive = true)" << endl;
    cout << "  GasDetROtime          val units (Mandatory if IsSensitive = true)" << endl;
    cout << "  GasDetEfficiency      val // Sensor Det-efficiency (Mandatory if IsSensitive = true)" << endl;
    cout << "  GasDetBkgRate         val units (Optional, default is zero)" << endl;
    cout << "  GasDetSystemName      string (Optional)" << endl;
    cout << "  GasDetLayerName       string (Optional)" << endl;
    cout << "  GasDetResolutionModel string (Optional)" << endl;
    cout << "  GasDetEfficiencyModel string (Optional)" << endl;
    cout << "EndGasDetector" << endl;
    cout << "Check you inputs. Exiting now!!!" << endl;
    cout << endl;
    
    assert(false);
  }
  
  return;
  
}
//====================================================================
void  Guariguanchi::ReadTrackFinderAlgo(TString  BlockName,ifstream* fin,GGeometry* aGeometry,const char* datacard, int &line_number)
{
  
  
  if(BlockName == TString("BeginFPCCDTrackFinderAlgo")) {
    bool IsEnd  = false;
    int counter = 0;
    int Max     = 10000;
    
    std::vector<TString>       FullSystemList;
    std::vector<SeedLayers_t>  SeedLayerConfigs;
    std::vector<TrackFindingRegion_t>  Regions;
    TString  Name               = Dummy_name;
    int      Nhits_min          = Dummy_value;
    double   PtMin_cut          = Dummy_value;
    double   PurityMin_cut      = 1.0;
    double   Chi2Ondf_Seed_cut  = 3.0;
    double   Chi2Ondf_Add_cut   = 3.0;
    bool     InwardTracking     = true;
    int      NfakesMax_seeding  = 0;
    int      Nmc_seed_effic     = Dummy_value;
    TVector3 CenterPosition     = Dummy_vector;
    
    FullSystemList.clear();
    SeedLayerConfigs.clear();
    Regions.clear();
    
    while(!IsEnd && counter <= Max && !fin->eof()) {
      char buf[MAX_CHARS_PER_LINE];
      fin->getline(buf, MAX_CHARS_PER_LINE);
      line_number++;
      int n = 0;
      const char* token[MAX_TOKENS_PER_LINE] = {};
      
      token[0] = strtok(buf, DELIMITER); // first token
      if(!token[0] || TString(token[0]).BeginsWith(COMMENT)) continue;
    
      for(n=1; n<MAX_TOKENS_PER_LINE; n++) {
        token[n] = strtok(0, DELIMITER);
	if(!token[n] || TString(token[n]).BeginsWith(COMMENT)) break; // no more tokens
      }
      
      if(n == 1 && TString(token[0]) == TString("EndFPCCDTrackFinderAlgo")) {
	bool GoodTrackFinderAlgo = true;
	
	if(FullSystemList.size() == 0) {
	  cout << endl;
	  cout << "ERROR inside FPCCDTrackFinderAlgo block: Systems list has zero elements!!!" << endl;
	  cout << endl;
	  GoodTrackFinderAlgo = false;
	}
	else if(SeedLayerConfigs.size() == 0) {
	  cout << endl;
	  cout << "ERROR inside FPCCDTrackFinderAlgo block: Seed configurations list has zero elements!!!" << endl;
	  cout << endl;
	  GoodTrackFinderAlgo = false;
	}
	else if(Regions.size() == 0) {
	  cout << endl;
	  cout << "ERROR inside FPCCDTrackFinderAlgo block: No track finding regions has been specified!!!" << endl;
	  cout << endl;
	  GoodTrackFinderAlgo = false;
	}
	else if(Name == Dummy_name) {
	  cout << endl;
	  cout << "ERROR inside FPCCDTrackFinderAlgo block: No track finding algorithm name has been specified!!!" << endl;
	  cout << endl;
	  GoodTrackFinderAlgo = false;
	}
	else if(Nhits_min == Dummy_value) {
	  cout << endl;
	  cout << "ERROR inside FPCCDTrackFinderAlgo block: Nhits minimum has been specified!!!" << endl;
	  cout << endl;
	  GoodTrackFinderAlgo = false;
	}
	else if(PtMin_cut == Dummy_value) {
	  cout << endl;
	  cout << "ERROR inside FPCCDTrackFinderAlgo block: minimum Pt cut has been specified!!!" << endl;
	  cout << endl;
	  GoodTrackFinderAlgo = false;
	}
	
	if(CenterPosition == Dummy_vector) {
	  CenterPosition = TVector3(0.0,0.0,0.0);
	}
	
	if(GoodTrackFinderAlgo) {
	  aGeometry->PushTrackFinderAlgoIntoGeometry(new  GTrackFinderAlgoFPCCD(Name,
										FullSystemList,
									        SeedLayerConfigs,
									        Regions[0],
									        global,
									        Nhits_min,
									        PtMin_cut,
									        PurityMin_cut,
									        Chi2Ondf_Seed_cut,
									        Chi2Ondf_Add_cut,
									        InwardTracking,
									        NfakesMax_seeding,
										CenterPosition));
	  if(Nmc_seed_effic != Dummy_value) {
	    int algo_index = aGeometry->GetNTrackFinderAlgos() - 1;
	    GTrackFinderAlgoFPCCD* anAlgo = dynamic_cast<GTrackFinderAlgoFPCCD*>(aGeometry->GetTrackFinderAlgo(algo_index));
	    anAlgo->SetNmcSeedEffic(Nmc_seed_effic);
	  }
	  
	  IsEnd = true;
	}
	
      }
      else if(n >= 2 && TString(token[0]) == TString("Name")) {
        Name = TString("");
        for(int i=1;i<n;i++) {
	  Name += TString(token[i]);
	  if(i < n-1) Name += TString(" ");
        }
      }
      else if(n == 3 && TString(token[0]) == TString("PtMin")) {
	if(!global->IsMomentumUnit(token[2])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	PtMin_cut = atof(token[1])*global->GetMomentumUnit(TString(token[2]));
      }
      else if(n == 2 && TString(token[0]) == TString("NhitsMin"))          Nhits_min      = atoi(token[1]);
      else if(n == 2 && TString(token[0]) == TString("NmcSeedEffic"))      Nmc_seed_effic = atoi(token[1]);
      else if(n == 2 && TString(token[0]) == TString("PurityMin"))         PurityMin_cut      = atof(token[1]);
      else if(n == 2 && TString(token[0]) == TString("Chi2OndfSeed"))      Chi2Ondf_Seed_cut  = atof(token[1]);
      else if(n == 2 && TString(token[0]) == TString("Chi2OndfAdd"))       Chi2Ondf_Add_cut   = atof(token[1]);
      else if(n == 2 && TString(token[0]) == TString("InwardTracking"))    InwardTracking     = global->SetBoolFromString(TString(token[1]));
      else if(n == 2 && TString(token[0]) == TString("NfakesMaxSeeding"))  NfakesMax_seeding  = atoi(token[1]);
      else if(n == 5 && TString(token[0]) == TString("CenterPosition")) {
	if(!global->IsDistanceUnit(token[4])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	double unit = global->GetDistanceUnit(TString(token[4]));
	CenterPosition = unit*TVector3(atof(token[1]),atof(token[2]),atof(token[3]));
      }
      else if(n == 1 && TString(token[0]) == TString("BeginSystems") && FullSystemList.size() == 0)  {
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
                    
          if(n1 == 1 && TString(token1[0]) == TString("EndSystems")) {
	    bool Good = true;
	    
	    if(FullSystemList.size() == 0) {
	      cout << endl;
	      cout << "  ERROR in Systems block: Systems list has zero elements!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(Good) IsEnd1 = true;
	  }
	  else if(n1 >= 2 && TString(token1[0]) == TString("System"))  {
	    TString SystemName = TString("");
            for(int i=1;i<n1;i++) {
	      SystemName += TString(token1[i]);
	      if(i < n1-1) SystemName += TString(" ");
	    }
	    if(SystemName != TString("")) FullSystemList.push_back(SystemName);
	  }
	  
	  counter1++;
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of FPCCDTrackFinderAlgo parameters inside the BeginSystems and EndSystems block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginSystems" << endl;
	  cout << "  System   string" << endl;
	  cout << "  ..." << endl;
	  cout << "EndSystems"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
	
      }
      else if(n == 1 && TString(token[0]) == TString("BeginSeedConfigs") && SeedLayerConfigs.size() == 0)  {
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 100;
      
        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
                    
          if(n1 == 1 && TString(token1[0]) == TString("EndSeedConfigs")) {
	    bool Good = true;
	    
	    if(SeedLayerConfigs.size() == 0) {
	      cout << endl;
	      cout << "  ERROR in SeedConfigs block: Seeding configurations list has zero elements!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(Good) IsEnd1 = true;
	  }
	  else if(n1 >= 2 && TString(token1[0]) == TString("SeedConfig"))  {
	    SeedLayers_t aSeedConfig;
	    aSeedConfig.SeedLayers.clear();
	    for(int i=1;i<n1;i++) aSeedConfig.SeedLayers.push_back(TString(token1[i]));
	    
	    if(aSeedConfig.SeedLayers.size() > 0) SeedLayerConfigs.push_back(aSeedConfig);
	  }
	  
	  counter1++;
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of FPCCDTrackFinderAlgo parameters inside the BeginSeedConfigs and EndSeedConfigs block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginSeedConfigs" << endl;
	  cout << "  SeedConfig   string1  string2  string3 ..." << endl;
	  cout << "  ..." << endl;
	  cout << "EndSeedConfigs"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
	
      }
      else if(n == 1 && TString(token[0]) == TString("BeginTrackFinderRegion") && Regions.size() == 0)  {
	bool IsEnd1  = false;
        int counter1 = 0;
        int Max1     = 200;
	
	TrackFindingRegion_t  aRegion;
	aRegion.posRange.Rx[0]     = Dummy_value;  aRegion.posRange.Rx[1]     = Dummy_value;
        aRegion.posRange.Ry[0]     = Dummy_value;  aRegion.posRange.Ry[1]     = Dummy_value;
        aRegion.posRange.Rz[0]     = Dummy_value;  aRegion.posRange.Rz[1]     = Dummy_value;
        aRegion.posRange.Rr[0]     = Dummy_value;  aRegion.posRange.Rr[1]     = Dummy_value;
        aRegion.posRange.Rtheta[0] = Dummy_value;  aRegion.posRange.Rtheta[1] = Dummy_value;
        aRegion.posRange.Rphi[0]   = Dummy_value;  aRegion.posRange.Rphi[1]   = Dummy_value;
        aRegion.momRange.Rp[0]     = Dummy_value;  aRegion.momRange.Rp[1]     = Dummy_value;
        aRegion.momRange.Rtheta[0] = Dummy_value;  aRegion.momRange.Rtheta[1] = Dummy_value;
        aRegion.momRange.Rphi[0]   = Dummy_value;  aRegion.momRange.Rphi[1]   = Dummy_value;

        while(!IsEnd1 && counter1 <= Max1 && !fin->eof()) {
          char buf1[MAX_CHARS_PER_LINE];
          fin->getline(buf1, MAX_CHARS_PER_LINE);
          line_number++;
          int n1 = 0;
          const char* token1[MAX_TOKENS_PER_LINE] = {};
      
          token1[0] = strtok(buf1, DELIMITER); // first token
          if(!token1[0] || TString(token1[0]).BeginsWith(COMMENT)) continue;
    
          for(n1=1; n1<MAX_TOKENS_PER_LINE; n1++) {
            token1[n1] = strtok(0, DELIMITER);
	    if(!token1[n1] || TString(token1[n1]).BeginsWith(COMMENT)) break; // no more tokens
          }
                    
          if(n1 == 1 && TString(token1[0]) == TString("EndTrackFinderRegion")) {
	    bool Good = true;
	    
	    if((aRegion.posRange.Rx[0]     == Dummy_value && aRegion.posRange.Rx[1] == Dummy_value) && 
	       (aRegion.posRange.Rx[0]     == Dummy_value && aRegion.posRange.Ry[1] == Dummy_value) && 
	       (aRegion.posRange.Rx[0]     == Dummy_value && aRegion.posRange.Rz[1] == Dummy_value) && 
	       (aRegion.posRange.Rr[0]     == Dummy_value && aRegion.posRange.Rr[1] == Dummy_value) && 
	       (aRegion.posRange.Rtheta[0] == Dummy_value && aRegion.posRange.Rtheta[1] == Dummy_value) && 
	       (aRegion.posRange.Rphi[0]   == Dummy_value && aRegion.posRange.Rphi[1] == Dummy_value) &&
	       (aRegion.momRange.Rp[0]     == Dummy_value && aRegion.momRange.Rp[1] == Dummy_value) && 
	       (aRegion.momRange.Rtheta[0] == Dummy_value && aRegion.momRange.Rtheta[1] == Dummy_value) && 
	       (aRegion.momRange.Rphi[0]   == Dummy_value && aRegion.momRange.Rphi[1] == Dummy_value)) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: no position or momentum ranges specficied!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(aRegion.posRange.Rx[0]     == Dummy_value) aRegion.posRange.Rx[0]     = -1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Rx[1]     == Dummy_value) aRegion.posRange.Rx[1]     = +1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Ry[0]     == Dummy_value) aRegion.posRange.Ry[0]     = -1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Ry[1]     == Dummy_value) aRegion.posRange.Ry[1]     = +1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Rz[0]     == Dummy_value) aRegion.posRange.Rz[0]     = -1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Rz[1]     == Dummy_value) aRegion.posRange.Rz[1]     = +1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Rr[0]     == Dummy_value) aRegion.posRange.Rr[0]     =     0.0*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Rr[1]     == Dummy_value) aRegion.posRange.Rr[1]     = +1.0e+3*global->GetDistanceUnit("m");
	    if(aRegion.posRange.Rtheta[0] == Dummy_value) aRegion.posRange.Rtheta[0] =     0.0*global->GetAngleUnit("deg");
	    if(aRegion.posRange.Rtheta[1] == Dummy_value) aRegion.posRange.Rtheta[1] =  +180.0*global->GetAngleUnit("deg");
	    if(aRegion.posRange.Rphi[0]   == Dummy_value) aRegion.posRange.Rphi[0]   =  -180.0*global->GetAngleUnit("deg");
	    if(aRegion.posRange.Rphi[1]   == Dummy_value) aRegion.posRange.Rphi[1]   =  +180.0*global->GetAngleUnit("deg");
	    if(aRegion.momRange.Rp[0]     == Dummy_value) aRegion.momRange.Rp[0]     =     0.0*global->GetMomentumUnit("GeV/c");
	    if(aRegion.momRange.Rp[1]     == Dummy_value) aRegion.momRange.Rp[1]     = +1.0e+3*global->GetMomentumUnit("TeV/c");
	    if(aRegion.momRange.Rtheta[0] == Dummy_value) aRegion.momRange.Rtheta[0] =     0.0*global->GetAngleUnit("deg");
	    if(aRegion.momRange.Rtheta[1] == Dummy_value) aRegion.momRange.Rtheta[1] =  +180.0*global->GetAngleUnit("deg");
	    if(aRegion.momRange.Rphi[0]   == Dummy_value) aRegion.momRange.Rphi[0]   =  -180.0*global->GetAngleUnit("deg");
	    if(aRegion.momRange.Rphi[1]   == Dummy_value) aRegion.momRange.Rphi[1]   =  +180.0*global->GetAngleUnit("deg");
	    
	    if(aRegion.posRange.Rx[0] >= aRegion.posRange.Rx[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for position X range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Ry[0] >= aRegion.posRange.Ry[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for position Y range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rz[0] >= aRegion.posRange.Rz[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for position Z range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rr[0] < 0.0) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min < 0 for position R range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rr[1] < 0.0) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: max < 0 for position R range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rr[0] >= aRegion.posRange.Rr[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for position R range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rtheta[0] < 0.0) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min < 0 for position Theta range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rtheta[0] > 180.0*global->GetAngleUnit("deg")) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: max > 180 deg for position Theta range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rtheta[0] >= aRegion.posRange.Rtheta[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for position Theta range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rphi[0] < -180.0*global->GetAngleUnit("deg")) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min < -180 deg for position Phi range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rphi[0] > 360.0*global->GetAngleUnit("deg")) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: max > 360 deg for position Phi range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.posRange.Rphi[0] >= aRegion.posRange.Rphi[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for position Phi range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rp[0] < 0.0) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min < 0 for momentum P range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rp[1] < 0.0) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: max < 0 for momentum P range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rp[0] >= aRegion.momRange.Rp[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for momentum P range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rtheta[0] < 0.0) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min < 0 for momentum Theta range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rtheta[0] > 180.0*global->GetAngleUnit("deg")) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: max > 180 deg for momentum Theta range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rtheta[0] >= aRegion.momRange.Rtheta[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for momentum Theta range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rphi[0] < -180.0*global->GetAngleUnit("deg")) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min < -180 deg for momentum Phi range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rphi[0] > 360.0*global->GetAngleUnit("deg")) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: max > 360 deg for momentum Phi range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    else if(aRegion.momRange.Rphi[0] >= aRegion.momRange.Rphi[1]) {
	      cout << endl;
	      cout << "  ERROR in TrackFinderRegion block: min >= max for momentum Phi range!!!"<< endl;
	      cout << endl;
	      Good = false;
	    }
	    
	    if(Good) {
	      Regions.push_back(aRegion);
	      IsEnd1 = true;
	    }
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("positionRangeX"))  {
	    if(!global->IsDistanceUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(TString(token1[3]));
	    aRegion.posRange.Rx[0] = atof(token1[1])*unit;
	    aRegion.posRange.Rx[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("positionRangeY"))  {
	    if(!global->IsDistanceUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(TString(token1[3]));
	    aRegion.posRange.Ry[0] = atof(token1[1])*unit;
	    aRegion.posRange.Ry[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("positionRangeZ"))  {
	    if(!global->IsDistanceUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(TString(token1[3]));
	    aRegion.posRange.Rz[0] = atof(token1[1])*unit;
	    aRegion.posRange.Rz[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("positionRangeR"))  {
	    if(!global->IsDistanceUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetDistanceUnit(TString(token1[3]));
	    aRegion.posRange.Rr[0] = atof(token1[1])*unit;
	    aRegion.posRange.Rr[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("positionRangeTheta"))  {
	    if(!global->IsAngleUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(TString(token1[3]));
	    aRegion.posRange.Rtheta[0] = atof(token1[1])*unit;
	    aRegion.posRange.Rtheta[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("positionRangePhi"))  {
	    if(!global->IsAngleUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(TString(token1[3]));
	    aRegion.posRange.Rphi[0] = atof(token1[1])*unit;
	    aRegion.posRange.Rphi[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("momentumRangeP"))  {
	    if(!global->IsMomentumUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetMomentumUnit(TString(token1[3]));
	    aRegion.momRange.Rp[0] = atof(token1[1])*unit;
	    aRegion.momRange.Rp[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("momentumRangeTheta"))  {
	    if(!global->IsAngleUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(TString(token1[3]));
	    aRegion.momRange.Rtheta[0] = atof(token1[1])*unit;
	    aRegion.momRange.Rtheta[1] = atof(token1[2])*unit;
	  }
	  else if(n1 >= 4 && TString(token1[0]) == TString("momentumRangePhi"))  {
	    if(!global->IsAngleUnit(token1[3])) cout << "ERROR in DataCard " << datacard << " line number " << line_number << endl;
	    double unit = global->GetAngleUnit(TString(token1[3]));
	    aRegion.momRange.Rphi[0] = atof(token1[1])*unit;
	    aRegion.momRange.Rphi[1] = atof(token1[2])*unit;
	  }
	  
	  counter1++;
	}
	
	if(!IsEnd1) {
	  cout << endl;
          cout << "Wrong specification of FPCCDTrackFinderAlgo parameters inside the BeginTrackFinderRegion and EndTrackFinderRegion block. Parameters have to be specified as follows:" << endl;
	  cout << "BeginTrackFinderRegion" << endl;
	  cout << "  positionRangeX      min  max  units // max > min (Optional)" << endl;
	  cout << "  positionRangeY      min  max  units // max > min (Optional)" << endl;
	  cout << "  positionRangeZ      min  max  units // max > min (Optional)" << endl;
	  cout << "  positionRangeR      min  max  units // max > min and both >= 0 (Optional)" << endl;
	  cout << "  positionRangeTheta  min  max  units // max > min and both in (0,180) deg (Optional)" << endl;
	  cout << "  positionRangePhi    min  max  units // max > min and both in (-180,180) or (0,360) deg (Optional)" << endl;
	  cout << "  momentumRangeP      min  max  units // max > min and both >= 0 (Optional)" << endl;
	  cout << "  momentumRangeTheta  min  max  units // max > min and both in (0,180) deg (Optional)" << endl;
	  cout << "  momentumRangePhi    min  max  units // max > min and both in (-180,180) or (0,360) deg (Optional)" << endl;
	  cout << "EndTrackFinderRegion"   << endl;
	  cout << "Check you inputs. Exiting now!!!" << endl;
	  cout << endl;
	  assert(false);
	}
	
      }
      
      counter++;
    }
      
    if(!IsEnd) {
      cout << endl;
      cout << "Wrong specification of FPCCDTrackFinderAlgo parameters inside the BeginFPCCDTrackFinderAlgo and EndFPCCDTrackFinderAlgo block. Parameters have to be specified as follows:" << endl;
      cout << "BeginFPCCDTrackFinderAlgo" << endl;
      cout << "  Name              string                (Mandatory)" << endl;
      cout << "  PtMin             val units // val > 0  (Mandatory)" << endl;
      cout << "  NhitsMin          val       // val > 0  (Mandatory)" << endl;
      cout << "  PurityMin         val       // val > 0  (Optional, defaul value is 1)" << endl;
      cout << "  Chi2OndfSeed      val       // val > 0  (Optional, defaul value is 3)" << endl;
      cout << "  Chi2OndfAdd       val       // val > 0  (Optional, defaul value is 3)" << endl;
      cout << "  InwardTracking    bool      //          (Optional, defaul value is true)" << endl;
      cout << "  NfakesMaxSeeding  val       // val >= 0 (Optional, defaul value is 0)" << endl;
      cout << "  BeginTrackFinderRegion" << endl;
      cout << "    positionRangeX      min  max  units // max > min (Optional)" << endl;
      cout << "    positionRangeY      min  max  units // max > min (Optional)" << endl;
      cout << "    positionRangeZ      min  max  units // max > min (Optional)" << endl;
      cout << "    positionRangeR      min  max  units // max > min and both >= 0 (Optional)" << endl;
      cout << "    positionRangeTheta  min  max  units // max > min and both in (0,180) deg (Optional)" << endl;
      cout << "    positionRangePhi    min  max  units // max > min and both in (-180,180) or (0,360) deg (Optional)" << endl;
      cout << "    momentumRangeP      min  max  units // max > min and both >= 0 (Optional)" << endl;
      cout << "    momentumRangeTheta  min  max  units // max > min and both in (0,180) deg (Optional)" << endl;
      cout << "    momentumRangePhi    min  max  units // max > min and both in (-180,180) or (0,360) deg (Optional)" << endl;
      cout << "  EndTrackFinderRegion"   << endl;
      cout << "  BeginSystems" << endl;
      cout << "    System   string" << endl;
      cout << "    ..." << endl;
      cout << "  EndSystems"   << endl;
      cout << "  BeginSeedConfigs" << endl;
      cout << "    SeedConfig   string1  string2  string3 ..." << endl;
      cout << "    ..." << endl;
      cout << "  EndSeedConfigs"   << endl;
      cout << "EndFPCCDTrackFinderAlgo"   << endl;
      cout << "Check you inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }
  
  return;
  
}
//====================================================================
