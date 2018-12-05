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
#include "include/GTrackFinderAlgoFPCCD.h"

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
GTrackFinderAlgoFPCCD::GTrackFinderAlgoFPCCD(TString   aName,
					     std::vector<TString>       aFullSystemList,
					     std::vector<SeedLayers_t>  aSeedLayerConfigs,
					     TrackFindingRegion_t       aRegion,
					     GGlobalTools* aglobal,
					     int       aNhits_min,
					     double    aPtMin_cut,
					     double    aPurityMin_cut,
					     double    aChi2Ondf_Seed_cut,
					     double    aChi2Ondf_Add_cut,
					     bool      aInwardTracking,
					     int       aNfakesMax_seeding,
					     TVector3  aCenterPosition)
                                             : GTrackFinderAlgo(aName,
								aRegion,
								aglobal)
{
  
  Nmc_seed_effic = NMC_SeedEffic;
  
  Type  = TString("FPCCDTrackFinder");
  
  FillSystemList(aFullSystemList);
  FillSeedConfigs(aSeedLayerConfigs);
  
  Nhits_min         = aNhits_min;
  PtMin_cut         = aPtMin_cut;
  PurityMin_cut     = aPurityMin_cut;
  Chi2Ondf_Seed_cut = aChi2Ondf_Seed_cut;
  Chi2Ondf_Add_cut  = aChi2Ondf_Add_cut;
  InwardTracking    = aInwardTracking;
  NfakesMax_seeding = aNfakesMax_seeding;
  CenterPosition    = aCenterPosition;
  
}
//====================================================================
GTrackFinderAlgoFPCCD::GTrackFinderAlgoFPCCD(const GTrackFinderAlgoFPCCD& other,TString aName)
                                             : GTrackFinderAlgo(aName,
								other.Region,
								other.global)
{
  
  Type = other.Type;
  Nmc_seed_effic = other.Nmc_seed_effic;
  
  FillSystemList(other.FullSystemList);
  FillSeedConfigs(other.SeedLayerConfigs);
  
  Nhits_min         = other.Nhits_min;
  PtMin_cut         = other.PtMin_cut;
  PurityMin_cut     = other.PurityMin_cut;
  Chi2Ondf_Seed_cut = other.Chi2Ondf_Seed_cut;
  Chi2Ondf_Add_cut  = other.Chi2Ondf_Add_cut;
  InwardTracking    = other.InwardTracking;
  NfakesMax_seeding = other.NfakesMax_seeding;
  CenterPosition    = other.CenterPosition;
  
}
//====================================================================
GTrackFinderAlgoFPCCD::~GTrackFinderAlgoFPCCD() 
{
  
}
//====================================================================
GTrackFinderAlgoFPCCD* GTrackFinderAlgoFPCCD::clone(TString aName) const
{
 
  return new GTrackFinderAlgoFPCCD(*this,aName);
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::FillSystemList(std::vector<TString>  aFullSystemList)
{
  
  //Remove duplicated system-names in input aFullSystemList
  FullSystemList.clear();
  for(int i=0;i<int(aFullSystemList.size());i++) {
    bool IsAlreadyInList = false;
    for(int j=0;j<int(FullSystemList.size());j++) {
      if(FullSystemList[j] == aFullSystemList[i]) {
	IsAlreadyInList = true;
	break;
      }
    }
    
    if(!IsAlreadyInList) FullSystemList.push_back(aFullSystemList[i]);
  }
  
  if(FullSystemList.size() == 0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::FillSystemList:: FullSystemList has zero entries. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::FillSeedConfigs(std::vector<SeedLayers_t>  aSeedLayerConfigs)
{
  
  SeedLayerConfigs.clear();
  for(int i=0;i<int(aSeedLayerConfigs.size());i++) {
    bool IsGoodConfig = true;
    if(aSeedLayerConfigs[i].SeedLayers.size() != 3) IsGoodConfig = false;
    for(int j=0;j<int(aSeedLayerConfigs[i].SeedLayers.size());j++) {
      for(int k=0;k<int(aSeedLayerConfigs[i].SeedLayers.size());k++) {
	if(j >= k) continue;
	if(aSeedLayerConfigs[i].SeedLayers[j] == aSeedLayerConfigs[i].SeedLayers[k]) {
	  IsGoodConfig = false;
	  break;
	}
      }
      if(!IsGoodConfig) break;
    }
    
    if(!IsGoodConfig) continue;
    
    SeedLayerConfigs.push_back(aSeedLayerConfigs[i]);
  }
  
  if(SeedLayerConfigs.size() == 0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::FillSeedConfigs:: SeedLayerConfigs has zero entries. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::CheckInputs()
{
  
  if(Nhits_min <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::CheckInputs:: Nhits_min is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(PtMin_cut <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::CheckInputs:: PtMin_cut is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(PurityMin_cut <= 0.0 || PurityMin_cut > 1.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::CheckInputs:: PurityMin_cut is outside allowed range (0,1]. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Chi2Ondf_Seed_cut <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::CheckInputs:: Chi2Ondf_Seed_cut is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(Chi2Ondf_Add_cut <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::CheckInputs:: Chi2Ondf_Add_cut is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  if(NfakesMax_seeding < 0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::CheckInputs:: NfakesMax_seeding is smaller than zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetNhitsMin(int aNhits_min) 
{ 
  
  if(aNhits_min <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetNhitsMin:: Nhits_min is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Nhits_min = aNhits_min; 
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetPtMinCut(double aPtMin_cut) 
{ 
  
  if(aPtMin_cut <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetPtMinCut:: PtMin_cut is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  PtMin_cut = aPtMin_cut; 
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetPurityMinCut(double aPurityMin_cut)
{ 
  
  if(aPurityMin_cut <= 0.0 || aPurityMin_cut > 1.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetPurityMinCut:: PurityMin_cut is outside allowed range (0,1]. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  PurityMin_cut = aPurityMin_cut;
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetChi2OndfSeedCut(double aChi2Ondf_Seed_cut)
{ 
  
  if(aChi2Ondf_Seed_cut <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetChi2OndfSeedCut:: Chi2Ondf_Seed_cut is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Chi2Ondf_Seed_cut = aChi2Ondf_Seed_cut;
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetChi2OndfAddCut(double aChi2Ondf_Add_cut) 
{ 
  
  if(aChi2Ondf_Add_cut <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetChi2OndfAddCut:: Chi2Ondf_Add_cut is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Chi2Ondf_Add_cut = aChi2Ondf_Add_cut; 
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetInwardTracking(bool aInwardTracking)
{ 
  
  InwardTracking = aInwardTracking;
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetNmcSeedEffic(int aNmc_seed_effic)
{ 
  
  if(aNmc_seed_effic <= 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetNmcSeedEffic:: Nmc_seed_effic is smaller than or equal to zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  Nmc_seed_effic = aNmc_seed_effic;
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::SetNfakesMaxSeeding(int aNfakesMax_seeding)
{ 
  
  if(aNfakesMax_seeding < 0.0) {
    cout << endl;
    cout << "ERROR in GTrackFinderAlgoFPCCD::SetNfakesMaxSeeding:: NfakesMax_seeding is smaller than zero. Check your inputs. Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  NfakesMax_seeding = aNfakesMax_seeding;
  
  return;
  
}
//====================================================================
void  GTrackFinderAlgoFPCCD::Print()
{
  
  cout << "    Begin FPCCD TrackFinderAlgo" << endl;
  cout << "      AlgoName          " << Name.Data() << endl;
  cout << "      AlgoType          " << Type.Data() << endl;
  cout << "      Region:           " << endl;
  cout << "        Position:" << endl;
  cout << "          X-Range     = (" << Region.posRange.Rx[0]/global->GetUnit("cm")      << "," << Region.posRange.Rx[1]/global->GetUnit("cm")      << ") cm" << endl;
  cout << "          Y-Range     = (" << Region.posRange.Ry[0]/global->GetUnit("cm")      << "," << Region.posRange.Ry[1]/global->GetUnit("cm")      << ") cm" << endl;
  cout << "          Z-Range     = (" << Region.posRange.Rz[0]/global->GetUnit("cm")      << "," << Region.posRange.Rz[1]/global->GetUnit("cm")      << ") cm" << endl;
  cout << "          R-Range     = (" << Region.posRange.Rr[0]/global->GetUnit("cm")      << "," << Region.posRange.Rr[1]/global->GetUnit("cm")      << ") cm" << endl;
  cout << "          Theta-Range = (" << Region.posRange.Rtheta[0]/global->GetUnit("deg") << "," << Region.posRange.Rtheta[1]/global->GetUnit("deg") << ") deg" << endl;
  cout << "          Phi-Range   = (" << Region.posRange.Rphi[0]/global->GetUnit("deg")   << "," << Region.posRange.Rphi[1]/global->GetUnit("deg")   << ") deg" << endl;
  cout << "        Momentum:" << endl;
  cout << "          P-Range     = (" << Region.momRange.Rp[0]/global->GetUnit("GeV/c")   << "," << Region.momRange.Rp[1]/global->GetUnit("GeV/c")   << ") GeV/c" << endl;
  cout << "          Theta-Range = (" << Region.momRange.Rtheta[0]/global->GetUnit("deg") << "," << Region.momRange.Rtheta[1]/global->GetUnit("deg") << ") deg"   << endl;
  cout << "          Phi-Range   = (" << Region.momRange.Rphi[0]/global->GetUnit("deg")   << "," << Region.momRange.Rphi[1]/global->GetUnit("deg")   << ") deg"   << endl;
  
  cout << "      Systems list:" << endl;
  for(int i=0;i<int(FullSystemList.size());i++) {
    cout << "        " << FullSystemList[i].Data() << endl;
  }
  cout << "      Seed configs list:" << endl;
  for(int i=0;i<int(SeedLayerConfigs.size());i++) {
    cout << "        ";
    for(int j=0;j<int(SeedLayerConfigs[i].SeedLayers.size());j++) {
      cout << SeedLayerConfigs[i].SeedLayers[j].Data();
      if(j+1 < int(SeedLayerConfigs[i].SeedLayers.size())) cout << "    ";
    }
    cout << endl;
  }
  
  cout << "      PtMin             " << PtMin_cut/global->GetUnit("GeV/c") << " GeV/c" << endl;
  cout << "      NhitsMin          " << Nhits_min << endl;
  cout << "      PurityMin         " << PurityMin_cut << endl;
  cout << "      Chi2OndfMinSeed   " << Chi2Ondf_Seed_cut << endl;
  cout << "      Chi2OndfMinAdd    " << Chi2Ondf_Add_cut << endl;
  cout << "      InwardTracking    ";
  if(InwardTracking) cout << "true";
  else               cout << "false";
  cout << endl;
  cout << "      NfakesSeedMax     " << NfakesMax_seeding << endl;
  cout << "      NmcEfficSeeding   " << Nmc_seed_effic << endl;
  cout << "    End   FPCCD TrackFinderAlgo" << endl;
  
  return;
  
}
//====================================================================

