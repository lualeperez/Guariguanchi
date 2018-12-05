////////////////////////////////////////////////////////
// In this file are handled all the analysis funtions //
////////////////////////////////////////////////////////

#include <TGraph.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
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
#include <TLatex.h>
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
void  Guariguanchi::DoAnalysis(void)
{
  
  //Analysis manager funtion. This function does the following,
  // - Calls for the reading of the input datacard
  // - Calls the filling of the geometries specified in the input datacard and calls for some checks
  // - Calls the functions for analysis depending on the specifications on the input datacard

  ReadDataCard();

  FillGeometries();
    
  FillTrackers();
    
  if(DoGeometryCheck) CheckGeometries();
  
  TString command = TString("rm ") + TheOutputFile + TString(".root");
  gSystem->Exec(command.Data());
    
  if(DoPrintGeometry) PrintGeometries();
  
  if(DoPrintGeometryWeight) PrintGeometryWeights();

  if(DoPlotGeometry)  PlotGeometries();
  
  //Material budget calculation vs momentum, theta and phi
  if(DoMaterialBugdetAnalysis) doMaterialBudgetAnalysis();
  
  // Doing track resolution parameters analysis
  if(DoTrkResolAnalysis || DoTelescopeAnalysis)  doTrkResolAnalysis();
  
  // Pseudo-efficiency vs momentum
  if(DoPseudoEfficVsMon) doTrackingPseudoEfficVsMomentum();

  BuildSingleFile();

  return;
  
}
//====================================================================
void  Guariguanchi::doMaterialBudgetAnalysis(void)
{

  //This function calculates the material budget that a particle of given momentum encounters
  //and show the results as a function of p, theta and phi for the different geometries

  cout << endl;
  cout << "Start doMaterialBudgetAnalysis  ";
  global->fWatch.Print();
  global->fWatch.Continue();

  GlobalFileCounter++;
  TString HistName,HistTitle;
  //char ytitle[300];
  
  TLatex* latex = new TLatex();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.08);
  latex->SetTextColor(1);
  
  TString MonUnits("GeV/c");
  TString ThetaUnits("deg");

  const int NGeometries(GeometryList.size());
  const int NthetaPoints(thetaArr.size());
  const int NphiPoints(phiArr.size());
  
  std::vector<float> momValuesList;
  momValuesList.clear();
  momValuesList.push_back(MatBudgetAnalysisPars[0].mom_min);
  momValuesList.push_back(0.5*(MatBudgetAnalysisPars[0].mom_min + MatBudgetAnalysisPars[0].mom_max));
  momValuesList.push_back(MatBudgetAnalysisPars[0].mom_max);
  const int NmomValues(momValuesList.size());
  
  const int MaxNSystems(50);
  for(int igeo=0;igeo<NGeometries;igeo++) {
    if(GeometryList[igeo]->GetNSystemNames() > MaxNSystems) {
      cout << endl;
      cout << "Number of Systems for geometry " << GeometryList[igeo]->GetName() << " is higher than maximum allowed number " << MaxNSystems << endl;
      cout << endl;
      
      return;
    }
  }
  
  TString MyPolarVariable;
  TString MyPolarVariableWithUnits;
  if(polarVariable == TString("theta")) {
    MyPolarVariable          = TString("#theta");
    MyPolarVariableWithUnits = TString("#theta (") + ThetaUnits + TString(")");
  }
  else if(polarVariable == TString("costheta")) {
    MyPolarVariable          = TString("cos(#theta)");
    MyPolarVariableWithUnits = TString("cos(#theta)");
  }
  else if(polarVariable == TString("eta")) {
    MyPolarVariable          = TString("#eta");
    MyPolarVariableWithUnits = TString("#eta");
  }
  
  TGraph*  gr_mat_budget_vs_polarVar[NGeometries][NmomValues][MaxNSystems];
  TGraph*  gr_mat_budget_vs_polarVar_all[NGeometries][NmomValues];
  TH1F*    h_mat_budget_vs_polarVar[NGeometries][NmomValues][MaxNSystems];
  TH1F*    h_mat_budget_vs_polarVar_all[NGeometries][NmomValues];
  
  for(int igeo=0;igeo<NGeometries;igeo++) { //begin loop over geometries
    for(int imom=0;imom<NmomValues;imom++) { //begin loop over momentum values
      char themom[300];
      sprintf(themom,"mom = %.1f GeV/c",momValuesList[imom]/global->GetMomentumUnit("GeV/c"));
      
      gr_mat_budget_vs_polarVar_all[igeo][imom] = new TGraph();
      HistName   = TString("gr_mat_budget_vs_polarVar_all_geo") + long(igeo+1) + TString("_mom") + long(imom+1);
      HistTitle  = TString("Total material budget vs ") + MyPolarVariable  + TString(", for ") + TString(themom) + TString(", and geo ") + GeometryList[igeo]->GetName();
      gr_mat_budget_vs_polarVar_all[igeo][imom]->SetName(HistName.Data());
      gr_mat_budget_vs_polarVar_all[igeo][imom]->SetTitle(HistTitle.Data());
      gr_mat_budget_vs_polarVar_all[igeo][imom]->SetMarkerColor(kBlack);
      gr_mat_budget_vs_polarVar_all[igeo][imom]->SetLineColor(kBlack);
      gr_mat_budget_vs_polarVar_all[igeo][imom]->SetLineWidth(2);
    } //end loop over momentum values
    
    int counter = 0;
    for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) { //begin loop over geometry systems
      int mycolor = (counter+1==10)?49:(counter+1);
      while(mycolor == 1 || mycolor == 5 || mycolor == 7) {
	mycolor++;
	counter++;
      }
      if(mycolor == 3) mycolor = kGreen+2;
      counter++;
      
      for(int imom=0;imom<NmomValues;imom++) { //begin loop over momentum values
        char themom[300];
        sprintf(themom,"mom = %.1f GeV/c",momValuesList[imom]/global->GetMomentumUnit("GeV/c"));
      
        gr_mat_budget_vs_polarVar[igeo][imom][isys] = new TGraph();
        HistName    = TString("gr_mat_budget_vs_polarVar_geo") + long(igeo+1) + TString("_mom") + long(imom+1);
        HistTitle   = TString("Material budget for system ") + GeometryList[igeo]->GetASystemName(isys);
	HistTitle  += TString("vs ") + MyPolarVariable  + TString(", for ") + TString(themom) + TString(", and geo ") + GeometryList[igeo]->GetName();
        gr_mat_budget_vs_polarVar[igeo][imom][isys]->SetName(HistName.Data());
        gr_mat_budget_vs_polarVar[igeo][imom][isys]->SetTitle(HistTitle.Data());
        gr_mat_budget_vs_polarVar[igeo][imom][isys]->SetMarkerColor(mycolor);
        gr_mat_budget_vs_polarVar[igeo][imom][isys]->SetLineColor(mycolor);
	gr_mat_budget_vs_polarVar[igeo][imom][isys]->SetFillColor(mycolor);
        gr_mat_budget_vs_polarVar[igeo][imom][isys]->SetLineWidth(2);
      } //end loop over momentum values
      
    } //end loop over geometry systems
    
  } //end loop over geometries
  
  double RpolarVar[2];
  RpolarVar[0] = +1.0e+20;
  RpolarVar[1] = -1.0e+20;
  double RmatBud[2];
  RmatBud[0] = +1.0e+20;
  RmatBud[1] = -1.0e+20;
  
  
  fPrintFreq = 5;
  for(int itheta=0;itheta<NthetaPoints;itheta++) { //Begin of loop over theta values
    int theta_kkk = itheta;
    if(polarVariable == TString("costheta") || polarVariable == TString("eta")) theta_kkk = NthetaPoints - itheta - 1;
    
    double polarValue = thetaArr[theta_kkk]/global->GetAngleUnit(ThetaUnits);
    if(polarVariable == TString("theta"))          polarValue = thetaArr[theta_kkk]/global->GetAngleUnit(ThetaUnits);
    else if(polarVariable == TString("costheta"))  polarValue = global->FromThetaToCosTheta(thetaArr[theta_kkk]);
    else if(polarVariable == TString("eta"))       polarValue = global->FromThetaToEta(thetaArr[theta_kkk]);
    
    if(RpolarVar[0] > polarValue) RpolarVar[0] = polarValue;
    if(RpolarVar[1] < polarValue) RpolarVar[1] = polarValue;
    
    for(int ip=0;ip<NmomValues;ip++) { //Begin of loop over momentum values
      
      for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
	double  Total_material_budget = 0.0;
	
	const int Nsystems(GeometryList[igeo]->GetNSystemNames());
	double  material_budget_per_systems[Nsystems];
	for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) material_budget_per_systems[isys] = 0.0;
	
	for(int iphi=0;iphi<NphiPoints;iphi++) { //Begin of loop over phi values
	  TVector3 momentum = global->GetMomentum(global->GetMomentumFromMomVar(momValuesList[ip],particle,momVariable),thetaArr[theta_kkk],phiArr[iphi]);
	  TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);
	  
	  double tot_mat_budget_temp = 0;
	  std::vector<double>  system_mat_budget_tmp;
	  TrackerList[igeo]->GetGeometrySystemsMaterialBudget(system_mat_budget_tmp, tot_mat_budget_temp);
	  
	  Total_material_budget += tot_mat_budget_temp*100;
	  for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) {
	    material_budget_per_systems[isys] += system_mat_budget_tmp[isys]*100;
	  }
	  
	} //end loop over phi values
	Total_material_budget /= NphiPoints;
	gr_mat_budget_vs_polarVar_all[igeo][ip]->SetPoint(gr_mat_budget_vs_polarVar_all[igeo][ip]->GetN(),polarValue,Total_material_budget);
	
	if(RmatBud[0] > Total_material_budget)  RmatBud[0] = Total_material_budget;
	if(RmatBud[1] < Total_material_budget)  RmatBud[1] = Total_material_budget;
	
	for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) {
	  material_budget_per_systems[isys] /= NphiPoints;
	  gr_mat_budget_vs_polarVar[igeo][ip][isys]->SetPoint(gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetN(),polarValue,material_budget_per_systems[isys]);
	}
	
      } //end loop over geometries
      
    } //end loop over momentum values
    
    if(!((itheta+1)%fPrintFreq)) {
      char thetaTitle[300];
      if(polarVariable == TString("theta"))         sprintf(thetaTitle,"theta = %.1f deg",thetaArr[theta_kkk]/global->GetAngleUnit("deg"));
      else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(theta) = %.3f",global->FromThetaToCosTheta(thetaArr[theta_kkk]));
      else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"eta = %.3f",global->FromThetaToEta(thetaArr[theta_kkk]));
    
      cout << itheta+1 << " " << thetaTitle << " processed (out of " << NthetaPoints << ")!!!  ";
      global->fWatch.Print();
      global->fWatch.Continue();
    }
  } //end loop over theta values

  cout << "End of loop in doMaterialBudgetAnalysis function";
  global->fWatch.Print();
  global->fWatch.Continue();
  cout << endl;
  
  int Nbins = 500;
  
  bool NoRange = false;  
  double delta,porcent,min,min_frac;
  min_frac = 0.6;
  
  porcent = 0.10;
  delta   = RpolarVar[1] - RpolarVar[0];
  if(TMath::Abs(delta) < 1.0e-6) {
    if(polarVariable == TString("theta"))         delta = 1.0*global->GetUnit("deg");
    else if(polarVariable == TString("costheta")) delta = 0.1;
    else if(polarVariable == TString("eta"))      delta = 0.1;
  }
  else delta *= porcent;
  RpolarVar[0] -= delta;
  RpolarVar[1] += delta;
  
  for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
    for(int ip=0;ip<NmomValues;ip++) { //Begin of loop over momentum values
      char themom[300];
      TString MomUnit_tmp("GeV/c");
      if(momVariable == TString("E") || momVariable == TString("Ekin") || momVariable == TString("EkinPerU")) MomUnit_tmp = TString("GeV");
      sprintf(themom,"%.1f %s",momValuesList[ip]/global->GetUnit(MomUnit_tmp),MomUnit_tmp.Data());
      
      TString Title_tmp("p = ");
      if(momVariable == TString("E"))             Title_tmp = TString("E = ");
      else if(momVariable == TString("Ekin"))     Title_tmp = TString("E_{kin} = ");
      else if(momVariable == TString("EkinPerU")) Title_tmp = TString("E_{kin}/u = ");
      
      HistName   = TString("h_mat_budget_vs_polarVar_all_geo") + long(igeo+1) + TString("_") + momVariable + long(ip+1);
      HistTitle  = TString("geo ") + GeometryList[igeo]->GetName() + TString(" (") + Title_tmp + TString(themom) + TString(")");
      h_mat_budget_vs_polarVar_all[igeo][ip] = new TH1F(HistName.Data(),HistTitle.Data(),
							Nbins,RpolarVar[0],RpolarVar[1]);
      HistName = MyPolarVariableWithUnits;
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetXTitle(HistName.Data());
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetXaxis()->CenterTitle(true);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetYTitle("Material Budget (%)");
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetYaxis()->CenterTitle(true);
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetYaxis()->SetTitleOffset(TitleOffSet);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetLineColor(1);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetLineWidth(2);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetFillColor(kWhite);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetMinimum(RmatBud[0]);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetMaximum(RmatBud[1]);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetNdivisions(5 + 100*5,"X");
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetNdivisions(5 + 100*5,"Y");
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetXaxis()->SetTitleSize(TheSize);
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetXaxis()->SetLabelSize(TheSize);
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetYaxis()->SetTitleSize(TheSize);
      h_mat_budget_vs_polarVar_all[igeo][ip]->GetYaxis()->SetLabelSize(TheSize);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetStats(false);
      
      for(int ix=0;ix<Nbins;ix++) {
	double x = h_mat_budget_vs_polarVar_all[igeo][ip]->GetBinCenter(ix+1);
	int idx = -999;
	for(int kkk=1;kkk<gr_mat_budget_vs_polarVar_all[igeo][ip]->GetN();kkk++) {
	  double x_kkkm1,x_kkk,y;
	  gr_mat_budget_vs_polarVar_all[igeo][ip]->GetPoint(kkk-1,x_kkkm1,y);
	  gr_mat_budget_vs_polarVar_all[igeo][ip]->GetPoint(kkk,  x_kkk,  y);
	  if(x >= x_kkkm1 && x <= x_kkk) {
	    idx = kkk-1;
	    break;
	  }
	}
	if(idx < 0) continue;
	
	double x1,y1,x2,y2,a,b;
	gr_mat_budget_vs_polarVar_all[igeo][ip]->GetPoint(idx,  x1,y1);
	gr_mat_budget_vs_polarVar_all[igeo][ip]->GetPoint(idx+1,x2,y2);
	a = (y2-y1)/(x2-x1);
	b = y2 - a*x2;
	h_mat_budget_vs_polarVar_all[igeo][ip]->SetBinContent(ix+1,a*x+b);
      }
      
      for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) { //begin loop over systems
	HistName   = TString("h_mat_budget_vs_polarVar_geo") + long(igeo+1) + TString("_mom") + long(ip+1) + TString("_sys") + long(isys+1);
	HistTitle   = TString("Material budget for system ") + GeometryList[igeo]->GetASystemName(isys);
	HistTitle  += TString("vs ") + MyPolarVariable  + TString(", for ") + TString(themom) + TString(", and geo ") + GeometryList[igeo]->GetName();
        h_mat_budget_vs_polarVar[igeo][ip][isys] = new TH1F(HistName.Data(),HistTitle.Data(),
							    Nbins,RpolarVar[0],RpolarVar[1]);
        HistName = MyPolarVariableWithUnits;
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetXTitle(HistName.Data());
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetXaxis()->CenterTitle(true);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetYTitle("Material Budget (%)");
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetYaxis()->CenterTitle(true);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetYaxis()->SetTitleOffset(TitleOffSet);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetLineColor(1);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetLineWidth(1);
	h_mat_budget_vs_polarVar[igeo][ip][isys]->SetFillColor(gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetLineColor());
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetMinimum(RmatBud[0]);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetMaximum(RmatBud[1]);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetNdivisions(5 + 100*5,"X");
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetNdivisions(5 + 100*5,"Y");
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetXaxis()->SetTitleSize(TheSize);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetXaxis()->SetLabelSize(TheSize);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetYaxis()->SetTitleSize(TheSize);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->GetYaxis()->SetLabelSize(TheSize);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetStats(false);
      }// end loop over systems
      
      for(int i=0;i<gr_mat_budget_vs_polarVar[igeo][ip][0]->GetN();i++) { //begin loop over polar-variable values
	for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) { //begin loop over systems
	  double x,y_current,y_prev;
	  
	  gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetPoint(i,x,y_current);
	  
	  
	  if(isys > 0) {
	    gr_mat_budget_vs_polarVar[igeo][ip][isys-1]->GetPoint(i,x,y_prev);
	    y_current += y_prev;
	    gr_mat_budget_vs_polarVar[igeo][ip][isys]->SetPoint(i,x,y_current);
	  }
	  
	  if(RmatBud[0] > y_current)  RmatBud[0] = y_current;
	  if(RmatBud[1] < y_current)  RmatBud[1] = y_current;
	} //end loop over systems
      } //end loop over polar-variable values
      
      for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) { //begin loop over systems
        for(int ix=0;ix<Nbins;ix++) {
	  double x = h_mat_budget_vs_polarVar[igeo][ip][isys]->GetBinCenter(ix+1);
	  int idx = -999;
	  for(int kkk=1;kkk<gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetN();kkk++) {
	    double x_kkkm1,x_kkk,y;
	    gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetPoint(kkk-1,x_kkkm1,y);
	    gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetPoint(kkk,  x_kkk,  y);
	    if(x >= x_kkkm1 && x <= x_kkk) {
	      idx = kkk-1;
	      break;
	    }
	  }
	  if(idx < 0) continue;
	
	  double x1,y1,x2,y2,a,b;
	  gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetPoint(idx,  x1,y1);
	  gr_mat_budget_vs_polarVar[igeo][ip][isys]->GetPoint(idx+1,x2,y2);
	  a = (y2-y1)/(x2-x1);
	  b = y2 - a*x2;
	  h_mat_budget_vs_polarVar[igeo][ip][isys]->SetBinContent(ix+1,a*x+b);
	}
      }
      
    } // end loop over momentum values
  } // end loop over geometries
  
  porcent = 0.10;
  delta = RpolarVar[1] - RpolarVar[0];
  if(TMath::Abs(delta) < 1.0e-6) {
    if(polarVariable == TString("theta"))         delta = 1.0*global->GetUnit("deg");
    else if(polarVariable == TString("costheta")) delta = 0.1;
    else if(polarVariable == TString("eta"))      delta = 0.1;
  }
  else delta *= porcent;
  RpolarVar[0] -= delta;
  RpolarVar[1] += delta;
  
  NoRange = false;
  if(RmatBud[0] == +1.0e+20 && RmatBud[1] == -1.0e+20) NoRange = true;
  if(NoRange) {
    RmatBud[0] = 1.0e-6;
    RmatBud[1] = 1.0;
  }
  else {
    porcent = 0.10;
    delta   = RmatBud[1] - RmatBud[0];
    if(TMath::Abs(delta) < 1.0e-8) delta = 1.0e-3;
    else                           delta *= porcent;
    min = RmatBud[0];
    RmatBud[0] -= delta;
    RmatBud[1] += delta;
    if(RmatBud[0] < 0.0) {
      if(min > 1.0e-6) RmatBud[0] = min*min_frac;
      else             RmatBud[0] = 1.0e-6;
    }
  }
  
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int ip=0;ip<NmomValues;ip++) {
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetMinimum(RmatBud[0]);
      h_mat_budget_vs_polarVar_all[igeo][ip]->SetMaximum(RmatBud[1]);
      for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) {
	h_mat_budget_vs_polarVar[igeo][ip][isys]->SetMinimum(RmatBud[0]);
        h_mat_budget_vs_polarVar[igeo][ip][isys]->SetMaximum(RmatBud[1]);
      }
    }
  }
  
  float leftMargin    = 0.20;
  float bottomMargin  = 0.15;
  float topMargin     = 0.10;
  
  TLegend* leg_syst_geo[NGeometries];
  
  TString EPSName  = TheOutputFile + TString("_") + long(GlobalFileCounter) + TString(".eps");
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");

  TCanvas* c1 = new TCanvas("c1","c1",LargeCanvasX,LargeCanvasY);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.10);
  c1->SetBottomMargin(0.10);
  //c1->SetRightMargin(0.10);

  c1->Print(EPSNameO.Data());
  
  for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
    leg_syst_geo[igeo] = new TLegend(0.1,0.1,0.9,0.9);
    leg_syst_geo[igeo]->SetFillColor(10);
    for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) {
      HistName = GeometryList[igeo]->GetASystemName(isys);
      leg_syst_geo[igeo]->AddEntry(h_mat_budget_vs_polarVar[igeo][0][isys],HistName.Data(),"f");
    }
    leg_syst_geo[igeo]->AddEntry(h_mat_budget_vs_polarVar_all[igeo][0],"Total","l");
    
    c1->Clear();
    c1->Divide(2,2);
    for(int ip=0;ip<NmomValues;ip++) {
      c1->cd(ip+1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(leftMargin);
      gPad->SetBottomMargin(bottomMargin);
      gPad->SetTopMargin(topMargin);
      h_mat_budget_vs_polarVar_all[igeo][ip]->Draw("l");
      for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) {
	int idx = GeometryList[igeo]->GetNSystemNames() - isys - 1;
	h_mat_budget_vs_polarVar[igeo][ip][idx]->Draw("lsame");
      }
      h_mat_budget_vs_polarVar_all[igeo][ip]->Draw("AXISsame");
    }
    c1->cd(4);
    leg_syst_geo[igeo]->Draw();
    
    c1->cd();
    PlotLogo(0.15,0.89,0.51,0.5);
    c1->Print(EPSName.Data());
  } //end loop over geometries
  
  c1->Print(EPSNameC.Data());

  if(SavePlots) {
    cout << endl;
    TString ROOTName = TheOutputFile + TString(".root");
    cout << "Saving track parameters resolution vs momentum to " << ROOTName.Data() << " file" << endl;
    TFile file(ROOTName.Data(),"UPDATE");

    TCanvas* c_material_budget[NGeometries];
    
    for(int igeo=0;igeo<NGeometries;igeo++) { //Begin of loop over geometries
      HistName  = TString("c_material_budget_geo") + long(igeo+1);
      HistTitle = TString("MAterial budget vs ") + MyPolarVariable + TString(" for geo ") + GeometryList[igeo]->GetName();
      c_material_budget[igeo] = new TCanvas(HistName.Data(),HistTitle.Data());
      c_material_budget[igeo]->SetFillColor(10);
      c_material_budget[igeo]->SetFrameFillColor(10);
      c_material_budget[igeo]->SetTickx(1);
      c_material_budget[igeo]->SetTicky(1);
      c_material_budget[igeo]->SetLeftMargin(0.10);
      c_material_budget[igeo]->SetBottomMargin(0.10);
    
      c_material_budget[igeo]->Clear();
      c_material_budget[igeo]->Divide(2,2);
      for(int ip=0;ip<NmomValues;ip++) {
        c_material_budget[igeo]->cd(ip+1);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(leftMargin);
        gPad->SetBottomMargin(bottomMargin);
        gPad->SetTopMargin(topMargin);
        h_mat_budget_vs_polarVar_all[igeo][ip]->Draw("l");
        for(int isys=0;isys<GeometryList[igeo]->GetNSystemNames();isys++) {
	  int idx = GeometryList[igeo]->GetNSystemNames() - isys - 1;
	  h_mat_budget_vs_polarVar[igeo][ip][idx]->Draw("lsame");
        }
        h_mat_budget_vs_polarVar_all[igeo][ip]->Draw("AXISsame");
      }
      c_material_budget[igeo]->cd(4);
      leg_syst_geo[igeo]->Draw();
      
      c_material_budget[igeo]->cd();
      PlotLogo(0.15,0.89,0.51,0.5);
      
      c_material_budget[igeo]->Write();
    } //end loop over geometries

    file.Close();
  }

  return;
  
}
//====================================================================
void  Guariguanchi::doTrkResolAnalysis(void) 
{

  //This function performs an analytical calculation of the resolution on the track parameters
  //as a function of particles momentum in the specified range as well as for the set polar 
  //angles specified in the input datacard
  //The set of track parameters depend on the configuration, i.e. if there is a magnetic field (helix) or not (stright line)

  CheckGeometriesComparability();
  CheckTelescopeConfigurations();
  
  cout << endl;
  cout << "Start doTrkResolAnalysis  ";
  global->fWatch.Print();
  global->fWatch.Continue();

  int Somebins = Nbins_mom_redu;
  int Thebins  = momArr.size();
  if(Thebins > Somebins) Thebins = Somebins;
  if(TrkResolAnalysisPars[0].UseAllMomVals) Thebins  = momArr.size();
  
  GlobalFileCounter++;
  TString HistName,HistTitle;
  char ytitle[300];
  
  TLatex* latex = new TLatex();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.08);
  latex->SetTextColor(1);
  
  const int Nmax_params(8);
  std::vector<TString> title;
  std::vector<TString> titleY;
  std::vector<TString> ParamUnits;
  title.clear();
  titleY.clear();
  ParamUnits.clear();
  
  TString MonUnits("GeV/c");
  if(momVariable == TString("E"))             MonUnits = TString("GeV");
  else if(momVariable == TString("Ekin"))     MonUnits = TString("GeV");
  else if(momVariable == TString("EkinPerU")) MonUnits = TString("GeV");
  TString ThetaUnits("deg");

  const int NGeometries(GeometryList.size());
  const int NthetaPoints(thetaArr.size());
  const int NphiPoints(phiArr.size());
  const int NmomPoints(momArr.size());
  const int MyNbins_mom_redu(Thebins);
  
  if(NthetaPoints == 1) TrkResolAnalysisPars[0].PlotPerformancesVsTheta = false;
  if(NphiPoints   == 1) {
    TrkResolAnalysisPars[0].PlotPhiAveraged     = false;
    TrkResolAnalysisPars[0].PlotOnlyPhiAveraged = false;    
  }
  
  TString  MyPolarVariable;
  if(polarVariable == TString("theta"))         MyPolarVariable = TString("#theta");
  else if(polarVariable == TString("costheta")) MyPolarVariable = TString("cos(#theta)");
  else if(polarVariable == TString("eta"))      MyPolarVariable = TString("#eta");
  
  //////////////////////////////
  ///// Defining the plots /////
  //////////////////////////////
  
  //Plotting a and b parameters
  TGraph* gr_a[NGeometries];
  TGraph* gr_a_vs_Phi[NGeometries][NthetaPoints];
  TGraph* gr_Nhits_vs_Phi[NGeometries][NthetaPoints];
  TGraph* gr_b[NGeometries];
  
  TGraph*       gr[NGeometries][NthetaPoints][NphiPoints][Nmax_params];
  TGraph*       gr_Nhits_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  TGraph*       gr_MatBudget_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  TGraphErrors* gr_aveVsPhi_ErrorRMS[NGeometries][NthetaPoints][Nmax_params];
  TGraphErrors* gr_aveVsPhi[NGeometries][NthetaPoints][Nmax_params];
  TGraph*       gr_aveVsPhi_NoErr[NGeometries][NthetaPoints][Nmax_params];
  TGraphErrors* gr_aveSigmaParVsTheta[NGeometries][MyNbins_mom_redu][Nmax_params];
  TGraphErrors* gr_aveNhitsVsTheta[NGeometries][MyNbins_mom_redu];
  TGraph*       gr_1stPorintDisToRef_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  TGraph*       gr_aveVsPhi1stPorintDisToRef_vs_mom[NGeometries][NthetaPoints];
  TGraph*       gr_aveVsPhiAndMom1stPorintDisToRef[NGeometries];
  
  //Graphs for telescope analyses
  const int     MaxDUTs(10);

  int**         NDUTS;
  TString***    DUTNames;
  TString*****  DUTLayerNames;
  float*****    TelResolUAtDUT;
  float*****    TelResolVAtDUT;
  float*****    TelResolUVAreaAtDUT;
  float*****    NbkgAtDUT;
  float*****    MomValue;
  
  TGraph*****   gr_TelResolUAtDUT;
  TGraph*****   gr_TelResolVAtDUT;
  TGraph*****   gr_TelResolUVAreaAtDUT;
  TGraph*****   gr_NbkgAtDUT;
  
  if(DoTelescopeAnalysis) {
    NDUTS    = new int*[NthetaPoints];
    DUTNames = new TString**[NthetaPoints];
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      NDUTS[itheta] = new int[NphiPoints];
      DUTNames[itheta] = new TString*[NphiPoints];
      for(int iphi=0;iphi<NphiPoints;iphi++) {
        DUTNames[itheta][iphi] = new TString[MaxDUTs];
      }
    }
  
    DUTLayerNames          = new TString****[NGeometries];
    TelResolUAtDUT         = new float****[NGeometries];
    TelResolVAtDUT         = new float****[NGeometries];
    TelResolUVAreaAtDUT    = new float****[NGeometries];
    NbkgAtDUT              = new float****[NGeometries];
    MomValue               = new float****[NGeometries];
    gr_TelResolUAtDUT      = new TGraph****[NGeometries];
    gr_TelResolVAtDUT      = new TGraph****[NGeometries];
    gr_TelResolUVAreaAtDUT = new TGraph****[NGeometries];
    gr_NbkgAtDUT           = new TGraph****[NGeometries];
    for(int igeo=0;igeo<NGeometries;igeo++) {
      DUTLayerNames[igeo]          = new TString***[NthetaPoints];
      TelResolUAtDUT[igeo]         = new float***[NthetaPoints];
      TelResolVAtDUT[igeo]         = new float***[NthetaPoints];
      TelResolUVAreaAtDUT[igeo]    = new float***[NthetaPoints];
      NbkgAtDUT[igeo]              = new float***[NthetaPoints];
      MomValue[igeo]               = new float***[NthetaPoints];
      gr_TelResolUAtDUT[igeo]      = new TGraph***[NthetaPoints];
      gr_TelResolVAtDUT[igeo]      = new TGraph***[NthetaPoints];
      gr_TelResolUVAreaAtDUT[igeo] = new TGraph***[NthetaPoints];
      gr_NbkgAtDUT[igeo]           = new TGraph***[NthetaPoints];
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
        DUTLayerNames[igeo][itheta]          = new TString**[NphiPoints];
        TelResolUAtDUT[igeo][itheta]         = new float**[NphiPoints];
        TelResolVAtDUT[igeo][itheta]         = new float**[NphiPoints];
        TelResolUVAreaAtDUT[igeo][itheta]    = new float**[NphiPoints];
        NbkgAtDUT[igeo][itheta]              = new float**[NphiPoints];
        MomValue[igeo][itheta]               = new float**[NphiPoints];
	gr_TelResolUAtDUT[igeo][itheta]      = new TGraph**[NphiPoints];
        gr_TelResolVAtDUT[igeo][itheta]      = new TGraph**[NphiPoints];
        gr_TelResolUVAreaAtDUT[igeo][itheta] = new TGraph**[NphiPoints];
        gr_NbkgAtDUT[igeo][itheta]           = new TGraph**[NphiPoints];
        for(int iphi=0;iphi<NphiPoints;iphi++) {
	  DUTLayerNames[igeo][itheta][iphi]          = new TString*[NmomPoints];
          TelResolUAtDUT[igeo][itheta][iphi]         = new float*[NmomPoints];
          TelResolVAtDUT[igeo][itheta][iphi]         = new float*[NmomPoints];
          TelResolUVAreaAtDUT[igeo][itheta][iphi]    = new float*[NmomPoints];
	  NbkgAtDUT[igeo][itheta][iphi]              = new float*[NmomPoints];
          MomValue[igeo][itheta][iphi]               = new float*[NmomPoints];
	  gr_TelResolUAtDUT[igeo][itheta][iphi]      = new TGraph*[MaxDUTs];
          gr_TelResolVAtDUT[igeo][itheta][iphi]      = new TGraph*[MaxDUTs];
          gr_TelResolUVAreaAtDUT[igeo][itheta][iphi] = new TGraph*[MaxDUTs];
          gr_NbkgAtDUT[igeo][itheta][iphi]           = new TGraph*[MaxDUTs];
	  for(int imom=0;imom<NmomPoints;imom++) {
	    DUTLayerNames[igeo][itheta][iphi][imom]       = new TString[MaxDUTs];
            TelResolUAtDUT[igeo][itheta][iphi][imom]      = new float[MaxDUTs];
            TelResolVAtDUT[igeo][itheta][iphi][imom]      = new float[MaxDUTs];
            TelResolUVAreaAtDUT[igeo][itheta][iphi][imom] = new float[MaxDUTs];
	    NbkgAtDUT[igeo][itheta][iphi][imom]           = new float[MaxDUTs];
            MomValue[igeo][itheta][iphi][imom]            = new float[MaxDUTs];
	  }
        }
      }
    }
  }

  TString MyMomentumVariable("p");
  if(momVariable == TString("E"))             MyMomentumVariable = TString("Energy");
  else if(momVariable == TString("Ekin"))     MyMomentumVariable = TString("E_{kin}");
  else if(momVariable == TString("EkinPerU")) MyMomentumVariable = TString("E_{kin}/u");
  
  int counter = 0;
  for(int igeo=0;igeo<NGeometries;igeo++) { //begin loop over geometries
    TrackerList[igeo]->GetTrajectory()->FillParametersNames(MonResolRepresentation);
    if(igeo == 0) {
      for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	HistName = TrackerList[igeo]->GetTrajectory()->GetParameterErrorName(ipar);
	HistName += TString(" vs ") + MyMomentumVariable;
	title.push_back(HistName);
	
	HistName =  TrackerList[igeo]->GetTrajectory()->GetParameterErrorName(ipar) + TString(" (");
	HistName += TrackerList[igeo]->GetTrajectory()->GetParameterErrorUnitTitle(ipar) + TString(")");
	titleY.push_back(HistName);
	
	ParamUnits.push_back(TrackerList[igeo]->GetTrajectory()->GetParameterErrorUnit(ipar));
      }
    }
    
    int mycolor = (counter+1==10)?49:(counter+1);
    if(mycolor == 5 || mycolor == 7) {
      mycolor++;
      counter++;
    }
    if(mycolor == 3) mycolor = kGreen+2;

    counter++;
    
    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //being loop over parameters
      for(int itheta=0;itheta<NthetaPoints;itheta++) { // Begin loop over theta values
	char thetaTitle[100];
	if(polarVariable == TString("theta"))         sprintf(thetaTitle,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
	else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(#theta) = %.1f",global->FromThetaToCosTheta(thetaArr[itheta]));
	else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"#eta = %.1f",global->FromThetaToEta(thetaArr[itheta]));
	  
	for(int iphi=0;iphi<NphiPoints;iphi++) { // begin loop over phi values
	  char phiTitle[100];
	  sprintf(phiTitle,"#phi = %.1f deg",phiArr[iphi]/global->GetAngleUnit("deg"));
	  
          gr[igeo][itheta][iphi][ipar] = new TGraph();
          HistName   = TString("gr_res_par") + long(ipar) + TString("_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_vs_mom_geo_") + long(igeo+1);
          HistTitle  = TString("sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	  gr[igeo][itheta][iphi][ipar]->SetName(HistName.Data());
          gr[igeo][itheta][iphi][ipar]->SetTitle(HistTitle.Data());
          gr[igeo][itheta][iphi][ipar]->SetMarkerColor(mycolor);
          gr[igeo][itheta][iphi][ipar]->SetLineColor(mycolor);
          gr[igeo][itheta][iphi][ipar]->SetLineWidth(2);
	  
	  if(ipar == 0) {
	    gr_Nhits_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_Nhits_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
	    sprintf(ytitle,"%.2f",thetaArr[itheta]/global->GetAngleUnit("deg"));
            HistTitle = TString("N_{hits} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_Nhits_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_Nhits_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_Nhits_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_Nhits_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_Nhits_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	    
	    gr_MatBudget_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_MatBudget_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
            HistTitle = TString("Mat-Budget vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	    
	    gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_1stPorintDisToRef_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
            HistTitle = TString("1^{st} point distance to reference vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	    
	    //TelescopeAnalysis
	    if(DoTelescopeAnalysis) { //begin if doing Telescope analysis
	      for(int idut=0;idut<MaxDUTs;idut++) { //begin loop over maximum number of DUTs
		for(int imom=0;imom<NmomPoints;imom++) {
	          DUTLayerNames[igeo][itheta][iphi][imom][idut]      = TString("");
		  TelResolUAtDUT[igeo][itheta][iphi][imom][idut]      = Dummy_value;
		  TelResolVAtDUT[igeo][itheta][iphi][imom][idut]      = Dummy_value;
		  TelResolUVAreaAtDUT[igeo][itheta][iphi][imom][idut] = Dummy_value;
		  NbkgAtDUT[igeo][itheta][iphi][imom][idut]           = Dummy_value;
		  MomValue[igeo][itheta][iphi][imom][idut]            = Dummy_value;
	        }
		
	        gr_TelResolUAtDUT[igeo][itheta][iphi][idut] = new TGraph();
                HistName = TString("gr_TelResolUAtDUT_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1) + TString("_dut") + long(idut+1);
                HistTitle = TString("Telescope pointing resolution U vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ")
		+ GeometryList[igeo]->GetName();
	        gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetName(HistName.Data());
                gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetTitle(HistTitle.Data());
                gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(mycolor);
                gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetLineColor(mycolor);
                gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetLineWidth(2);
		
		gr_TelResolVAtDUT[igeo][itheta][iphi][idut] = new TGraph();
                HistName = TString("gr_TelResolVAtDUT_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1) + TString("_dut") + long(idut+1);
                HistTitle = TString("Telescope pointing resolution V vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ")
		+ GeometryList[igeo]->GetName();
	        gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetName(HistName.Data());
                gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetTitle(HistTitle.Data());
                gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(mycolor);
                gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetLineColor(mycolor);
                gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetLineWidth(2);
		
		gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut] = new TGraph();
                HistName = TString("gr_TelResolUVAreaAtDUT_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1) + TString("_dut") + long(idut+1);
                HistTitle = TString("area of Telescope pointing resolution vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ")
		+ GeometryList[igeo]->GetName();
	        gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetName(HistName.Data());
                gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetTitle(HistTitle.Data());
                gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(mycolor);
                gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetLineColor(mycolor);
                gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetLineWidth(2);
		
		gr_NbkgAtDUT[igeo][itheta][iphi][idut] = new TGraph();
                HistName = TString("gr_NbkgAtDUT_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1) + TString("_dut") + long(idut+1);
                HistTitle = TString("area of Telescope pointing resolution vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") 
		+ GeometryList[igeo]->GetName();
	        gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetName(HistName.Data());
                gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetTitle(HistTitle.Data());
                gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(mycolor);
                gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetLineColor(mycolor);
                gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetLineWidth(2);
	      } //end of loop over maximum number of DUTs
	      
	    } //end of if doing telescope analysis
	    
	  }
	  
	  if(iphi == 0) {
	    gr_aveVsPhi[igeo][itheta][ipar] = new TGraphErrors();
            HistName = TString("gr_aveVsPhi_res_par") + long(ipar) + TString("_theta") + long(itheta+1) + TString("_vs_mom_geo_") + long(igeo+1);
            HistTitle = TString("#phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
            gr_aveVsPhi[igeo][itheta][ipar]->SetName(HistName.Data());
            gr_aveVsPhi[igeo][itheta][ipar]->SetTitle(HistTitle.Data());
            gr_aveVsPhi[igeo][itheta][ipar]->SetMarkerColor(mycolor);
            gr_aveVsPhi[igeo][itheta][ipar]->SetLineColor(mycolor);
            gr_aveVsPhi[igeo][itheta][ipar]->SetLineWidth(2);
	    
	    gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar] = new TGraphErrors();
            HistName = TString("gr_aveVsPhi_ErrorRMS_res_par") + long(ipar) + TString("_theta") + long(itheta+1) + TString("_vs_mom_geo_") + long(igeo+1);
	    HistTitle = TString("#phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
            gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetName(HistName.Data());
            gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetTitle(HistTitle.Data());
            gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetMarkerColor(mycolor);
            gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetLineColor(mycolor);
            gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetLineWidth(2);
	    
	    gr_aveVsPhi_NoErr[igeo][itheta][ipar] = new TGraph();
            HistName = TString("gr_aveVsPhi_NoErr_res_par") + long(ipar) + TString("_theta") + long(itheta+1) + TString("_vs_mom_geo_") + long(igeo+1);
	    HistTitle = TString("#phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
            gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetName(HistName.Data());
            gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetTitle(HistTitle.Data());
            gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetMarkerColor(mycolor);
            gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetLineColor(mycolor);
            gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetLineWidth(2);
	    
	    if(ipar == 0) {
	      gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta] = new TGraph();
              HistName = TString("gr_aveVsPhi1stPorintDisToRef_vs_mom_theta") + long(itheta+1) + TString("_geo_") + long(igeo+1);
	      HistTitle = TString("#phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
              gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetName(HistName.Data());
              gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetTitle(HistTitle.Data());
              gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetMarkerColor(mycolor);
              gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetLineColor(mycolor);
              gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetLineWidth(2);
	    }
	  }
	  
	} // end loop over phi values
      } // end loop over theta values
      
      for(int ip=0;ip<MyNbins_mom_redu;ip++) { // begin loop over momentum values
	double pi;
	if(TrkResolAnalysisPars[0].UseAllMomVals) pi = momArr[ip];
	else                                      pi = momArr[0] + (ip+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
        pi /= global->GetUnit(MonUnits);
	gr_aveSigmaParVsTheta[igeo][ip][ipar] = new TGraphErrors();
	HistName = TString("gr_aveSigmaParVsTheta_res_par") + long(ipar) + TString("_mom") + long(ip+1) + TString("_geo_") + long(igeo+1);
        sprintf(ytitle,"%.3f",pi);
        HistTitle = TString("#phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyPolarVariable +  
        TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ") + GeometryList[igeo]->GetName();
	gr_aveSigmaParVsTheta[igeo][ip][ipar]->SetName(HistName.Data());
	gr_aveSigmaParVsTheta[igeo][ip][ipar]->SetTitle(HistTitle.Data());
	gr_aveSigmaParVsTheta[igeo][ip][ipar]->SetMarkerColor(mycolor);
	gr_aveSigmaParVsTheta[igeo][ip][ipar]->SetLineColor(mycolor);
	gr_aveSigmaParVsTheta[igeo][ip][ipar]->SetLineWidth(2);
      } // end loop over momentum values
      
    } //end loop over parameters
    
    for(int ip=0;ip<MyNbins_mom_redu;ip++) { //begin loop over momentum values
      double pi;
      if(TrkResolAnalysisPars[0].UseAllMomVals) pi = momArr[ip];
      else                                      pi = momArr[0] + (ip+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
      pi /= global->GetUnit(MonUnits);
      gr_aveNhitsVsTheta[igeo][ip] = new TGraphErrors();
      HistName = TString("gr_aveNhitsVsTheta_mom") + long(ip+1) + TString("_vs_mom_geo_") + long(igeo+1);
      sprintf(ytitle,"%.3f",pi);
      HistTitle = TString("#phi average N_{hits} vs ") + MyPolarVariable + TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ")
      + GeometryList[igeo]->GetName();
      gr_aveNhitsVsTheta[igeo][ip]->SetName(HistName.Data());
      gr_aveNhitsVsTheta[igeo][ip]->SetTitle(HistTitle.Data());
      gr_aveNhitsVsTheta[igeo][ip]->SetMarkerColor(mycolor);
      gr_aveNhitsVsTheta[igeo][ip]->SetLineColor(mycolor);
      gr_aveNhitsVsTheta[igeo][ip]->SetLineWidth(2);
    } //end loop over momentum values
    
    gr_a[igeo] = new TGraph();
    HistName  = TString("gr_a_vs_theta_for_geo_") + long(igeo+1);
    HistTitle = TString("a vs ") + MyPolarVariable + TString(" for geo ") + GeometryList[igeo]->GetName();
    gr_a[igeo]->SetName(HistName.Data());
    gr_a[igeo]->SetTitle(HistTitle.Data());
    gr_a[igeo]->SetMarkerColor(mycolor);
    gr_a[igeo]->SetMarkerStyle(20);
    gr_a[igeo]->SetLineColor(mycolor);
    gr_a[igeo]->SetLineWidth(2);

    gr_b[igeo] = new TGraph();
    HistName = TString("gr_b_vs_theta_for_and_geo_") + long(igeo+1);
    HistTitle = TString("b vs ") + MyPolarVariable + TString(" for geo ") + GeometryList[igeo]->GetName();
    gr_b[igeo]->SetName(HistName.Data());
    gr_b[igeo]->SetTitle(HistTitle.Data());
    gr_b[igeo]->SetMarkerColor(mycolor);
    gr_b[igeo]->SetMarkerStyle(20);
    gr_b[igeo]->SetLineColor(mycolor);
    gr_b[igeo]->SetLineWidth(2);
    
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo] = new TGraph();
    HistName = TString("gr_aveVsPhiAndMom1stPorintDisToRef_geo_") + long(igeo+1);
    HistTitle = TString("Ave 1^{st} point distance to reference vs ") + MyPolarVariable + TString(" for geo ") + GeometryList[igeo]->GetName();
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetName(HistName.Data());
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetTitle(HistTitle.Data());
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetMarkerColor(mycolor);
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetMarkerStyle(20);
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetLineColor(mycolor);
    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetLineWidth(2);
    
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      char thetaTitle[100];
      if(polarVariable == TString("theta"))         sprintf(thetaTitle,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
      else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(#theta) = %.1f",global->FromThetaToCosTheta(thetaArr[itheta]));
      else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"#eta = %.1f",global->FromThetaToEta(thetaArr[itheta]));
      
      gr_a_vs_Phi[igeo][itheta] = new TGraph();
      HistName  = TString("gr_a_vs_phi_for_itheta_") + long(itheta) + TString("_and_geo_") + long(igeo+1);
      HistTitle = TString("a vs #phi for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
      gr_a_vs_Phi[igeo][itheta]->SetName(HistName.Data());
      gr_a_vs_Phi[igeo][itheta]->SetTitle(HistTitle.Data());
      gr_a_vs_Phi[igeo][itheta]->SetMarkerColor(mycolor);
      //gr_a_vs_Phi[igeo][itheta]->SetMarkerStyle(20);
      gr_a_vs_Phi[igeo][itheta]->SetLineColor(mycolor);
      gr_a_vs_Phi[igeo][itheta]->SetLineWidth(2);
      
      gr_Nhits_vs_Phi[igeo][itheta] = new TGraph();
      HistName = TString("gr_Nhits_vs_phi_for_itheta_") + long(itheta) + TString("_and_geo_") + long(igeo+1);
      HistTitle = TString("N_{hits} vs #phi for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
      gr_Nhits_vs_Phi[igeo][itheta]->SetName(HistName.Data());
      gr_Nhits_vs_Phi[igeo][itheta]->SetTitle(HistTitle.Data());
      gr_Nhits_vs_Phi[igeo][itheta]->SetMarkerColor(mycolor);
      //gr_Nhits_vs_Phi[igeo][itheta]->SetMarkerStyle(20);
      gr_Nhits_vs_Phi[igeo][itheta]->SetLineColor(mycolor);
      gr_Nhits_vs_Phi[igeo][itheta]->SetLineWidth(2);
    }
    
  } //end loop over geometries
  
  
  //////////////////////////////////////////////////////////////////
  /// Looping over geometries and momentum, theta and phi values ///
  //////////////////////////////////////////////////////////////////
  
  fPrintFreq = 10;
  //fPrintFreq = 1;
  for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
    char thetaTitle[100];
    if(polarVariable == TString("theta"))         sprintf(thetaTitle,"theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
    else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    char thetaVal[100];
    if(polarVariable == TString("theta"))         sprintf(thetaVal,"%.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
    else if(polarVariable == TString("costheta")) sprintf(thetaVal,"%.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(thetaVal,"%.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
      for(int iP=0; iP<int(momArr.size()); iP++) { // loop on momentum

        double pi = momArr[iP];
	TVector3 momentum = global->GetMomentum(global->GetMomentumFromMomVar(pi,particle,momVariable),thetaArr[itheta],phiArr[iphi]);
        if(verbose) momentum.Print();
        //double MyMon = momentum.Mag();

        for(int igeo=0;igeo<NGeometries;igeo++) { //loop over geometries
	  TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);

	  int     Nhits              = 0;
	  double  Material_budget    = 0.0;
	  double  dist1stPointToRef  = 0.0;
	  TMatrixD  FitCovMatrix;
	  std::vector<TelResolAtDUT_t> TelResolAtDUTList;
	  TelResolAtDUTList.clear();
	  
	  bool IsGoodFit = TrackerList[igeo]->doTrkResolAnalysis(TrkResolAnalysisPars[0].NhitsMin,
								 DoTelescopeAnalysis,
								 Nhits,Material_budget,dist1stPointToRef,
							         FitCovMatrix,
							         TelResolAtDUTList);

	  if(std::isnan(Material_budget)) {
	    cout << "igeo = " << igeo+1 << ", " << thetaTitle << ", phi = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg, " 
	         << MyMomentumVariable.Data() << " = " << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << ", "
	         << "Material budget = " << Material_budget*100 << " %, "
	         << "Nhits = " << Nhits << ", "
	         << endl;
	  }
	  
	  gr_Nhits_vs_mom[igeo][itheta][iphi]->SetPoint(iP,pi/global->GetUnit(MonUnits),Nhits);
	  gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetPoint(iP,pi/global->GetUnit(MonUnits),Material_budget*100);
	  
	  //if(!IsGoodFit) continue;

	  if(IsGoodFit) gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetPoint(gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetN(),pi/global->GetUnit(MonUnits),dist1stPointToRef/global->GetDistanceUnit("cm"));

          bool IsNAN = false;
	  for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //loop over track parameters
	    if(IsGoodFit) {
	      if(std::isnan(FitCovMatrix(ipar,ipar))) {
		if(verbose) {
	          cout << "WARNNIN:: parameter " << ipar << " (" << titleY[ipar] << ") for " << MyMomentumVariable.Data() << " = " << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << ", "
	               << thetaTitle << ", and phi = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg, Is nan!!!!"
		       << endl;
		}
	        IsNAN = true;
	        break;
	      }
	    }
	  }

          //if(IsGoodFit && !IsNAN) { //begin if good fit
	    if(DoTelescopeAnalysis) { //begin if doing telescope analysis
	      if(int(TelResolAtDUTList.size()) > MaxDUTs) {
		cout << endl;
		cout << "Number of DUT systems for geometry " << GeometryList[igeo]->GetName().Data() 
		     << " and (p," << polarVariable.Data() << ",phi) = (" 
		     << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << ","
		     << thetaVal << ","
		     << phiArr[iphi]/global->GetUnit("deg") << " deg) "
		     << " is higher than maximum allowed number " << MaxDUTs << ". Exiting now!!!"
		     << endl;
	        cout << endl;
	        assert(false);
	      }
	    
	      //Check if two different DUTs have the same system name
	      for(int idut1=0;idut1<int(TelResolAtDUTList.size());idut1++) {
	        if(TelResolAtDUTList[idut1].TelResolU < 0) continue;
	        for(int idut2=0;idut2<int(TelResolAtDUTList.size());idut2++) {
		  if(TelResolAtDUTList[idut2].TelResolU < 0) continue;
		  if(idut1 >= idut2) continue;
		  if(GeometryList[igeo]->GetGeometryElement(TelResolAtDUTList[idut1].geoElement_idx)->GetLayerName() == GeometryList[igeo]->GetGeometryElement(TelResolAtDUTList[idut2].geoElement_idx)->GetLayerName()) {
		    cout << endl;
		    cout << "For geometry " << GeometryList[igeo]->GetName().Data() 
		         << " and (p," << polarVariable.Data() << ",phi) = (" 
		         << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << ","
			 << thetaVal << ","
		         << phiArr[iphi]/global->GetUnit("deg") << " deg) "
		         << "there are two DUT geometry elements with the same System name "
		         << GeometryList[igeo]->GetGeometryElement(TelResolAtDUTList[idut1].geoElement_idx)->GetSystemName().Data() << ". Check your inputs. Exiting now!!!"
		         << endl;
		    cout << endl;
		    assert(false);
		  }
	        }
	      }
	  
	      for(int idut=0;idut<int(TelResolAtDUTList.size());idut++) { //begin loop over DUTs
	        if(TelResolAtDUTList[idut].TelResolU < 0) continue;
		
		int p_idx = gr[igeo][itheta][iphi][0]->GetN();
		if(IsGoodFit && !IsNAN) { //begin if good fit
                  TelResolUAtDUT[igeo][itheta][iphi][p_idx][idut]      = TelResolAtDUTList[idut].TelResolU/global->GetUnit("um");
		  TelResolVAtDUT[igeo][itheta][iphi][p_idx][idut]      = TelResolAtDUTList[idut].TelResolV/global->GetUnit("um");
		  TelResolUVAreaAtDUT[igeo][itheta][iphi][p_idx][idut] = TelResolAtDUTList[idut].Shadow_area/global->GetUnit("um2");
		  NbkgAtDUT[igeo][itheta][iphi][p_idx][idut]           = TelResolAtDUTList[idut].Nbkg;
		}
		else {
                  TelResolUAtDUT[igeo][itheta][iphi][p_idx][idut]      = -1;
		  TelResolVAtDUT[igeo][itheta][iphi][p_idx][idut]      = -1;
		  TelResolUVAreaAtDUT[igeo][itheta][iphi][p_idx][idut] = -1;
		  NbkgAtDUT[igeo][itheta][iphi][p_idx][idut]           = -1;
		}
		DUTLayerNames[igeo][itheta][iphi][p_idx][idut]       = GeometryList[igeo]->GetGeometryElement(TelResolAtDUTList[idut].geoElement_idx)->GetLayerName();
		MomValue[igeo][itheta][iphi][p_idx][idut]            = pi/global->GetUnit(MonUnits);
		
	      } //end of loop over DUTs
	    } //end if doing Telescope analysis
	    
	    
	    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //loop over track parameters
	      //filling up the TGraph with the parameter resolution vs momentum
	      
	      if(IsGoodFit && !IsNAN) { //begin if good fit
	        double y = TrackerList[igeo]->GetTrajectory()->GetParamError(ipar,FitCovMatrix);
	        gr[igeo][itheta][iphi][ipar]->SetPoint(gr[igeo][itheta][iphi][ipar]->GetN(),pi/global->GetUnit(MonUnits),y);
	      }
	      else {
		gr[igeo][itheta][iphi][ipar]->SetPoint(gr[igeo][itheta][iphi][ipar]->GetN(),pi/global->GetUnit(MonUnits),-1);
	      }
            }  //end of loop over track parameters
	  //} // end if good fit

          if(verbose) {
	    cout << "igeo = " << igeo+1 << ", parameters resoluton : ";
	    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //loop over track parameters
	      cout << titleY[ipar].Data() << " = " << TrackerList[igeo]->GetTrajectory()->GetParamError(ipar,FitCovMatrix);
	      if(ipar < TrackerList[igeo]->GetTrajectory()->GetNParameters() - 1) cout << ", ";
	      else                                                                cout << "  ";
              cout << "for " << MyMomentumVariable.Data() << " (bin = " << iP+1 << " of " << momArr.size() << ") = " << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << " ";
	      cout << ", " << thetaTitle << "   ";
	      cout << "and phi value = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg   ";
	    }
            global->fWatch.Print();
            global->fWatch.Continue();
          }

          if(verbose || true) {
            if(!((iP+1)%fPrintFreq)) {
              cout << iP+1 << " " << MyMomentumVariable.Data() << " bins processed (out of " << NmomPoints << ") for polar value " << itheta+1 << " (out of " << thetaArr.size() << ")  " << thetaTitle
                   << " , phi value " << iphi+1 << " (out of " << phiArr.size() << ") phi = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg"
                   << " and geometry " << GeometryList[igeo]->GetName().Data() << " !!!  ";
              global->fWatch.Print();
              global->fWatch.Continue();
            }
          }

          TelResolAtDUTList.clear();
        } //end of loop over geometries
      } //end of loop on momentum
    } //end of loop on phi values
  } //end of loop on theta values

  cout << "End of loop in doTrkResolAnalysis function  ";
  global->fWatch.Print();
  global->fWatch.Continue();
  cout << endl;

  // Check if there was at least one scanned momentum/theta/phi/geometry that 
  // produced a good track parameters covariance matrix
  bool NoGoodCovMatrix = true;
  bool GoodRangeParam[Nmax_params];
  for(int ipar=0;ipar<Nmax_params;ipar++) GoodRangeParam[ipar] = false;

  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
      for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	  for(int i=0;i<gr[igeo][itheta][iphi][ipar]->GetN();i++) {
	    double x,y;
            gr[igeo][itheta][iphi][ipar]->GetPoint(i, x, y);
	    if(y >= 0.0) {
	      NoGoodCovMatrix      = false;
	      GoodRangeParam[ipar] = true;
	      break;
	    }
	  }
	  
	  //if(gr[igeo][itheta][iphi][ipar]->GetN() > 0) {
	  //  NoGoodCovMatrix      = false;
	  //  GoodRangeParam[ipar] = true;
	  //}
	}
      }
    }
  }
  if(NoGoodCovMatrix) {
    cout << endl;
    cout << "Guariguanchi::doTrackingResolutionVsMomentum::WANNING  None of the geometries for the scanned values of momentum/theta/phi produced a good Track param cov matrix." << endl;
    cout << endl;
    return;
  }

  double RX[2];
  RX[0] = +1.0e+20;
  RX[1] = -1.0e+20;
  double RY[Nmax_params][2];
  for(int ipar=0;ipar<Nmax_params;ipar++) {
    RY[ipar][0] = +1.0e+20;
    RY[ipar][1] = -1.0e+20;
  }

  double RNhits_vs_theta_phi_mom[2];
  double RMatBudget_vs_theta_phi_mom[2];
  RNhits_vs_theta_phi_mom[0] = RMatBudget_vs_theta_phi_mom[0] = +1.0e+20;
  RNhits_vs_theta_phi_mom[1] = RMatBudget_vs_theta_phi_mom[1] = -1.0e+20;
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
      for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	  for(int iP=0; iP<gr[igeo][itheta][iphi][ipar]->GetN(); iP++) {
            double x,y;
            gr[igeo][itheta][iphi][ipar]->GetPoint(iP, x, y);
	  
            if(RX[0] > x) RX[0] = x;
            if(RX[1] < x) RX[1] = x;

	    if(y < 0.0) continue;
	  
            if(RY[ipar][0] > y) RY[ipar][0] = y;
            if(RY[ipar][1] < y) RY[ipar][1] = y;
          }
        }
        
        if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
          for(int i=0;i<gr_Nhits_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	    double x,y;
	    gr_Nhits_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	    if(RNhits_vs_theta_phi_mom[0] > y) RNhits_vs_theta_phi_mom[0] = y;
	    if(RNhits_vs_theta_phi_mom[1] < y) RNhits_vs_theta_phi_mom[1] = y;
	  
	    gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	    if(RMatBudget_vs_theta_phi_mom[0] > y) RMatBudget_vs_theta_phi_mom[0] = y;
	    if(RMatBudget_vs_theta_phi_mom[1] < y) RMatBudget_vs_theta_phi_mom[1] = y;
	  }
        }
        
      }
    }
  }

  double min_frac = 0.6;
  double porcent,delta,min;
  delta = RX[1] - RX[0];
  porcent = 0.1;
  if(TMath::Abs(delta*global->GetUnit(MonUnits)) < 1.0e-6*global->GetUnit(MonUnits)) {
    delta   = 1.0;
    porcent = 1.0;
  }
  RX[0] -= porcent*delta;
  RX[1] += porcent*delta;
  if(RX[0] < 0.0) RX[0] = 0.0;
  
  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    delta = RY[ipar][1] - RY[ipar][0];
    if(TMath::Abs(delta) < 1.0e-6) delta = RY[ipar][0];
    porcent = 0.1;
    TString ParErrorUnit = TrackerList[0]->GetTrajectory()->GetParameterErrorUnit(ipar);
    if(TMath::Abs(delta*global->GetUnit(ParErrorUnit)) < 1.0e-6*global->GetUnit(ParErrorUnit)) {
      delta   = 1.0;
      porcent = 1.0;
    }
    
    min = RY[ipar][0];
    RY[ipar][0] -= porcent*delta;
    RY[ipar][1] += porcent*delta;
    if(RY[ipar][0] < 0.0) RY[ipar][0] = min*min_frac;
  }
  
  if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
    min = RNhits_vs_theta_phi_mom[0];
    delta = RNhits_vs_theta_phi_mom[1] - RNhits_vs_theta_phi_mom[0];
    porcent = 0.1;
    RNhits_vs_theta_phi_mom[0] -= porcent*delta;
    RNhits_vs_theta_phi_mom[1] += porcent*delta;
    if(RNhits_vs_theta_phi_mom[0] < 0.0) RNhits_vs_theta_phi_mom[0] = min*min_frac;
  
    min = RMatBudget_vs_theta_phi_mom[0];
    delta = RMatBudget_vs_theta_phi_mom[1] - RMatBudget_vs_theta_phi_mom[0];
    porcent = 0.1;
    RMatBudget_vs_theta_phi_mom[0] -= porcent*delta;
    RMatBudget_vs_theta_phi_mom[1] += porcent*delta;
    if(RMatBudget_vs_theta_phi_mom[0] < 0.0) RMatBudget_vs_theta_phi_mom[0] = min*min_frac;
  }

  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    if(!GoodRangeParam[ipar]) {
      RY[ipar][0] =  1.0e-6;
      RY[ipar][1] =  1.0;
    }
  }


  if(DoTelescopeAnalysis) { //begin if doing telescope analysis
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
      for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	NDUTS[itheta][iphi] = 0;
	
	//Getting the list of systems for a given theta and phi
	std::vector<TString> LayersList;
	LayersList.clear();
	for(int imom=0;imom<NmomPoints;imom++) {
	  for(int idut=0;idut<MaxDUTs;idut++) {
	    for(int igeo=0;igeo<NGeometries;igeo++) {
	      if(DUTLayerNames[igeo][itheta][iphi][imom][idut] == TString("")) continue;
	      bool IsinList = false;
	      for(int kkk=0;kkk<int(LayersList.size());kkk++) {
	        if(LayersList[kkk] == DUTLayerNames[igeo][itheta][iphi][imom][idut]) {
		  IsinList = true;
		  break;
	        }
	      }
	      if(!IsinList) LayersList.push_back(DUTLayerNames[igeo][itheta][iphi][imom][idut]);
	    }
	  }
	}
	
	NDUTS[itheta][iphi] = LayersList.size();
	for(int kkk=0;kkk<int(LayersList.size());kkk++) { //begin systems loop
	  DUTNames[itheta][iphi][kkk] = LayersList[kkk];
	  
	  for(int igeo=0;igeo<NGeometries;igeo++) { //begin of loop over geometries
	    HistTitle  = TString(gr_TelResolUAtDUT[igeo][itheta][iphi][kkk]->GetTitle());
	    HistTitle += TString(" ") + LayersList[kkk];
	    gr_TelResolUAtDUT[igeo][itheta][iphi][kkk]->SetTitle(HistTitle.Data());
	    
	    HistTitle  = TString(gr_TelResolVAtDUT[igeo][itheta][iphi][kkk]->GetTitle());
	    HistTitle += TString(" ") + LayersList[kkk];
	    gr_TelResolVAtDUT[igeo][itheta][iphi][kkk]->SetTitle(HistTitle.Data());
	    
	    HistTitle  = TString(gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][kkk]->GetTitle());
	    HistTitle += TString(" ") + LayersList[kkk];
	    gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][kkk]->SetTitle(HistTitle.Data());
	    
	    HistTitle  = TString(gr_NbkgAtDUT[igeo][itheta][iphi][kkk]->GetTitle());
	    HistTitle += TString(" ") + LayersList[kkk];
	    gr_NbkgAtDUT[igeo][itheta][iphi][kkk]->SetTitle(HistTitle.Data());

	    std::vector<float> mom_values;
	    std::vector<float> TelResolU_values;
	    std::vector<float> TelResolV_values;
	    std::vector<float> TelResolUVArea_values;
	    std::vector<float> Nbkg_values;
	    mom_values.clear();
	    TelResolU_values.clear();
	    TelResolV_values.clear();
	    TelResolUVArea_values.clear();
	    Nbkg_values.clear();
	    for(int idut=0;idut<MaxDUTs;idut++) {
	      for(int imom=0;imom<NmomPoints;imom++) {
		if(LayersList[kkk] != DUTLayerNames[igeo][itheta][iphi][imom][idut]) continue;
		if(MomValue[igeo][itheta][iphi][imom][idut] == Dummy_value) continue;

		mom_values.push_back(MomValue[igeo][itheta][iphi][imom][idut]);
		TelResolU_values.push_back(TelResolUAtDUT[igeo][itheta][iphi][imom][idut]);
		TelResolV_values.push_back(TelResolVAtDUT[igeo][itheta][iphi][imom][idut]);
		TelResolUVArea_values.push_back(TelResolUVAreaAtDUT[igeo][itheta][iphi][imom][idut]);
		Nbkg_values.push_back(NbkgAtDUT[igeo][itheta][iphi][imom][idut]);
	      }
	    }
	    
	    for(int imom=0;imom<int(mom_values.size());imom++) {
	      gr_TelResolUAtDUT[igeo][itheta][iphi][kkk]->SetPoint(imom,mom_values[imom],TelResolU_values[imom]);
	      gr_TelResolVAtDUT[igeo][itheta][iphi][kkk]->SetPoint(imom,mom_values[imom],TelResolV_values[imom]);
	      gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][kkk]->SetPoint(imom,mom_values[imom],TelResolUVArea_values[imom]);
	      gr_NbkgAtDUT[igeo][itheta][iphi][kkk]->SetPoint(imom,mom_values[imom],Nbkg_values[imom]);
	    }
	    
	  } //end of loop over geometries
	} //end systems loop
	
      } //end of phi loop
    } // end of theta loop
    
  } // end if doing telescope analysis
  
  
  //Calculating the average of the resolution on track parameters vs phi
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	std::vector<double> momList;
	std::vector<double> ParList;
	std::vector<double> Par2List;
	std::vector<int>    pointsList;
	momList.clear();
	ParList.clear();
	Par2List.clear();
	pointsList.clear();
	for(int iphi=0;iphi<NphiPoints;iphi++) {
	  for(int i=0;i<gr[igeo][itheta][iphi][ipar]->GetN();i++) {
	    double x,y;
	    gr[igeo][itheta][iphi][ipar]->GetPoint(i,x,y);
	    //if(y < 0) continue;
	    
	    bool IsInlist = false;
	    for(int imom=0;imom<int(momList.size());imom++) {
	      if(TMath::Abs(x - momList[imom]) < 1.0e-6) {
		IsInlist = true;
		break;
	      }
	    }
	    if(!IsInlist) {
	      momList.push_back(x);
	      ParList.push_back(0);
	      Par2List.push_back(0);
	      pointsList.push_back(0);
	    }
	  }
	}
	
	int counter_mom = 0;
	for(int imom=0;imom<int(momList.size());imom++) {
	  for(int iphi=0;iphi<NphiPoints;iphi++) {
	    double x,y;
	    for(int i=0;i<gr[igeo][itheta][iphi][ipar]->GetN();i++) {
	      gr[igeo][itheta][iphi][ipar]->GetPoint(i,x,y);
	      if(TMath::Abs(momList[imom] - x) < 1.0e-6) break;
	    }
	    if(x == -999.0) continue;
	    if(y < 0.0) continue;
	    ParList[imom]  += y;
	    Par2List[imom] += pow(y,2);
	    pointsList[imom]++;
	  }
	  if(pointsList[imom] > 0) {
	    ParList[imom]  /= pointsList[imom];
	    Par2List[imom] /= pointsList[imom];
	    Par2List[imom] -= pow(ParList[imom],2);
	    Par2List[imom]  = sqrt(Par2List[imom]);
	  }
	  else {
	    ParList[imom] = -1.0;
	    Par2List[imom] = 1.0e-20;
	  }
	  
	  gr_aveVsPhi[igeo][itheta][ipar]->SetPoint(counter_mom,momList[imom],ParList[imom]);
	  gr_aveVsPhi[igeo][itheta][ipar]->SetPointError(counter_mom,1.0e-8,Par2List[imom]/sqrt(pointsList[imom]));
	  
	  gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetPoint(counter_mom,momList[imom],ParList[imom]);
	  gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->SetPointError(counter_mom,1.0e-8,Par2List[imom]);
	  
	  gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetPoint(counter_mom,momList[imom],ParList[imom]);
	  
	  counter_mom++;
	}
	
      }
    }
  }

  double Rdist1stPoint[2];
  Rdist1stPoint[0] = +1.0e+20;
  Rdist1stPoint[1] = -1.0e+20;
  if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
    for(int igeo=0;igeo<NGeometries;igeo++) {
      for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
        for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	  for(int iP=0; iP<gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetN(); iP++) {
	    double x,y;
	    gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetPoint(iP, x, y);
	    if(Rdist1stPoint[0] > y) Rdist1stPoint[0] = y;
	    if(Rdist1stPoint[1] < y) Rdist1stPoint[1] = y;
          }
        }
      }
    }
  }
  
  if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
    //Calculate average value of distance of 1st point to track origin vs phi
    for(int igeo=0;igeo<NGeometries;igeo++) {
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
        double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
        if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
        else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
        else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
      
        std::vector<double> momList;
        std::vector<double> ParList;
        std::vector<int>    pointsList;
        momList.clear();
        ParList.clear();
        pointsList.clear();
        for(int iphi=0;iphi<NphiPoints;iphi++) {
	  for(int i=0;i<gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	    double x,y;
	    gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	  
	    bool IsInlist = false;
	    for(int imom=0;imom<int(momList.size());imom++) {
	      if(TMath::Abs(x - momList[imom]) < 1.0e-6) {
	        IsInlist = true;
	        break;
	      }
	    }
	    if(!IsInlist) {
	      momList.push_back(x);
	      ParList.push_back(0);
	      pointsList.push_back(0);
	    }
  	  }
        }
      
        int counter_mom = 0;
        for(int imom=0;imom<int(momList.size());imom++) {
	  for(int iphi=0;iphi<NphiPoints;iphi++) {
	    double x,y;
	    for(int i=0;i<gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(TMath::Abs(momList[imom] - x) < 1.0e-6) break;
	    }
	    if(x == -999.0) continue;
	    ParList[imom]  += y;
	    pointsList[imom]++;
	  }
	  if(pointsList[imom] == 0) continue;
	  ParList[imom]  /= pointsList[imom];
	
	  gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetPoint(counter_mom,momList[imom],ParList[imom]);
	  counter_mom++;
        }
      
        double AveMon = 0.0;
        for(int i=0;i<gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetN();i++) {
	  double x,y;
	  gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetPoint(i,x,y);
	  AveMon += y;
        }
        if(gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetN() > 0) AveMon /= gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetN();
  
        gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->SetPoint(itheta,thetaValue,AveMon);
      }
    }
  }

  //Estimate the a and b parameters for each geometry
  double Rtheta[2];
  double Rphi[2];
  double Ra[2];
  double Ra_vs_phi[2];
  double Rb[2];
  double RNhits[2];
  Rtheta[0] = Rphi[0] = Ra[0] = Ra_vs_phi[0] = Rb[0] = RNhits[0] = +1.0e+20;
  Rtheta[1] = Rphi[1] = Ra[1] = Ra_vs_phi[1] = Rb[1] = RNhits[1] = -1.0e+20;
  double a_param[NGeometries][NthetaPoints];
  double b_param[NGeometries][NthetaPoints];
  double exp_param[NGeometries][NthetaPoints];
  double a_param_ave[NGeometries];
  double b_param_ave[NGeometries];
  double a_param_ave_err[NGeometries];
  double b_param_ave_err[NGeometries];
  double exp_param_ave[NGeometries];
  int par_idx = TrackerList[0]->GetTrajectory()->GetTheDOCAParIndex();

  for(int itheta=0;itheta<NthetaPoints;itheta++) { //begin loop over theta
    double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
    if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
    else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
    else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
    
    if(Rtheta[0] > thetaValue) Rtheta[0] = thetaValue;
    if(Rtheta[1] < thetaValue) Rtheta[1] = thetaValue;
  }
  
  //Calculate a and b parameters for each geometry and for each theta value
  //Calculate as well the a and b average value for all theta
  if(TrkResolAnalysisPars[0].PlotDOCAatHighMom || TrkResolAnalysisPars[0].PlotMaterialBudget) { //begin if do DOCA or MatBudget
    for(int igeo=0;igeo<NGeometries;igeo++) { //begin loop over geometries
      a_param_ave[igeo]      = 0.0;
      b_param_ave[igeo]      = 0.0;
      int counter_theta_a    = 0;
      int counter_theta_b    = 0;
      for(int itheta=0;itheta<NthetaPoints;itheta++) { //begin loop over theta
        double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
        if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
        else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
        else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
      
        a_param[igeo][itheta]   = 0.0;
        b_param[igeo][itheta]   = 0.0;
        exp_param[igeo][itheta] = 0.0;
      
        int Ngood_phi_points_a = 0;
        for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	  double mon_ppp = 1000.0*global->GetMomentumUnit("GeV/c");
	  TVector3 momentum =  global->GetMomentum(mon_ppp,thetaArr[itheta],phiArr[iphi]);
	  if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
	    if(mon_ppp > momArr[momArr.size()-1]) {
	      TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);
	      TMatrixD  FitCovMatrix;
	      bool GoodFit = TrackerList[igeo]->GetFitTrackParsCovMatrix(TrkResolAnalysisPars[0].NhitsMin,FitCovMatrix);
	  
	      if(GoodFit) {
	        double the_a = sqrt(FitCovMatrix(par_idx,par_idx));
	        a_param[igeo][itheta] += the_a;
	        Ngood_phi_points_a++;
	  
	        if(Ra_vs_phi[0] > the_a) Ra_vs_phi[0] = the_a;
                if(Ra_vs_phi[1] < the_a) Ra_vs_phi[1] = the_a;
	        if(Ra[0]        > the_a) Ra[0]        = the_a;
                if(Ra[1]        < the_a) Ra[1]        = the_a;
	  
	        gr_a_vs_Phi[igeo][itheta]->SetPoint(iphi,phiArr[iphi]/global->GetAngleUnit("deg"),the_a/global->GetDistanceUnit("um"));
	      }
	    }
	    else {
	      int Npoints_gr = gr[igeo][itheta][iphi][par_idx]->GetN();
	      if(Npoints_gr > 0) {
	        double x,the_a;
	        gr[igeo][itheta][iphi][par_idx]->GetPoint(Npoints_gr-1,x,the_a);
	        the_a *= global->GetDistanceUnit("um");
	  
	        a_param[igeo][itheta] += the_a;
	        Ngood_phi_points_a++;
	  
	        if(Ra_vs_phi[0] > the_a) Ra_vs_phi[0] = the_a;
                if(Ra_vs_phi[1] < the_a) Ra_vs_phi[1] = the_a;
	        if(Ra[0]        > the_a) Ra[0]        = the_a;
                if(Ra[1]        < the_a) Ra[1]        = the_a;
	  
	        gr_a_vs_Phi[igeo][itheta]->SetPoint(iphi,phiArr[iphi]/global->GetAngleUnit("deg"),the_a/global->GetDistanceUnit("um"));
	      }
	    }
	  }

	  if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
	    int Nhits = 0;
	    if(mon_ppp > momArr[momArr.size()-1]) {
	      TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);
	      double  Material_budget   = 0.0;
	      double  dist1stPointToRef = 0.0;
	      TrackerList[igeo]->GetMaterialBudget(Material_budget,Nhits,dist1stPointToRef);
	    }
	    else {
	      int Npoints_gr = gr_Nhits_vs_mom[igeo][itheta][iphi]->GetN();
	      if(Npoints_gr > 0) {
	        double x,y;
	        gr_Nhits_vs_mom[igeo][itheta][iphi]->GetPoint(Npoints_gr-1,x,y);
	        Nhits = int(y);
	      }
	    }

	    if(RNhits[0] > Nhits) RNhits[0] = Nhits;
            if(RNhits[1] < Nhits) RNhits[1] = Nhits;
	    
	    gr_Nhits_vs_Phi[igeo][itheta]->SetPoint(iphi,phiArr[iphi]/global->GetAngleUnit("deg"),Nhits);
	  }
	
	  if(Rphi[0] > phiArr[iphi]) Rphi[0] = phiArr[iphi];
          if(Rphi[1] < phiArr[iphi]) Rphi[1] = phiArr[iphi];
        } //end of loop over phi values

        if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
          if(Ngood_phi_points_a > 0) {
	    a_param[igeo][itheta] /= Ngood_phi_points_a;
	    if(Ra[0] > a_param[igeo][itheta]) Ra[0] = a_param[igeo][itheta];
            if(Ra[1] < a_param[igeo][itheta]) Ra[1] = a_param[igeo][itheta];
	    gr_a[igeo]->SetPoint(counter_theta_a,thetaValue,a_param[igeo][itheta]/global->GetDistanceUnit("um"));
	    a_param_ave[igeo]   += a_param[igeo][itheta];
	    counter_theta_a++;
          }

          if(gr_aveVsPhi[igeo][itheta][par_idx]->GetN() >= 5) {
            double a;
            global->GetGeometryImpacParameterParameters(gr_aveVsPhi[igeo][itheta][par_idx],thetaArr[itheta],ParamUnits[par_idx],MonUnits,FitPowerForImpactParam,a,b_param[igeo][itheta],exp_param[igeo][itheta]);
      
            if(Rb[0] > b_param[igeo][itheta]) Rb[0] = b_param[igeo][itheta];
            if(Rb[1] < b_param[igeo][itheta]) Rb[1] = b_param[igeo][itheta];

            gr_b[igeo]->SetPoint(counter_theta_b,thetaValue,b_param[igeo][itheta]/(global->GetDistanceUnit("um")*global->GetMomentumUnit("GeV/c")));

            b_param_ave[igeo]   += b_param[igeo][itheta];
            exp_param_ave[igeo] += exp_param[igeo][itheta];
	    counter_theta_b++;
          }
        }
        
      } //end loop over theta

      if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
        if(counter_theta_a > 0) a_param_ave[igeo]   /= counter_theta_a;
        if(counter_theta_b > 0) {
          b_param_ave[igeo]   /= counter_theta_b;
          exp_param_ave[igeo] /= counter_theta_b;
        }
      }

    } //end loop over geometries
  } //end if do DOCA or MatBudget

  if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
    for(int igeo=0;igeo<NGeometries;igeo++) {
      if(NthetaPoints == 1) {
        a_param_ave_err[igeo] = 0.0;
        b_param_ave_err[igeo] = 0.0;
      }
      else {
        a_param_ave_err[igeo] = -1.0e+20;
        b_param_ave_err[igeo] = -1.0e+20;
        for(int itheta=0;itheta<NthetaPoints;itheta++) {
	  for(int jtheta=0;jtheta<NthetaPoints;jtheta++) {
	    if(itheta >= jtheta) continue;
	  
	    double delta_a = TMath::Abs(a_param[igeo][itheta] - a_param[igeo][jtheta]);
	    double delta_b = TMath::Abs(b_param[igeo][itheta] - b_param[igeo][jtheta]);
	  
	    if(a_param_ave_err[igeo] < delta_a) a_param_ave_err[igeo] = delta_a;
	    if(b_param_ave_err[igeo] < delta_b) b_param_ave_err[igeo] = delta_b;
	  }
        }
      }
    
      a_param_ave_err[igeo] *= 0.5;
      b_param_ave_err[igeo] *= 0.5;
    }
  }

  double RNhits_vs_theta[2];
  RNhits_vs_theta[0] = +1.0e+20;
  RNhits_vs_theta[1] = -1.0e+20;
  if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
    for(int igeo=0;igeo<NGeometries;igeo++) {
      for(int iP=0; iP<MyNbins_mom_redu; iP++) {
        double mom;
        if(TrkResolAnalysisPars[0].UseAllMomVals) mom = momArr[iP];
        else                                      mom = momArr[0] + (iP+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
      
        for(int itheta=0;itheta<NthetaPoints;itheta++) {
	  double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
          if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
          else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
          else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
      
	  double Nhits_ave = 0;
	  int    Npoints   = 0;
	  for(int iphi=0;iphi<NphiPoints;iphi++) {
	    TVector3 momentum = global->GetMomentum(mom,thetaArr[itheta],phiArr[iphi]);
	    TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);
	  
	    int     Nhits             = 0;
	    double  Material_budget   = 0.0;
	    double  dist1stPointToRef = 0.0;
	    TrackerList[igeo]->GetMaterialBudget(Material_budget,Nhits,dist1stPointToRef);
	  
	    Nhits_ave += Nhits;
	    Npoints++;
	  }
	  if(Npoints > 0) Nhits_ave /= Npoints;
	
	  gr_aveNhitsVsTheta[igeo][iP]->SetPoint(itheta,thetaValue,Nhits_ave);
	
	  if(RNhits_vs_theta[0] > Nhits_ave) RNhits_vs_theta[0] = Nhits_ave;
	  if(RNhits_vs_theta[1] < Nhits_ave) RNhits_vs_theta[1] = Nhits_ave;
        }
      }
    }
  }

  double Rpar_vs_theta[Nmax_params][2];
  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters() ;ipar++) {
    Rpar_vs_theta[ipar][0] = +1.0e+20;
    Rpar_vs_theta[ipar][1] = -1.0e+20;
  }
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
      for(int iP=0; iP<MyNbins_mom_redu; iP++) {
	double pi;
	if(TrkResolAnalysisPars[0].UseAllMomVals) pi = momArr[iP];
	else                                      pi = momArr[0] + (iP+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
	pi /= global->GetUnit(MonUnits);
	int counter_theta = 0;
	for(int itheta=0;itheta<NthetaPoints;itheta++) {
	  double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
          if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
          else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
          else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
	  
	  //if(gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetN() == 0) continue;
	  int idx = -999;
	  double mom1,mom2,y1,y2;
	  for(int i=0;i<gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetN()-1;i++) {
	    gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetPoint(i,  mom1,y1);
	    gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetPoint(i+1,mom2,y2);
	    if(pi >= mom1 && pi <= mom2) {
	      idx = i;
	      break;
	    }
	  }
	  if(idx == -999) continue;
	  
	  double a = (y2-y1)/(mom2 - mom1);
	  double b = y2 - a*mom2;
	  
	  double value = a*pi + b;
	  gr_aveSigmaParVsTheta[igeo][iP][ipar]->SetPoint(counter_theta,thetaValue,value);
	  if(Rpar_vs_theta[ipar][0] > value) Rpar_vs_theta[ipar][0] = value;
	  if(Rpar_vs_theta[ipar][1] < value) Rpar_vs_theta[ipar][1] = value;
	  
	  counter_theta++;
	}
      }
    }
  }

  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    delta = Rpar_vs_theta[ipar][1] - Rpar_vs_theta[ipar][0];
    if(TMath::Abs(delta) < 1.0e-6) delta = Rpar_vs_theta[ipar][0];
    porcent = 0.1;
    
    min = Rpar_vs_theta[ipar][0];
    Rpar_vs_theta[ipar][0] -= porcent*delta;
    Rpar_vs_theta[ipar][1] += porcent*delta;
    if(Rpar_vs_theta[ipar][0] < 0.0) Rpar_vs_theta[ipar][0] = min*min_frac;
  }

  if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
    for(int igeo=0;igeo<NGeometries;igeo++) {
      cout << endl;
      cout << "===========================================================================================" << endl;
      cout << "Printing out a and b parameters for geometry " << GeometryList[igeo]->GetName().Data() << " :" << endl;
      cout << "===========================================================================================" << endl;
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
	char thetaTitle[100];
        if(polarVariable == TString("theta"))         sprintf(thetaTitle,"theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
        else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
        else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
	
        cout << " - a = " << a_param[igeo][itheta]/global->GetDistanceUnit("um") << " um, b = " << b_param[igeo][itheta]/(global->GetDistanceUnit("um")*global->GetMomentumUnit("GeV/c")) << " um x GeV/c, "
             << "exp = " << exp_param[igeo][itheta] << ", "
	     << "for " << thetaTitle << " "
	     << endl;
      }
      cout << "-------------------------------------------------------------------------------------------" << endl;
      cout << "-------------------------------------------------------------------------------------------" << endl;
      cout << "a = (" << a_param_ave[igeo]/global->GetDistanceUnit("um") << " +/- " << a_param_ave_err[igeo]/global->GetDistanceUnit("um") << ") um, " 
           << "b = (" << b_param_ave[igeo]/(global->GetDistanceUnit("um")*global->GetMomentumUnit("GeV/c")) << " +/- " << b_param_ave_err[igeo]/(global->GetDistanceUnit("um")*global->GetMomentumUnit("GeV/c")) << ") um x GeV/c, "
           << "Average for all " << polarVariable.Data() << " values" 
	   << endl;
      cout << "===========================================================================================" << endl;
      cout << endl;
    }
  }
  
  delta = Rtheta[1] - Rtheta[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Rtheta[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = 1.0;
  porcent = 0.1;
  Rtheta[0] -= porcent*delta;
  Rtheta[1] += porcent*delta;
  if(Rtheta[0] < 0.0) Rtheta[0] = 0.0;
  
  delta = Rphi[1] - Rphi[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Rphi[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = 1.0;
  porcent = 0.1;
  Rphi[0] -= porcent*delta;
  Rphi[1] += porcent*delta;
  //if(Rphi[0] < 0.0) Rphi[0] = 0.0;
  
  delta = Ra[1] - Ra[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Ra[0];
  porcent = 0.1;
  min = Ra[0];
  Ra[0] -= porcent*delta;
  Ra[1] += porcent*delta;
  if(Ra[0] < 0.0) Ra[0] = min*min_frac;
  
  delta = Ra_vs_phi[1] - Ra_vs_phi[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Ra_vs_phi[0];
  porcent = 0.1;
  min = Ra_vs_phi[0];
  Ra_vs_phi[0] -= porcent*delta;
  Ra_vs_phi[1] += porcent*delta;
  if(Ra_vs_phi[0] < 0.0) Ra_vs_phi[0] = min*min_frac;
  
  delta = Rb[1] - Rb[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Rb[0];
  porcent = 0.1;
  min = Rb[0];
  Rb[0] -= porcent*delta;
  Rb[1] += porcent*delta;
  if(Rb[0] < 0.0) Rb[0] = min*min_frac;
  
  delta = RNhits[1] - RNhits[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = RNhits[0];
  porcent = 0.1;
  min = RNhits[0];
  RNhits[0] -= porcent*delta;
  RNhits[1] += porcent*delta;
  if(RNhits[0] < 0.0) RNhits[0] = min*min_frac;
  
  delta = RNhits_vs_theta[1] - RNhits_vs_theta[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = RNhits_vs_theta[0];
  porcent = 0.1;
  min = RNhits_vs_theta[0];
  RNhits_vs_theta[0] -= porcent*delta;
  RNhits_vs_theta[1] += porcent*delta;
  if(RNhits_vs_theta[0] < 0.0) RNhits_vs_theta[0] = min*min_frac;
  
  Rphi[0]      /= global->GetAngleUnit("deg");
  Rphi[1]      /= global->GetAngleUnit("deg");
  Ra[0]        /= global->GetDistanceUnit("um");
  Ra[1]        /= global->GetDistanceUnit("um");
  Ra_vs_phi[0] /= global->GetDistanceUnit("um");
  Ra_vs_phi[1] /= global->GetDistanceUnit("um");
  Rb[0]        /= global->GetDistanceUnit("um")*global->GetMomentumUnit("GeV/c");
  Rb[1]        /= global->GetDistanceUnit("um")*global->GetMomentumUnit("GeV/c");

  TString MyPolarVarName;
  TString MyPolarVarNameWithUnits;
  if(polarVariable == TString("theta")) {
    MyPolarVarName          = TString("#theta");
    MyPolarVarNameWithUnits = TString("#theta (") + ThetaUnits + TString(")");
  }
  else if(polarVariable == TString("costheta")) {
    MyPolarVarName          = TString("cos(#theta)");
    MyPolarVarNameWithUnits = TString("cos(#theta)");
  }
  else if(polarVariable == TString("eta")) {
    MyPolarVarName          = TString("#eta");
    MyPolarVarNameWithUnits = TString("#eta");
  }
  
  TH1F* href[Nmax_params];
  TH1F* href_vs_theta[Nmax_params];
  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    HistName = TString("href_par") + long(ipar);
    href[ipar] = new TH1F(HistName.Data(),title[ipar].Data(),
			  100,RX[0],RX[1]);
    HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
    href[ipar]->SetXTitle(HistName.Data());
    href[ipar]->GetXaxis()->CenterTitle(true);
    href[ipar]->SetYTitle(titleY[ipar].Data());
    href[ipar]->GetYaxis()->CenterTitle(true);
    href[ipar]->GetYaxis()->SetTitleOffset(TitleOffSet);
    href[ipar]->SetLineColor(1);
    href[ipar]->SetLineWidth(1);
    href[ipar]->SetMinimum(RY[ipar][0]);
    href[ipar]->SetMaximum(RY[ipar][1]);
    href[ipar]->SetNdivisions(5 + 100*5,"X");
    href[ipar]->SetNdivisions(5 + 100*5,"Y");
    href[ipar]->GetXaxis()->SetTitleSize(TheSize);
    href[ipar]->GetXaxis()->SetLabelSize(TheSize);
    href[ipar]->GetYaxis()->SetTitleSize(TheSize);
    href[ipar]->GetYaxis()->SetLabelSize(TheSize);
    href[ipar]->SetStats(false);
    
    HistName  = TString("href_par") + long(ipar) + TString("_vs_theta");
    HistTitle = titleY[ipar] + TString(" vs ") + MyPolarVarName;
    href_vs_theta[ipar] = new TH1F(HistName.Data(),
				   HistTitle.Data(),
				   100,Rtheta[0],Rtheta[1]);
    HistName = MyPolarVarNameWithUnits;
    href_vs_theta[ipar]->SetXTitle(HistName.Data());
    HistName = TString("#phi Ave ") + titleY[ipar];
    href_vs_theta[ipar]->GetXaxis()->CenterTitle(true);
    href_vs_theta[ipar]->SetYTitle(HistName.Data());
    href_vs_theta[ipar]->GetYaxis()->CenterTitle(true);
    href_vs_theta[ipar]->GetYaxis()->SetTitleOffset(TitleOffSet);
    href_vs_theta[ipar]->SetLineColor(1);
    href_vs_theta[ipar]->SetLineWidth(1);
    href_vs_theta[ipar]->SetMinimum(Rpar_vs_theta[ipar][0]);
    href_vs_theta[ipar]->SetMaximum(Rpar_vs_theta[ipar][1]);
    href_vs_theta[ipar]->SetNdivisions(5 + 100*5,"X");
    href_vs_theta[ipar]->SetNdivisions(5 + 100*5,"Y");
    href_vs_theta[ipar]->GetXaxis()->SetTitleSize(TheSize);
    href_vs_theta[ipar]->GetXaxis()->SetLabelSize(TheSize);
    href_vs_theta[ipar]->GetYaxis()->SetTitleSize(TheSize);
    href_vs_theta[ipar]->GetYaxis()->SetLabelSize(TheSize);
    href_vs_theta[ipar]->SetStats(false);
  }
  
  //Reference hist for a and b parameters
  HistName  = TString("href_a");
  HistTitle = TString("a parameters vs ") + MyPolarVarName;
  TH1F* href_a = new TH1F(HistName.Data(),
			  HistTitle.Data(),
			  100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_a->SetXTitle(HistName.Data());
  href_a->GetXaxis()->CenterTitle(true);
  href_a->SetYTitle("a (#mum)");
  href_a->GetYaxis()->CenterTitle(true);
  href_a->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_a->SetLineColor(1);
  href_a->SetLineWidth(1);
  href_a->SetMinimum(Ra[0]);
  href_a->SetMaximum(Ra[1]);
  href_a->SetNdivisions(5 + 100*5,"X");
  href_a->SetNdivisions(5 + 100*5,"Y");
  href_a->GetXaxis()->SetTitleSize(TheSize);
  href_a->GetXaxis()->SetLabelSize(TheSize);
  href_a->GetYaxis()->SetTitleSize(TheSize);
  href_a->GetYaxis()->SetLabelSize(TheSize);
  href_a->SetStats(false);
  
  HistName  = TString("href_a_vs_phi");
  HistTitle = TString("a parameters vs #phi");
  TH1F* href_a_vs_phi = new TH1F(HistName.Data(),HistTitle.Data(),
				 100,Rphi[0],Rphi[1]);
  HistName = TString("#phi (") + ThetaUnits + TString(")");
  href_a_vs_phi->SetXTitle(HistName.Data());
  href_a_vs_phi->GetXaxis()->CenterTitle(true);
  href_a_vs_phi->SetYTitle("a (#mum)");
  href_a_vs_phi->GetYaxis()->CenterTitle(true);
  href_a_vs_phi->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_a_vs_phi->SetLineColor(1);
  href_a_vs_phi->SetLineWidth(1);
  href_a_vs_phi->SetMinimum(Ra[0]);
  href_a_vs_phi->SetMaximum(Ra[1]);
  href_a_vs_phi->SetNdivisions(5 + 100*5,"X");
  href_a_vs_phi->SetNdivisions(5 + 100*5,"Y");
  href_a_vs_phi->GetXaxis()->SetTitleSize(TheSize);
  href_a_vs_phi->GetXaxis()->SetLabelSize(TheSize);
  href_a_vs_phi->GetYaxis()->SetTitleSize(TheSize);
  href_a_vs_phi->GetYaxis()->SetLabelSize(TheSize);
  href_a_vs_phi->SetStats(false);
  
  HistName = TString("href_Nhits_vs_phi");
  TH1F* href_Nhits_vs_phi = new TH1F(HistName.Data(),
				     "N_{hits} vs #phi",
				     100,Rphi[0],Rphi[1]);
  HistName = TString("#phi (") + ThetaUnits + TString(")");
  href_Nhits_vs_phi->SetXTitle(HistName.Data());
  href_Nhits_vs_phi->GetXaxis()->CenterTitle(true);
  href_Nhits_vs_phi->SetYTitle("N_{hits}");
  href_Nhits_vs_phi->GetYaxis()->CenterTitle(true);
  href_Nhits_vs_phi->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Nhits_vs_phi->SetLineColor(1);
  href_Nhits_vs_phi->SetLineWidth(1);
  href_Nhits_vs_phi->SetMinimum(RNhits[0]);
  href_Nhits_vs_phi->SetMaximum(RNhits[1]);
  href_Nhits_vs_phi->SetNdivisions(5 + 100*5,"X");
  href_Nhits_vs_phi->SetNdivisions(5 + 100*5,"Y");
  href_Nhits_vs_phi->GetXaxis()->SetTitleSize(TheSize);
  href_Nhits_vs_phi->GetXaxis()->SetLabelSize(TheSize);
  href_Nhits_vs_phi->GetYaxis()->SetTitleSize(TheSize);
  href_Nhits_vs_phi->GetYaxis()->SetLabelSize(TheSize);
  href_Nhits_vs_phi->SetStats(false);
  
  HistName  = TString("href_AveNhits_vs_theta");
  HistTitle = TString("#phi Ave N_{hits} vs ") + MyPolarVarName;
  TH1F* href_AveNhits_vs_theta = new TH1F(HistName.Data(),
					  HistTitle.Data(),
					  100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_AveNhits_vs_theta->SetXTitle(HistName.Data());
  href_AveNhits_vs_theta->GetXaxis()->CenterTitle(true);
  href_AveNhits_vs_theta->SetYTitle("#phi Ave N_{hits}");
  href_AveNhits_vs_theta->GetYaxis()->CenterTitle(true);
  href_AveNhits_vs_theta->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_AveNhits_vs_theta->SetLineColor(1);
  href_AveNhits_vs_theta->SetLineWidth(1);
  href_AveNhits_vs_theta->SetMinimum(RNhits_vs_theta[0]);
  href_AveNhits_vs_theta->SetMaximum(RNhits_vs_theta[1]);
  href_AveNhits_vs_theta->SetNdivisions(5 + 100*5,"X");
  href_AveNhits_vs_theta->SetNdivisions(5 + 100*5,"Y");
  href_AveNhits_vs_theta->GetXaxis()->SetTitleSize(TheSize);
  href_AveNhits_vs_theta->GetXaxis()->SetLabelSize(TheSize);
  href_AveNhits_vs_theta->GetYaxis()->SetTitleSize(TheSize);
  href_AveNhits_vs_theta->GetYaxis()->SetLabelSize(TheSize);
  href_AveNhits_vs_theta->SetStats(false);
  
  HistName  = TString("href_b");
  HistTitle = TString("b parameters vs ") + MyPolarVarName;
  TH1F* href_b = new TH1F(HistName.Data(),HistTitle.Data(),
			  100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_b->SetXTitle(HistName.Data());
  href_b->GetXaxis()->CenterTitle(true);
  href_b->SetYTitle("b (#mum #times GeV/c)");
  href_b->GetYaxis()->CenterTitle(true);
  href_b->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_b->SetLineColor(1);
  href_b->SetLineWidth(1);
  href_b->SetMinimum(Rb[0]);
  href_b->SetMaximum(Rb[1]);
  href_b->SetNdivisions(5 + 100*5,"X");
  href_b->SetNdivisions(5 + 100*5,"Y");
  href_b->GetXaxis()->SetTitleSize(TheSize);
  href_b->GetXaxis()->SetLabelSize(TheSize);
  href_b->GetYaxis()->SetTitleSize(TheSize);
  href_b->GetYaxis()->SetLabelSize(TheSize);
  href_b->SetStats(false);
  
  HistName = TString("href_Nhits_vs_theta_phi_mom");
  TH1F* href_Nhits_vs_theta_phi_mom = new TH1F(HistName.Data(),
					       "N_{hits} vs momentum",
					       100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_Nhits_vs_theta_phi_mom->SetXTitle(HistName.Data());
  href_Nhits_vs_theta_phi_mom->GetXaxis()->CenterTitle(true);
  href_Nhits_vs_theta_phi_mom->SetYTitle("N_{hits}");
  href_Nhits_vs_theta_phi_mom->GetYaxis()->CenterTitle(true);
  href_Nhits_vs_theta_phi_mom->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Nhits_vs_theta_phi_mom->SetLineColor(1);
  href_Nhits_vs_theta_phi_mom->SetLineWidth(1);
  href_Nhits_vs_theta_phi_mom->SetMinimum(RNhits_vs_theta_phi_mom[0]);
  href_Nhits_vs_theta_phi_mom->SetMaximum(RNhits_vs_theta_phi_mom[1]);
  href_Nhits_vs_theta_phi_mom->SetNdivisions(5 + 100*5,"X");
  href_Nhits_vs_theta_phi_mom->SetNdivisions(5 + 100*5,"Y");
  href_Nhits_vs_theta_phi_mom->GetXaxis()->SetTitleSize(TheSize);
  href_Nhits_vs_theta_phi_mom->GetXaxis()->SetLabelSize(TheSize);
  href_Nhits_vs_theta_phi_mom->GetYaxis()->SetTitleSize(TheSize);
  href_Nhits_vs_theta_phi_mom->GetYaxis()->SetLabelSize(TheSize);
  href_Nhits_vs_theta_phi_mom->SetStats(false);
  
  HistName = TString("href_MatBudget_vs_theta_phi_mom");
  TH1F* href_MatBudget_vs_theta_phi_mom = new TH1F(HistName.Data(),
						   "Material budget vs momentum",
						   100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_MatBudget_vs_theta_phi_mom->SetXTitle(HistName.Data());
  href_MatBudget_vs_theta_phi_mom->GetXaxis()->CenterTitle(true);
  href_MatBudget_vs_theta_phi_mom->SetYTitle("Material Budget (%)");
  href_MatBudget_vs_theta_phi_mom->GetYaxis()->CenterTitle(true);
  href_MatBudget_vs_theta_phi_mom->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_MatBudget_vs_theta_phi_mom->SetLineColor(1);
  href_MatBudget_vs_theta_phi_mom->SetLineWidth(1);
  href_MatBudget_vs_theta_phi_mom->SetMinimum(RMatBudget_vs_theta_phi_mom[0]);
  href_MatBudget_vs_theta_phi_mom->SetMaximum(RMatBudget_vs_theta_phi_mom[1]);
  href_MatBudget_vs_theta_phi_mom->SetNdivisions(5 + 100*5,"X");
  href_MatBudget_vs_theta_phi_mom->SetNdivisions(5 + 100*5,"Y");
  href_MatBudget_vs_theta_phi_mom->GetXaxis()->SetTitleSize(TheSize);
  href_MatBudget_vs_theta_phi_mom->GetXaxis()->SetLabelSize(TheSize);
  href_MatBudget_vs_theta_phi_mom->GetYaxis()->SetTitleSize(TheSize);
  href_MatBudget_vs_theta_phi_mom->GetYaxis()->SetLabelSize(TheSize);
  href_MatBudget_vs_theta_phi_mom->SetStats(false);
  
  HistName = TString("href_1stPointDist_vs_theta_phi_mom");
  TH1F* href_1stPointDist_vs_theta_phi_mom = new TH1F(HistName.Data(),
						      "distance of 1^{st} hit to origin vs momentum",
						      100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_1stPointDist_vs_theta_phi_mom->SetXTitle(HistName.Data());
  href_1stPointDist_vs_theta_phi_mom->GetXaxis()->CenterTitle(true);
  href_1stPointDist_vs_theta_phi_mom->SetYTitle("dist of 1^{st} hit to origin (cm)");
  href_1stPointDist_vs_theta_phi_mom->GetYaxis()->CenterTitle(true);
  href_1stPointDist_vs_theta_phi_mom->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_1stPointDist_vs_theta_phi_mom->SetLineColor(1);
  href_1stPointDist_vs_theta_phi_mom->SetLineWidth(1);
  href_1stPointDist_vs_theta_phi_mom->SetMinimum(Rdist1stPoint[0]);
  href_1stPointDist_vs_theta_phi_mom->SetMaximum(Rdist1stPoint[1]);
  href_1stPointDist_vs_theta_phi_mom->SetNdivisions(5 + 100*5,"X");
  href_1stPointDist_vs_theta_phi_mom->SetNdivisions(5 + 100*5,"Y");
  href_1stPointDist_vs_theta_phi_mom->GetXaxis()->SetTitleSize(TheSize);
  href_1stPointDist_vs_theta_phi_mom->GetXaxis()->SetLabelSize(TheSize);
  href_1stPointDist_vs_theta_phi_mom->GetYaxis()->SetTitleSize(TheSize);
  href_1stPointDist_vs_theta_phi_mom->GetYaxis()->SetLabelSize(TheSize);
  href_1stPointDist_vs_theta_phi_mom->SetStats(false);
  
  HistName  = TString("href_Ave1stPointDist_vs_theta");
  HistTitle = TString("Ave distance of 1^{st} hit to origin vs ") + MyPolarVarName;
  TH1F* href_Ave1stPointDist_vs_theta = new TH1F(HistName.Data(),
						 HistTitle.Data(),
						 100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_Ave1stPointDist_vs_theta->SetXTitle(HistName.Data());
  href_Ave1stPointDist_vs_theta->GetXaxis()->CenterTitle(true);
  href_Ave1stPointDist_vs_theta->SetYTitle("dist of 1^{st} hit to origin (cm)");
  href_Ave1stPointDist_vs_theta->GetYaxis()->CenterTitle(true);
  href_Ave1stPointDist_vs_theta->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Ave1stPointDist_vs_theta->SetLineColor(1);
  href_Ave1stPointDist_vs_theta->SetLineWidth(1);
  href_Ave1stPointDist_vs_theta->SetMinimum(Rdist1stPoint[0]);
  href_Ave1stPointDist_vs_theta->SetMaximum(Rdist1stPoint[1]);
  href_Ave1stPointDist_vs_theta->SetNdivisions(5 + 100*5,"X");
  href_Ave1stPointDist_vs_theta->SetNdivisions(5 + 100*5,"Y");
  href_Ave1stPointDist_vs_theta->GetXaxis()->SetTitleSize(TheSize);
  href_Ave1stPointDist_vs_theta->GetXaxis()->SetLabelSize(TheSize);
  href_Ave1stPointDist_vs_theta->GetYaxis()->SetTitleSize(TheSize);
  href_Ave1stPointDist_vs_theta->GetYaxis()->SetLabelSize(TheSize);
  href_Ave1stPointDist_vs_theta->SetStats(false);

  int MyNhitsMin_tmp = TrkResolAnalysisPars[0].NhitsMin;
  if(TrkResolAnalysisPars[0].NhitsMin < 0) MyNhitsMin_tmp = TrackerList[0]->GetTrajectory()->GetMinHits();
  else {
    MyNhitsMin_tmp = TrackerList[0]->GetTrajectory()->GetMinHits();    
    if(TrkResolAnalysisPars[0].NhitsMin >= MyNhitsMin_tmp) MyNhitsMin_tmp = TrkResolAnalysisPars[0].NhitsMin;
  }

  TGraph* l_NhitsMin_vs_mom = new TGraph();
  l_NhitsMin_vs_mom->SetPoint(0,RX[0],MyNhitsMin_tmp);
  l_NhitsMin_vs_mom->SetPoint(1,RX[1],MyNhitsMin_tmp);
  l_NhitsMin_vs_mom->SetLineColor(kRed);
  l_NhitsMin_vs_mom->SetLineWidth(1);
  l_NhitsMin_vs_mom->SetLineStyle(2);
  
  TGraph* l_NhitsMin_vs_theta = new TGraph();
  l_NhitsMin_vs_theta->SetPoint(0,Rtheta[0],MyNhitsMin_tmp);
  l_NhitsMin_vs_theta->SetPoint(1,Rtheta[1],MyNhitsMin_tmp);
  l_NhitsMin_vs_theta->SetLineColor(kRed);
  l_NhitsMin_vs_theta->SetLineWidth(1);
  l_NhitsMin_vs_theta->SetLineStyle(2);
  
  TLegend* leg = NULL;
  leg = new TLegend(0.10,0.10,0.90,0.80);
  leg->SetFillColor(10);
  
  TString EPSName  = TheOutputFile + TString("_") + long(GlobalFileCounter) + TString(".eps");
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");

  TCanvas* c1 = new TCanvas("c1","c1",LargeCanvasX,LargeCanvasY);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.20);
  c1->SetBottomMargin(0.15);
  //c1->SetRightMargin(0.10);

  c1->Print(EPSNameO.Data());

  for(int igeo=0;igeo<NGeometries;igeo++) leg->AddEntry(gr[igeo][0][0][0],GeometryList[igeo]->GetName().Data(),"lp");
  
  TLine* l_a_ave_for_theta[NGeometries][NthetaPoints];
  for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
    char thetaTitle[300];
    if(polarVariable == TString("theta"))         sprintf(thetaTitle,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit(ThetaUnits));
    else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(#theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"#eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    for(int iphi=0;iphi<NphiPoints;iphi++) {
      char phiTitle[300];
      sprintf(phiTitle,"#phi = %.1f %s",phiArr[iphi]/global->GetUnit(ThetaUnits),ThetaUnits.Data());
      
      TLine* l_ref_mom[Nmax_params];
      
      if(!TrkResolAnalysisPars[0].PlotOnlyPhiAveraged) {
        c1->Clear();
        c1->Divide(3,2);
        int counter_pad = 0;
        for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
          counter_pad++;
          if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++; 
          c1->cd(counter_pad);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(0.20);
          gPad->SetBottomMargin(0.15);
	  gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
	  double R[2];
	  R[0] = +1.0e+20;
	  R[1] = -1.0e+20;
	  if(!TrkResolAnalysisPars[0].SameRange) {
	    for(int igeo=0;igeo<NGeometries;igeo++) {
	      for(int i=0;i<gr[igeo][itheta][iphi][ipar]->GetN();i++) {
	        double x,y;
	        gr[igeo][itheta][iphi][ipar]->GetPoint(i,x,y);
		if(y < 0.0) continue;
		
	        if(R[0] > y) R[0] = y;
	        if(R[1] < y) R[1] = y;
	      }
	    }
	    bool NoRange = false;
	    if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	    if(NoRange) {
	      R[0] = 1.0e-6;
	      R[1] = 1.0;
	    }
	    else {
	      delta = R[1] - R[0];
              if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
              porcent = 0.1;
	      min = R[0];
	      R[0] -= porcent*delta;
	      R[1] += porcent*delta;
	      if(R[0] < 0.0) R[0] = min*min_frac;
	    }
	    
	    href[ipar]->SetMinimum(R[0]);
	    href[ipar]->SetMaximum(R[1]);
	  }
	  l_ref_mom[ipar] = new TLine(momArr[0]/global->GetUnit(MonUnits),href[ipar]->GetMinimum(),
				      momArr[0]/global->GetUnit(MonUnits),href[ipar]->GetMaximum());
	  l_ref_mom[ipar]->SetLineColor(2);
	  l_ref_mom[ipar]->SetLineWidth(1);
	  l_ref_mom[ipar]->SetLineStyle(2);
	
          href[ipar]->Draw();
          for(int igeo=0;igeo<NGeometries;igeo++) {
            if(gr[igeo][itheta][iphi][ipar]->GetN() > 0) {
	      if(gr[igeo][itheta][iphi][ipar]->GetN() == 1) {
		gr[igeo][itheta][iphi][ipar]->SetMarkerColor(gr[igeo][itheta][iphi][ipar]->GetLineColor());
		gr[igeo][itheta][iphi][ipar]->SetMarkerStyle(SinglePointMarkerSytle);
		gr[igeo][itheta][iphi][ipar]->SetMarkerSize(SinglePointMarkerSize);
	      }
	      gr[igeo][itheta][iphi][ipar]->Draw("PEL");
	    }
          }
          l_ref_mom[ipar]->Draw();
        }
        c1->cd(6);
        HistName  = TString(thetaTitle) + TString(", ") + TString(phiTitle);
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.10,0.88,HistName.Data());
        leg->Draw("same");
	
	c1->cd();
	PlotLogo(0.15,0.89,0.35,0.5);
        c1->Print(EPSName.Data());
      }
      
      if(DoTelescopeAnalysis) {
	for(int idut=0;idut<NDUTS[itheta][iphi];idut++) { //begin loop over DUTs
	  double R[2];
	  bool NoRange;
	  
	  float MyTitleOffSet = 2.0;
	  float leftMargin    = 0.25;
	  float bottomMargin  = 0.15;
	  float topMargin     = 0.10;
	  
	  float TitleStartPositionX = gStyle->GetTitleX();
	  float TitleStartPositionY = gStyle->GetTitleY();
	  
	  gStyle->SetTitleX(0.6);
	  gStyle->SetTitleY(1.0);
	  
	  //Reference histos
	  HistName  = TString("href_TelResolUAtDUT");
	  HistTitle = TString("Tel-ResolU @ ") + DUTNames[itheta][iphi][idut];
          TH1F href_TelResolUAtDUT(HistName.Data(),
				   HistTitle.Data(),
				   100,RX[0],RX[1]);
	  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
	  href_TelResolUAtDUT.SetXTitle(HistName.Data());
	  href_TelResolUAtDUT.GetXaxis()->CenterTitle(true);
	  href_TelResolUAtDUT.SetYTitle("Tel-ResolU (#mum)");
	  href_TelResolUAtDUT.GetYaxis()->CenterTitle(true);
	  href_TelResolUAtDUT.GetYaxis()->SetTitleOffset(MyTitleOffSet);
	  href_TelResolUAtDUT.SetLineColor(1);
	  href_TelResolUAtDUT.SetLineWidth(1);
	  href_TelResolUAtDUT.SetNdivisions(5 + 100*5,"X");
	  href_TelResolUAtDUT.SetNdivisions(5 + 100*5,"Y");
	  href_TelResolUAtDUT.GetXaxis()->SetTitleSize(TheSize);
	  href_TelResolUAtDUT.GetXaxis()->SetLabelSize(TheSize);
	  href_TelResolUAtDUT.GetYaxis()->SetTitleSize(TheSize);
	  href_TelResolUAtDUT.GetYaxis()->SetLabelSize(TheSize);
	  href_TelResolUAtDUT.SetStats(false);
	  
	  HistName  = TString("href_TelResolVAtDUT");
	  HistTitle = TString("Tel-ResolV @ ") + DUTNames[itheta][iphi][idut];
          TH1F href_TelResolVAtDUT(HistName.Data(),
				   HistTitle.Data(),
				   100,RX[0],RX[1]);
	  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
	  href_TelResolVAtDUT.SetXTitle(HistName.Data());
	  href_TelResolVAtDUT.GetXaxis()->CenterTitle(true);
	  href_TelResolVAtDUT.SetYTitle("Tel-ResolV (#mum)");
	  href_TelResolVAtDUT.GetYaxis()->CenterTitle(true);
	  href_TelResolVAtDUT.GetYaxis()->SetTitleOffset(MyTitleOffSet);
	  href_TelResolVAtDUT.SetLineColor(1);
	  href_TelResolVAtDUT.SetLineWidth(1);
	  href_TelResolVAtDUT.SetNdivisions(5 + 100*5,"X");
	  href_TelResolVAtDUT.SetNdivisions(5 + 100*5,"Y");
	  href_TelResolVAtDUT.GetXaxis()->SetTitleSize(TheSize);
	  href_TelResolVAtDUT.GetXaxis()->SetLabelSize(TheSize);
	  href_TelResolVAtDUT.GetYaxis()->SetTitleSize(TheSize);
	  href_TelResolVAtDUT.GetYaxis()->SetLabelSize(TheSize);
	  href_TelResolVAtDUT.SetStats(false);
	  
	  HistName  = TString("href_TelResolUVAreaAtDUT");
	  HistTitle = TString("Tel-Resol S @ ") + DUTNames[itheta][iphi][idut];
          TH1F href_TelResolUVAreaAtDUT(HistName.Data(),
				   HistTitle.Data(),
				   100,RX[0],RX[1]);
	  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
	  href_TelResolUVAreaAtDUT.SetXTitle(HistName.Data());
	  href_TelResolUVAreaAtDUT.GetXaxis()->CenterTitle(true);
	  href_TelResolUVAreaAtDUT.SetYTitle("Tel-Resol 1-#sigma S (#mum^{2})");
	  href_TelResolUVAreaAtDUT.GetYaxis()->CenterTitle(true);
	  href_TelResolUVAreaAtDUT.GetYaxis()->SetTitleOffset(MyTitleOffSet);
	  href_TelResolUVAreaAtDUT.SetLineColor(1);
	  href_TelResolUVAreaAtDUT.SetLineWidth(1);
	  href_TelResolUVAreaAtDUT.SetNdivisions(5 + 100*5,"X");
	  href_TelResolUVAreaAtDUT.SetNdivisions(5 + 100*5,"Y");
	  href_TelResolUVAreaAtDUT.GetXaxis()->SetTitleSize(TheSize);
	  href_TelResolUVAreaAtDUT.GetXaxis()->SetLabelSize(TheSize);
	  href_TelResolUVAreaAtDUT.GetYaxis()->SetTitleSize(TheSize);
	  href_TelResolUVAreaAtDUT.GetYaxis()->SetLabelSize(TheSize);
	  href_TelResolUVAreaAtDUT.SetStats(false);
	  
	  HistName  = TString("href_NbkgAtDUT");
	  HistTitle = TString("N_{bkg} @ ") + DUTNames[itheta][iphi][idut];
          TH1F href_NbkgAtDUT(HistName.Data(),
				   HistTitle.Data(),
				   100,RX[0],RX[1]);
	  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
	  href_NbkgAtDUT.SetXTitle(HistName.Data());
	  href_NbkgAtDUT.GetXaxis()->CenterTitle(true);
	  href_NbkgAtDUT.SetYTitle("N_{bkg} hits");
	  href_NbkgAtDUT.GetYaxis()->CenterTitle(true);
	  href_NbkgAtDUT.GetYaxis()->SetTitleOffset(MyTitleOffSet);
	  href_NbkgAtDUT.SetLineColor(1);
	  href_NbkgAtDUT.SetLineWidth(1);
	  href_NbkgAtDUT.SetNdivisions(5 + 100*5,"X");
	  href_NbkgAtDUT.SetNdivisions(5 + 100*5,"Y");
	  href_NbkgAtDUT.GetXaxis()->SetTitleSize(TheSize);
	  href_NbkgAtDUT.GetXaxis()->SetLabelSize(TheSize);
	  href_NbkgAtDUT.GetYaxis()->SetTitleSize(TheSize);
	  href_NbkgAtDUT.GetYaxis()->SetLabelSize(TheSize);
	  href_NbkgAtDUT.SetStats(false);
	  
	  c1->Clear();
          c1->Divide(3,2);
          c1->cd(1);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(leftMargin);
          gPad->SetBottomMargin(bottomMargin);
	  gPad->SetTopMargin(topMargin);
	  gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
	  R[0] = +1.0e+20;
	  R[1] = -1.0e+20;
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->GetN();i++) {
	      double x,y;
	      gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href_TelResolUAtDUT.SetMinimum(R[0]);
	  href_TelResolUAtDUT.SetMaximum(R[1]);
          href_TelResolUAtDUT.Draw();
          for(int igeo=0;igeo<NGeometries;igeo++) {
            if(gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->GetN() > 0) {
	      if(gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->GetN() == 1) {
		gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->GetLineColor());
		gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetMarkerStyle(SinglePointMarkerSytle);
		gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->SetMarkerSize(SinglePointMarkerSize);
	      }
	      gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->Draw("PEL");
	    }
          }
          c1->cd(2);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(leftMargin);
          gPad->SetBottomMargin(bottomMargin);
	  gPad->SetTopMargin(topMargin);
	  gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
	  R[0] = +1.0e+20;
	  R[1] = -1.0e+20;
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->GetN();i++) {
	      double x,y;
	      gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href_TelResolVAtDUT.SetMinimum(R[0]);
	  href_TelResolVAtDUT.SetMaximum(R[1]);
          href_TelResolVAtDUT.Draw();
          for(int igeo=0;igeo<NGeometries;igeo++) {
            if(gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->GetN() > 0) {
	      if(gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->GetN() == 1) {
		gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->GetLineColor());
		gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetMarkerStyle(SinglePointMarkerSytle);
		gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->SetMarkerSize(SinglePointMarkerSize);
	      }
	      gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->Draw("PEL");
	    }
          }
          c1->cd(4);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(leftMargin);
          gPad->SetBottomMargin(bottomMargin);
	  gPad->SetTopMargin(topMargin);
	  gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
	  R[0] = +1.0e+20;
	  R[1] = -1.0e+20;
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->GetN();i++) {
	      double x,y;
	      gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href_TelResolUVAreaAtDUT.SetMinimum(R[0]);
	  href_TelResolUVAreaAtDUT.SetMaximum(R[1]);
          href_TelResolUVAreaAtDUT.Draw();
          for(int igeo=0;igeo<NGeometries;igeo++) {
            if(gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->GetN() > 0) {
	      if(gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->GetN() == 1) {
		gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->GetLineColor());
		gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetMarkerStyle(SinglePointMarkerSytle);
		gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->SetMarkerSize(SinglePointMarkerSize);
	      }
	      gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->Draw("PEL");
	    }
          }
          c1->cd(5);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(leftMargin);
          gPad->SetBottomMargin(bottomMargin);
	  gPad->SetTopMargin(topMargin);
	  gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
	  R[0] = +1.0e+20;
	  R[1] = -1.0e+20;
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_NbkgAtDUT[igeo][itheta][iphi][idut]->GetN();i++) {
	      double x,y;
	      gr_NbkgAtDUT[igeo][itheta][iphi][idut]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href_NbkgAtDUT.SetMinimum(R[0]);
	  href_NbkgAtDUT.SetMaximum(R[1]);
          href_NbkgAtDUT.Draw();
          for(int igeo=0;igeo<NGeometries;igeo++) {
            if(gr_NbkgAtDUT[igeo][itheta][iphi][idut]->GetN() > 0) {
	      if(gr_NbkgAtDUT[igeo][itheta][iphi][idut]->GetN() == 1) {
		gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetMarkerColor(gr_NbkgAtDUT[igeo][itheta][iphi][idut]->GetLineColor());
		gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetMarkerStyle(SinglePointMarkerSytle);
		gr_NbkgAtDUT[igeo][itheta][iphi][idut]->SetMarkerSize(SinglePointMarkerSize);
	      }
	      gr_NbkgAtDUT[igeo][itheta][iphi][idut]->Draw("PEL");
	    }
          }
          c1->cd(6);
	  HistName  = TString(thetaTitle) + TString(", ") + TString(phiTitle);
          latex->SetTextColor(kBlack);
          latex->DrawLatex(0.10,0.88,HistName.Data());
          leg->Draw("same");
	  
	  c1->cd();
	  PlotLogo(0.20,0.89,0.80,0.78);
          c1->Print(EPSName.Data());
	  
	  gStyle->SetTitleX(TitleStartPositionX);
	  gStyle->SetTitleY(TitleStartPositionY);
	} // end loop over DUTs
      }
      
      if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
        c1->Clear();
        c1->Divide(2,2);
        c1->cd(1);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        double Rtmp[2];
        Rtmp[0] = +1.0e+20;
        Rtmp[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
          for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_Nhits_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_Nhits_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(Rtmp[0] > y) Rtmp[0] = y;
	      if(Rtmp[1] < y) Rtmp[1] = y;
	    }
	  }
          bool NoRange = false;
          if(Rtmp[0] == +1.0e+20 && Rtmp[1] == -1.0e+20) NoRange = true;
          if(NoRange) {
	    Rtmp[0] = 1.0e-6;
	    Rtmp[1] = 1.0;
          }
          else {
	    delta = Rtmp[1] - Rtmp[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = Rtmp[0];
	    porcent = 0.1;
	    min = Rtmp[0];
	    Rtmp[0] -= porcent*delta;
	    Rtmp[1] += porcent*delta;
	    if(Rtmp[0] < 0.0) Rtmp[0] = min*min_frac;
          }
          href_Nhits_vs_theta_phi_mom->SetMinimum(Rtmp[0]);
          href_Nhits_vs_theta_phi_mom->SetMaximum(Rtmp[1]);
        }
        href_Nhits_vs_theta_phi_mom->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_Nhits_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_Nhits_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_Nhits_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_Nhits_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_Nhits_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_Nhits_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_Nhits_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
        }
        l_NhitsMin_vs_mom->Draw("PEL");
        c1->cd(3);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        Rtmp[0] = +1.0e+20;
        Rtmp[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
          for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(Rtmp[0] > y) Rtmp[0] = y;
	      if(Rtmp[1] < y) Rtmp[1] = y;
	    }
	  }
          bool NoRange = false;
          if(Rtmp[0] == +1.0e+20 && Rtmp[1] == -1.0e+20) NoRange = true;
          if(NoRange) {
	    Rtmp[0] = 1.0e-6;
	    Rtmp[1] = 1.0;
          }
          else {
	    delta = Rtmp[1] - Rtmp[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = Rtmp[0];
	    porcent = 0.1;
	    min = Rtmp[0];
	    Rtmp[0] -= porcent*delta;
	    Rtmp[1] += porcent*delta;
	    if(Rtmp[0] < 0.0) Rtmp[0] = min*min_frac;
          }
          href_MatBudget_vs_theta_phi_mom->SetMinimum(Rtmp[0]);
          href_MatBudget_vs_theta_phi_mom->SetMaximum(Rtmp[1]);
        }
        href_MatBudget_vs_theta_phi_mom->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_MatBudget_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_MatBudget_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
        }
        c1->cd(2);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        Rtmp[0] = +1.0e+20;
        Rtmp[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
          for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(Rtmp[0] > y) Rtmp[0] = y;
	      if(Rtmp[1] < y) Rtmp[1] = y;
	    }
	  }
          bool NoRange = false;
          if(Rtmp[0] == +1.0e+20 && Rtmp[1] == -1.0e+20) NoRange = true;
          if(NoRange) {
	    Rtmp[0] = 1.0e-6;
	    Rtmp[1] = 1.0;
          }
          else {
	    delta = Rtmp[1] - Rtmp[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = Rtmp[0];
	    porcent = 0.1;
	    min = Rtmp[0];
	    Rtmp[0] -= porcent*delta;
	    Rtmp[1] += porcent*delta;
	    if(Rtmp[0] < 0.0) Rtmp[0] = min*min_frac;
          }
          href_1stPointDist_vs_theta_phi_mom->SetMinimum(Rtmp[0]);
          href_1stPointDist_vs_theta_phi_mom->SetMaximum(Rtmp[1]);
        }
        href_1stPointDist_vs_theta_phi_mom->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_MatBudget_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_1stPorintDisToRef_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
        }
        c1->cd(4);
	HistName  = TString(thetaTitle) + TString(", ") + TString(phiTitle);
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.10,0.88,HistName.Data());
        leg->Draw("same");
	
	c1->cd();
	PlotLogo(0.17,0.89,0.5,0.5);
        c1->Print(EPSName.Data());
      }
	
    }
    
    if(TrkResolAnalysisPars[0].PlotDOCAatHighMom && TrkResolAnalysisPars[0].PlotDOCAvsMonFit) {
      for(int igeo=0;igeo<NGeometries;igeo++) {
        c1->Clear();
	c1->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
        TF1 func("func","sqrt( pow([0],2) + pow([1]/(x*pow(TMath::Sin([2]),1.5)),2) )",href[par_idx]->GetXaxis()->GetXmin(),href[par_idx]->GetXaxis()->GetXmax());
        func.SetParameters(a_param[igeo][itheta],b_param[igeo][itheta],thetaArr[itheta]);
        func.SetLineColor(gr_aveVsPhi[igeo][itheta][par_idx]->GetLineColor());
        func.SetLineWidth(2);
        func.SetLineStyle(2);
      
        TLine l_a(href[par_idx]->GetXaxis()->GetXmin(),a_param[igeo][itheta]/global->GetDistanceUnit("um"),href[par_idx]->GetXaxis()->GetXmax(),a_param[igeo][itheta]/global->GetDistanceUnit("um"));
        l_a.SetLineColor(gr_aveVsPhi[igeo][itheta][par_idx]->GetLineColor());
        l_a.SetLineWidth(2);
        l_a.SetLineStyle(2);
      
        href[par_idx]->Draw();
        if(gr_aveVsPhi[igeo][itheta][par_idx]->GetN()) gr_aveVsPhi[igeo][itheta][par_idx]->Draw("PEL");
        func.Draw("same");
        l_a.Draw();
        c1->Print(EPSName.Data());
	c1->SetLogy(false);
      }
    }

    if(TrkResolAnalysisPars[0].PlotOnlyPhiAveraged || TrkResolAnalysisPars[0].PlotPhiAveraged) {
      c1->Clear();
      c1->Divide(3,2);
      int counter_pad = 0;
      TString Titles[Nmax_params];
      TString YTitles[Nmax_params];
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
        counter_pad++;
        if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++; 
        c1->cd(counter_pad);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
      
        Titles[ipar] = TString(href[ipar]->GetTitle());
        HistName = TString("Ave ") + TString(href[ipar]->GetTitle());
        href[ipar]->SetTitle(HistName.Data());
      
        YTitles[ipar] = TString(href[ipar]->GetYaxis()->GetTitle());
        HistName = TString("Ave ") + TString(href[ipar]->GetYaxis()->GetTitle());
        href[ipar]->SetYTitle(HistName.Data());
      
        double R[2];
        R[0] = +1.0e+20;
        R[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetN();i++) {
	      double x,y;
	      gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href[ipar]->SetMinimum(R[0]);
	  href[ipar]->SetMaximum(R[1]);
        }
      
        href[ipar]->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
	  if(gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetN() > 0) {
	    if(gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetN() == 1) {
	      gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetMarkerColor(gr_aveVsPhi_NoErr[igeo][itheta][ipar]->GetLineColor());
	      gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_aveVsPhi_NoErr[igeo][itheta][ipar]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_aveVsPhi_NoErr[igeo][itheta][ipar]->Draw("PEL");
	  }
        }
      }
      c1->cd(6);
      HistName  = TString(thetaTitle);
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      
      c1->cd();
      PlotLogo(0.15,0.89,0.35,0.5);
      c1->Print(EPSName.Data());
    
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
        href[ipar]->SetTitle(Titles[ipar].Data());
        href[ipar]->SetYTitle(YTitles[ipar].Data());
      }
    }
  
    TString MMMTitles;
    TString MMMYTitles;
    MMMTitles = TString(href_1stPointDist_vs_theta_phi_mom->GetTitle());
    HistName  = TString("Ave ") + TString(href_1stPointDist_vs_theta_phi_mom->GetTitle());
    href_1stPointDist_vs_theta_phi_mom->SetTitle(HistName.Data());
    
    MMMYTitles = TString(href_1stPointDist_vs_theta_phi_mom->GetYaxis()->GetTitle());
    HistName   = TString("Ave ") + TString(href_1stPointDist_vs_theta_phi_mom->GetYaxis()->GetTitle());
    href_1stPointDist_vs_theta_phi_mom->SetYTitle(HistName.Data());
    
    if(polarVariable == TString("theta"))         sprintf(ytitle,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit(ThetaUnits));
    else if(polarVariable == TString("costheta")) sprintf(ytitle,"cos(#theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(ytitle,"#eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    TString OriginalTitle(href_a_vs_phi->GetTitle());
    HistTitle = OriginalTitle + TString(" (") + TString(ytitle) + TString(")");
    href_a_vs_phi->SetTitle(HistTitle.Data());
    
    TString OriginalTitle_Nhits(href_Nhits_vs_phi->GetTitle());
    HistTitle = OriginalTitle_Nhits + TString(" (") + TString(ytitle) + TString(")");
    href_Nhits_vs_phi->SetTitle(HistTitle.Data());

    if((TrkResolAnalysisPars[0].PlotDOCAatHighMom || TrkResolAnalysisPars[0].PlotMaterialBudget) && 
       (TrkResolAnalysisPars[0].PlotOnlyPhiAveraged || TrkResolAnalysisPars[0].PlotPhiAveraged)) {
      c1->Clear();
      c1->Divide(2,2);
      if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
        c1->cd(2);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        double Rmmm[2];
        Rmmm[0] = +1.0e+20;
        Rmmm[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
          for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetN();i++) {
	      double x,y;
	      gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(Rmmm[0] > y) Rmmm[0] = y;
	      if(Rmmm[1] < y) Rmmm[1] = y;
	    }
          }
          bool NoRange = false;
          if(Rmmm[0] == +1.0e+20 && Rmmm[1] == -1.0e+20) NoRange = true;
          if(NoRange) {
	    Rmmm[0] = 1.0e-6;
	    Rmmm[1] = 1.0;
          }
          else {
	    delta = Rmmm[1] - Rmmm[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = Rmmm[0];
	    porcent = 0.1;
	    min = Rmmm[0];
	    Rmmm[0] -= porcent*delta;
	    Rmmm[1] += porcent*delta;
	    if(Rmmm[0] < 0.0) Rmmm[0] = min*min_frac;
          }
          href_1stPointDist_vs_theta_phi_mom->SetMinimum(Rmmm[0]);
          href_1stPointDist_vs_theta_phi_mom->SetMaximum(Rmmm[1]);
        }
        href_1stPointDist_vs_theta_phi_mom->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetN() > 0) {
	    if(gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetN() == 1) {
	      gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetMarkerColor(gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->GetLineColor());
	      gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_aveVsPhi1stPorintDisToRef_vs_mom[igeo][itheta]->Draw("PEL");
	  }
        }      
        c1->cd(1);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        href_Nhits_vs_phi->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_Nhits_vs_Phi[igeo][itheta]->GetN() > 0) gr_Nhits_vs_Phi[igeo][itheta]->Draw("PEL");
      
          double Nhits_ave = 0;
          for(int i=0;i<gr_Nhits_vs_Phi[igeo][itheta]->GetN();i++) {
	    double x,y;
	    gr_Nhits_vs_Phi[igeo][itheta]->GetPoint(i,x,y);
	    Nhits_ave += y;
          }
          if(gr_Nhits_vs_Phi[igeo][itheta]->GetN() > 0) Nhits_ave /= gr_Nhits_vs_Phi[igeo][itheta]->GetN();
      
          TLine* l_hits_ave = new TLine(href_Nhits_vs_phi->GetXaxis()->GetXmin(),Nhits_ave,
				        href_Nhits_vs_phi->GetXaxis()->GetXmax(),Nhits_ave);
          l_hits_ave->SetLineColor(gr_Nhits_vs_Phi[igeo][itheta]->GetLineColor());
          l_hits_ave->SetLineWidth(2);
          l_hits_ave->SetLineStyle(2);
          l_hits_ave->Draw();
        }
      }
      if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
        c1->cd(3);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        double R[2];
        R[0] = +1.0e+20;
        R[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
          for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_a_vs_Phi[igeo][itheta]->GetN();i++) {
	      double x,y;
	      gr_a_vs_Phi[igeo][itheta]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	
	    if(R[0] > a_param[igeo][itheta]/global->GetDistanceUnit("um")) R[0] = a_param[igeo][itheta]/global->GetDistanceUnit("um");
	    if(R[1] < a_param[igeo][itheta]/global->GetDistanceUnit("um")) R[1] = a_param[igeo][itheta]/global->GetDistanceUnit("um");
          }
          bool NoRange = false;
          if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
          if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
          }
          else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
          }
          href_a_vs_phi->SetMinimum(R[0]);
          href_a_vs_phi->SetMaximum(R[1]);
        }
    
        href_a_vs_phi->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_a_vs_Phi[igeo][itheta]->GetN() > 0) {
	    if(gr_a_vs_Phi[igeo][itheta]->GetN() == 1) {
	      gr_a_vs_Phi[igeo][itheta]->SetMarkerColor(gr_a_vs_Phi[igeo][itheta]->GetLineColor());
	      gr_a_vs_Phi[igeo][itheta]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_a_vs_Phi[igeo][itheta]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_a_vs_Phi[igeo][itheta]->Draw("PEL");
	  }
      
          l_a_ave_for_theta[igeo][itheta] = new TLine(href_a_vs_phi->GetXaxis()->GetXmin(),a_param[igeo][itheta]/global->GetDistanceUnit("um"),
	                                              href_a_vs_phi->GetXaxis()->GetXmax(),a_param[igeo][itheta]/global->GetDistanceUnit("um"));
          l_a_ave_for_theta[igeo][itheta]->SetLineColor(gr_a_vs_Phi[igeo][itheta]->GetLineColor());
          l_a_ave_for_theta[igeo][itheta]->SetLineWidth(2);
          l_a_ave_for_theta[igeo][itheta]->SetLineStyle(2);
          l_a_ave_for_theta[igeo][itheta]->Draw();
        }
      }
      c1->cd(4);      
      HistName  = TString(thetaTitle);
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      
      c1->cd();
      PlotLogo(0.17,0.89,0.5,0.5);
      c1->Print(EPSName.Data());
    
      href_1stPointDist_vs_theta_phi_mom->SetTitle(MMMTitles.Data());
      href_1stPointDist_vs_theta_phi_mom->SetYTitle(MMMYTitles.Data());
      href_a_vs_phi->SetTitle(OriginalTitle.Data());
      href_Nhits_vs_phi->SetTitle(OriginalTitle_Nhits.Data());
    }
    
  }

  if(TrkResolAnalysisPars[0].PlotPerformancesVsTheta && (TrkResolAnalysisPars[0].PlotDOCAatHighMom || TrkResolAnalysisPars[0].PlotMaterialBudget)) {
    c1->Clear();
    c1->Divide(2,2);
    if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      //gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
      href_a->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) gr_a[igeo]->Draw("PEL");
      c1->cd(3);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      //gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
      href_b->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) gr_b[igeo]->Draw("PEL");
    }
    if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
      c1->cd(4);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      double Rnnn[2];
      Rnnn[0] = +1.0e+20;
      Rnnn[1] = -1.0e+20;
      if(!TrkResolAnalysisPars[0].SameRange) {
        for(int igeo=0;igeo<NGeometries;igeo++) {
          for(int i=0;i<gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->GetN();i++) {
	    double x,y;
	    gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->GetPoint(i,x,y);
	    if(y < 0.0) continue;
	    if(Rnnn[0] > y) Rnnn[0] = y;
	    if(Rnnn[1] < y) Rnnn[1] = y;
          }
        }
        bool NoRange = false;
        if(Rnnn[0] == +1.0e+20 && Rnnn[1] == -1.0e+20) NoRange = true;
        if(NoRange) {
          Rnnn[0] = 1.0e-6;
          Rnnn[1] = 1.0;
        }
        else {
          delta = Rnnn[1] - Rnnn[0];
          if(TMath::Abs(delta) < 1.0e-8) delta = Rnnn[0];
          porcent = 0.1;
	  min = Rnnn[0];
          Rnnn[0] -= porcent*delta;
          Rnnn[1] += porcent*delta;
          if(Rnnn[0] < 0.0) Rnnn[0] = min*min_frac;
        }
        href_Ave1stPointDist_vs_theta->SetMinimum(Rnnn[0]);
        href_Ave1stPointDist_vs_theta->SetMaximum(Rnnn[1]);
      }
      href_Ave1stPointDist_vs_theta->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) gr_aveVsPhiAndMom1stPorintDisToRef[igeo]->Draw("PEL");
    }
    c1->cd(2);
    leg->Draw("same");
    
    c1->cd();
    PlotLogo(0.17,0.89,0.5,0.5);
    c1->Print(EPSName.Data());
  }

  if(TrkResolAnalysisPars[0].PlotPerformancesVsTheta) {
    for(int iP=0; iP<MyNbins_mom_redu; iP++) {
      double mom;
      if(TrkResolAnalysisPars[0].UseAllMomVals) mom = momArr[iP];
      else                                      mom = momArr[0] + (iP+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
      mom /= global->GetUnit(MonUnits);
    
      c1->Clear();
      c1->Divide(3,2);
      int counter_pad = 0;
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
        counter_pad++;
        if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++;
        c1->cd(counter_pad);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	//gPad->SetLogy(TrkResolAnalysisPars[0].UseLogYAxes);
        double R[2];
        R[0] = +1.0e+20;
        R[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_aveSigmaParVsTheta[igeo][iP][ipar]->GetN();i++) {
	      double x,y;
	      gr_aveSigmaParVsTheta[igeo][iP][ipar]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href_vs_theta[ipar]->SetMinimum(R[0]);
	  href_vs_theta[ipar]->SetMaximum(R[1]);
        }
      
        href_vs_theta[ipar]->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
	  if(gr_aveSigmaParVsTheta[igeo][iP][ipar]->GetN() > 0) gr_aveSigmaParVsTheta[igeo][iP][ipar]->Draw("PEL");
        }
      }
      c1->cd(6);
      sprintf(ytitle,"%.3f",mom);
      HistName  = TString("p = ") + TString(ytitle) + TString(" ") + MonUnits;
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      
      c1->cd();
      PlotLogo(0.15,0.89,0.35,0.5);
      c1->Print(EPSName.Data());
    
      if(TrkResolAnalysisPars[0].PlotMaterialBudget) {
        c1->Clear();
        c1->Clear();
        c1->Divide(2,2);
        c1->cd(1);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
        double R[2];
        R[0] = +1.0e+20;
        R[1] = -1.0e+20;
        if(!TrkResolAnalysisPars[0].SameRange) {
          for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_aveNhitsVsTheta[igeo][iP]->GetN();i++) {
	      double x,y;
	      gr_aveNhitsVsTheta[igeo][iP]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
          }
          bool NoRange = false;
          if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
          if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
          }
          else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
          }
          href_AveNhits_vs_theta->SetMinimum(R[0]);
          href_AveNhits_vs_theta->SetMaximum(R[1]);
        }
        href_AveNhits_vs_theta->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_aveNhitsVsTheta[igeo][iP]->GetN() > 0) gr_aveNhitsVsTheta[igeo][iP]->Draw("PEL");
        }
        l_NhitsMin_vs_theta->Draw("PEL");
        c1->cd(2);
        sprintf(ytitle,"%.3f",mom);
        HistName  = TString("p = ") + TString(ytitle) + TString(" ") + MonUnits;
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.10,0.88,HistName.Data());
        leg->Draw("same");
	
	c1->cd();
	PlotLogo(0.17,0.89,0.5,0.5);
        c1->Print(EPSName.Data());
      }
    }
    
  }

  c1->Print(EPSNameC.Data());

  if(SavePlots) {
    cout << endl;
    TString ROOTName = TheOutputFile + TString(".root");
    cout << "Saving track parameters resolution vs momentum to " << ROOTName.Data() << " file" << endl;
    TFile file(ROOTName.Data(),"UPDATE");

    TCanvas* c_trk_resol_vs_mom[NthetaPoints][NphiPoints];
    TCanvas* c_trk_resol_vs_mom_PhiAve[NthetaPoints];
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
      char theta[300];      
      if(polarVariable == TString("theta"))         sprintf(theta,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit(ThetaUnits));
      else if(polarVariable == TString("costheta")) sprintf(theta,"cos(#theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
      else if(polarVariable == TString("eta"))      sprintf(theta,"#eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
      
      if(!TrkResolAnalysisPars[0].PlotOnlyPhiAveraged) {
        for(int iphi=0;iphi<NphiPoints;iphi++) {
          char phi[300];
          sprintf(phi,  "%.1f %s",phiArr[iphi]/global->GetUnit(ThetaUnits),ThetaUnits.Data());
	  
          HistName  = TString("c_trk_resol_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1);
          HistTitle = TString("Track parameters resolution vs monmentum for ") + TString(theta) + TString(",") + TString(phi);
          c_trk_resol_vs_mom[itheta][iphi] = new TCanvas(HistName.Data(),HistTitle.Data());
    
          c_trk_resol_vs_mom[itheta][iphi]->Clear();
          c_trk_resol_vs_mom[itheta][iphi]->Divide(3,2);
          int counter_pad = 0;
          for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
            counter_pad++;
            if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++; 
            c_trk_resol_vs_mom[itheta][iphi]->cd(counter_pad);
            gPad->SetFillColor(10);
            gPad->SetFrameFillColor(10);
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLeftMargin(0.20);
            gPad->SetBottomMargin(0.15);
            href[ipar]->Draw();
            for(int igeo=0;igeo<NGeometries;igeo++) gr[igeo][itheta][iphi][ipar]->Draw("PEL");
          }
          c_trk_resol_vs_mom[itheta][iphi]->cd(6);
	  HistName  = TString(theta) + TString(" ,") + TString(phi);
          latex->SetTextColor(kBlack);
          latex->DrawLatex(0.10,0.88,HistName.Data());
          leg->Draw("same");
	  
	  c_trk_resol_vs_mom[itheta][iphi]->cd();
	  PlotLogo(0.15,0.89,0.35,0.5);
          c_trk_resol_vs_mom[itheta][iphi]->Write();
        }
      }
      
      if(TrkResolAnalysisPars[0].PlotOnlyPhiAveraged || TrkResolAnalysisPars[0].PlotPhiAveraged) {
        HistName  = TString("c_trk_resol_PhiAve_vs_mom_theta") + long(itheta+1);
        HistTitle = TString("#phi-ave Track parameters resolution vs monmentum for ") + TString(theta);
        c_trk_resol_vs_mom_PhiAve[itheta] = new TCanvas(HistName.Data(),HistTitle.Data());
    
        c_trk_resol_vs_mom_PhiAve[itheta]->Clear();
        c_trk_resol_vs_mom_PhiAve[itheta]->Divide(3,2);
        int counter_pad = 0;
        TString Titles[Nmax_params];
        TString YTitles[Nmax_params];
        for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
          counter_pad++;
          if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++; 
          c_trk_resol_vs_mom_PhiAve[itheta]->cd(counter_pad);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(0.20);
          gPad->SetBottomMargin(0.15);
      
          Titles[ipar] = TString(href[ipar]->GetTitle());
          HistName = TString("Ave ") + TString(href[ipar]->GetTitle());
          href[ipar]->SetTitle(HistName.Data());
      
          YTitles[ipar] = TString(href[ipar]->GetYaxis()->GetTitle());
          HistName = TString("Ave ") + TString(href[ipar]->GetYaxis()->GetTitle());
          href[ipar]->SetYTitle(HistName.Data());
      
          href[ipar]->Draw();
          //for(int igeo=0;igeo<NGeometries;igeo++) gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->Draw("PEL");
	  //for(int igeo=0;igeo<NGeometries;igeo++) gr_aveVsPhi[igeo][itheta][ipar]->Draw("PEL");
	  for(int igeo=0;igeo<NGeometries;igeo++) gr_aveVsPhi_NoErr[igeo][itheta][ipar]->Draw("PEL");
        }
        c_trk_resol_vs_mom_PhiAve[itheta]->cd(6);
	HistName  = TString(theta);
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.10,0.88,HistName.Data());
        leg->Draw("same");
	
	c_trk_resol_vs_mom_PhiAve[itheta]->cd();
	PlotLogo(0.15,0.89,0.35,0.5);
        c_trk_resol_vs_mom_PhiAve[itheta]->Write();
	
        for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
          href[ipar]->SetTitle(Titles[ipar].Data());
          href[ipar]->SetYTitle(YTitles[ipar].Data());
        }
      }
      
    }
    
    //Saving the TGraph objects
    for(int igeo=0;igeo<NGeometries;igeo++) {
      if(TrkResolAnalysisPars[0].PlotDOCAatHighMom) {
        gr_a[igeo]->Write();
	gr_b[igeo]->Write();
      }
      
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
	for(int iphi=0;iphi<NphiPoints;iphi++) {
	  if(!TrkResolAnalysisPars[0].PlotOnlyPhiAveraged) {
	    for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) gr[igeo][itheta][iphi][ipar]->Write();
	  }
	  
	  if(DoTelescopeAnalysis) {
	    for(int idut=0;idut<NDUTS[itheta][iphi];idut++) { //begin loop over DUTs
	      gr_TelResolUAtDUT[igeo][itheta][iphi][idut]->Write();
	      gr_TelResolVAtDUT[igeo][itheta][iphi][idut]->Write();
	      gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut]->Write();
	      gr_NbkgAtDUT[igeo][itheta][iphi][idut]->Write();
	    }
	  }
	  
        }
        
        if(TrkResolAnalysisPars[0].PlotOnlyPhiAveraged || TrkResolAnalysisPars[0].PlotPhiAveraged) {
          for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
	    gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar]->Write();
	    gr_aveVsPhi[igeo][itheta][ipar]->Write();
	    gr_aveVsPhi_NoErr[igeo][itheta][ipar]->Write();
	  }
	}
	
      }
      
      if(TrkResolAnalysisPars[0].PlotPerformancesVsTheta) {
        for(int ip=0;ip<MyNbins_mom_redu;ip++) {
	  if(TrkResolAnalysisPars[0].PlotMaterialBudget)  gr_aveNhitsVsTheta[igeo][ip]->Write();
	
	  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) gr_aveSigmaParVsTheta[igeo][ip][ipar]->Write();
	}
      }
      
    }
  
    file.Close();
  }

  //Freeing the memory
  for(int igeo=0;igeo<NGeometries;igeo++) {
    delete gr_a[igeo];
    delete gr_b[igeo];
      
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      for(int iphi=0;iphi<NphiPoints;iphi++) {
        for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) delete gr[igeo][itheta][iphi][ipar];
	
	if(DoTelescopeAnalysis) {
	  for(int idut=0;idut<MaxDUTs;idut++) {
	    delete  gr_TelResolUAtDUT[igeo][itheta][iphi][idut];
	    delete  gr_TelResolVAtDUT[igeo][itheta][iphi][idut];
	    delete  gr_TelResolUVAreaAtDUT[igeo][itheta][iphi][idut];
	    delete  gr_NbkgAtDUT[igeo][itheta][iphi][idut];
	  }
	}
	
      }
      
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
	delete gr_aveVsPhi_ErrorRMS[igeo][itheta][ipar];
	delete gr_aveVsPhi[igeo][itheta][ipar];
      }
      
    }
  }

  if(DoTelescopeAnalysis) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      delete [] NDUTS[itheta];
      for(int iphi=0;iphi<NphiPoints;iphi++) {
        delete [] DUTNames[itheta][iphi];
      }
      delete [] DUTNames[itheta];
    }
    delete [] NDUTS;
    delete [] DUTNames;
  
    for(int igeo=0;igeo<NGeometries;igeo++) {
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
        for(int iphi=0;iphi<NphiPoints;iphi++) {
	  for(int imom=0;imom<NmomPoints;imom++) {
	    delete [] DUTLayerNames[igeo][itheta][iphi][imom];
            delete [] TelResolUAtDUT[igeo][itheta][iphi][imom];
            delete [] TelResolVAtDUT[igeo][itheta][iphi][imom];
            delete [] TelResolUVAreaAtDUT[igeo][itheta][iphi][imom];
	    delete [] NbkgAtDUT[igeo][itheta][iphi][imom];
            delete [] MomValue[igeo][itheta][iphi][imom];
	  }
	  delete [] DUTLayerNames[igeo][itheta][iphi];
          delete [] TelResolUAtDUT[igeo][itheta][iphi];
          delete [] TelResolVAtDUT[igeo][itheta][iphi];
          delete [] TelResolUVAreaAtDUT[igeo][itheta][iphi];
	  delete [] NbkgAtDUT[igeo][itheta][iphi];
          delete [] MomValue[igeo][itheta][iphi];
	  delete [] gr_TelResolUAtDUT[igeo][itheta][iphi];
          delete [] gr_TelResolVAtDUT[igeo][itheta][iphi];
          delete [] gr_TelResolUVAreaAtDUT[igeo][itheta][iphi];
          delete [] gr_NbkgAtDUT[igeo][itheta][iphi];
        }
        delete [] DUTLayerNames[igeo][itheta];
        delete [] TelResolUAtDUT[igeo][itheta];
        delete [] TelResolVAtDUT[igeo][itheta];
        delete [] TelResolUVAreaAtDUT[igeo][itheta];
        delete [] NbkgAtDUT[igeo][itheta];
        delete [] MomValue[igeo][itheta];
	delete [] gr_TelResolUAtDUT[igeo][itheta];
        delete [] gr_TelResolVAtDUT[igeo][itheta];
        delete [] gr_TelResolUVAreaAtDUT[igeo][itheta];
        delete [] gr_NbkgAtDUT[igeo][itheta];
      }
      delete [] DUTLayerNames[igeo];
      delete [] TelResolUAtDUT[igeo];
      delete [] TelResolVAtDUT[igeo];
      delete [] TelResolUVAreaAtDUT[igeo];
      delete [] NbkgAtDUT[igeo];
      delete [] MomValue[igeo];
      delete [] gr_TelResolUAtDUT[igeo];
      delete [] gr_TelResolVAtDUT[igeo];
      delete [] gr_TelResolUVAreaAtDUT[igeo];
      delete [] gr_NbkgAtDUT[igeo];
    }
    delete [] DUTLayerNames;
    delete [] TelResolUAtDUT;
    delete [] TelResolVAtDUT;
    delete [] TelResolUVAreaAtDUT;
    delete [] NbkgAtDUT;
    delete [] MomValue;
    delete [] gr_TelResolUAtDUT;
    delete [] gr_TelResolVAtDUT;
    delete [] gr_TelResolUVAreaAtDUT;
    delete [] gr_NbkgAtDUT;
  }
  
  return;
  
}
//====================================================================
void Guariguanchi::doTrackingPseudoEfficVsMomentum(void)
{

  //This function performs an analytical tracking pseudo-efficiency calculation
  //as a function of the particle momentum for specified values of polar angles
  
  CheckGeometriesComparability();
  CheckTelescopeConfigurations();
  
  cout << endl;
  cout << "Start doTrackingPseudoEffic  ";
  global->fWatch.Print();
  global->fWatch.Continue();

  double MinEffic = 1.0e-6;
  
  const int NGeometries(GeometryList.size());
  const int NthetaPoints(thetaArr.size());
  const int NphiPoints(phiArr.size());
  const int NmomPoints(momArr.size());
  
  int Somebins = Nbins_mom_redu;
  int Thebins  = NmomPoints;
  if(NmomPoints > Somebins) Thebins = Somebins;
  if(EfficAnalysisPars[0].UseAllMomVals) Thebins  = momArr.size();
  const int MyNbins_mom_redu(Thebins);

  GlobalFileCounter++;
  TString HistName,HistTitle;
  char ytitle[300];
  
  TLatex* latex = new TLatex();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.08);
  latex->SetTextColor(1);
  
  const int Nmax_params(8);
  std::vector<TString> title;
  std::vector<TString> titleY;
  std::vector<TString> ParamUnits;
  title.clear();
  titleY.clear();
  ParamUnits.clear();
  
  TString MonUnits("GeV/c");
  if(momVariable == TString("E"))             MonUnits = TString("GeV");
  else if(momVariable == TString("Ekin"))     MonUnits = TString("GeV");
  else if(momVariable == TString("EkinPerU")) MonUnits = TString("GeV");
  TString ThetaUnits("deg");
  
  if(NthetaPoints == 1) EfficAnalysisPars[0].PlotPerformancesVsTheta = false;
  if(NphiPoints   == 1) {
    EfficAnalysisPars[0].PlotPhiAveraged     = false;
    EfficAnalysisPars[0].PlotOnlyPhiAveraged = false;
  }

  TString  MyPolarVariable;
  if(polarVariable == TString("theta"))         MyPolarVariable = TString("#theta");
  else if(polarVariable == TString("costheta")) MyPolarVariable = TString("cos(#theta)");
  else if(polarVariable == TString("eta"))      MyPolarVariable = TString("#eta");
  
  TString MyMomentumVariable("p");
  if(momVariable == TString("E"))             MyMomentumVariable = TString("Energy");
  else if(momVariable == TString("Ekin"))     MyMomentumVariable = TString("E_{kin}");
  else if(momVariable == TString("EkinPerU")) MyMomentumVariable = TString("E_{kin}/u");
  
  //////////////////////////////
  ///// Defining the plots /////
  //////////////////////////////

  TGraph* gr_res_ave_vs_mom[NGeometries][NthetaPoints][NphiPoints][Nmax_params];
  TGraph* gr_effic_tot_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  TGraph* gr_effic_NoFakes_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  TGraph* gr_effic_1Fake_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  TGraph* gr_effic_2orMoreFakes_vs_mom[NGeometries][NthetaPoints][NphiPoints];
  
  TGraph* gr_res_ave_aveVsPhi_vs_mom[NGeometries][NthetaPoints][Nmax_params];
  TGraph* gr_effic_tot_aveVsPhi_vs_mom[NGeometries][NthetaPoints];
  TGraph* gr_effic_NoFakes_aveVsPhi_vs_mom[NGeometries][NthetaPoints];
  TGraph* gr_effic_1Fake_aveVsPhi_vs_mom[NGeometries][NthetaPoints];
  TGraph* gr_effic_2orMoreFakes_aveVsPhi_vs_mom[NGeometries][NthetaPoints];
  
  TGraphErrors* gr_PhiAveSigmaParVsTheta[NGeometries][MyNbins_mom_redu][Nmax_params];
  TGraphErrors* gr_PhiAveEffic_tot_VsTheta[NGeometries][MyNbins_mom_redu];
  TGraphErrors* gr_PhiAveEffic_NoFakes_VsTheta[NGeometries][MyNbins_mom_redu];
  TGraphErrors* gr_PhiAveEffic_1Fake_VsTheta[NGeometries][MyNbins_mom_redu];
  TGraphErrors* gr_PhiAveEffic_2orMoreFakes_VsTheta[NGeometries][MyNbins_mom_redu];

  int counter = 0;
  for(int igeo=0;igeo<NGeometries;igeo++) { //begin loop over geometries
    TrackerList[igeo]->GetTrajectory()->FillParametersNames(MonResolRepresentation);
    if(igeo == 0) {
      for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	HistName = TrackerList[igeo]->GetTrajectory()->GetParameterErrorName(ipar);
	HistName += TString(" vs ") + MyMomentumVariable;
	title.push_back(HistName);
	
	HistName =  TrackerList[igeo]->GetTrajectory()->GetParameterErrorName(ipar) + TString(" (");
	HistName += TrackerList[igeo]->GetTrajectory()->GetParameterErrorUnitTitle(ipar) + TString(")");
	titleY.push_back(HistName);
	
	ParamUnits.push_back(TrackerList[igeo]->GetTrajectory()->GetParameterErrorUnit(ipar));
      }
    }
    
    int mycolor = (counter+1==10)?49:(counter+1);
    if(mycolor == 5 || mycolor == 7) {
      mycolor++;
      counter++;
    }
    if(mycolor == 3) mycolor = kGreen+2;

    counter++;
    
    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //being loop over parameters
      for(int itheta=0;itheta<NthetaPoints;itheta++) { // Begin loop over theta values
	char thetaTitle[100];
	if(polarVariable == TString("theta"))         sprintf(thetaTitle,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
	else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(#theta) = %.1f",global->FromThetaToCosTheta(thetaArr[itheta]));
	else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"#eta = %.1f",global->FromThetaToEta(thetaArr[itheta]));
	  
	for(int iphi=0;iphi<NphiPoints;iphi++) { // begin loop over phi values
	  char phiTitle[100];
	  sprintf(phiTitle,"#phi = %.1f deg",phiArr[iphi]/global->GetAngleUnit("deg"));
	  
          gr_res_ave_vs_mom[igeo][itheta][iphi][ipar] = new TGraph();
          HistName   = TString("gr_res_ave_par") + long(ipar) + TString("_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_vs_mom_geo_") + long(igeo+1);
          HistTitle  = TString("Effic-average sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	  gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetName(HistName.Data());
          gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetTitle(HistTitle.Data());
          gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetMarkerColor(mycolor);
          gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetLineColor(mycolor);
          gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetLineWidth(2);
	  
	  if(ipar == 0) {
	    gr_effic_tot_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_effic_tot_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
            HistTitle = TString("Effic_{tot} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	    
	    gr_effic_NoFakes_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_effic_NoFakes_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
            HistTitle = TString("Effic_{NoFakes} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	    
	    gr_effic_1Fake_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_effic_1Fake_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
            HistTitle = TString("Effic_{1-Fake} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	    
	    gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi] = new TGraph();
            HistName = TString("gr_effic_2orMoreFakes_vs_mom_theta") + long(itheta+1) + TString("_phi") + long(iphi+1) + TString("_geo") + long(igeo+1);
            HistTitle = TString("Effic_{#geq 2 Fakes} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(", ") + TString(phiTitle) + TString(", and geo ") + GeometryList[igeo]->GetName();
	    gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetName(HistName.Data());
            gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetTitle(HistTitle.Data());
            gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetMarkerColor(mycolor);
            gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetLineColor(mycolor);
            gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetLineWidth(2);
	  }
	  
	  if(iphi == 0) {
	    gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar] = new TGraphErrors();
            HistName = TString("gr_res_ave_aveVsPhi_res_par") + long(ipar) + TString("_vs_mom_theta") + long(itheta+1) + TString("_vs_mom_geo_") + long(igeo+1);
            HistTitle = TString("Effic & #phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + TString(" and geo ") + GeometryList[igeo]->GetName();
            gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetName(HistName.Data());
            gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetTitle(HistTitle.Data());
            gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetMarkerColor(mycolor);
            gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetLineColor(mycolor);
            gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetLineWidth(2);
	    
	    
	    if(ipar == 0) {
	      gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta] = new TGraph();
              HistName = TString("gr_effic_tot_aveVsPhi_vs_mom_theta") + long(itheta+1) + TString("_geo") + long(igeo+1);
              HistTitle = TString("#phi average Effic_{tot} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + 
	                  TString(", and geo ") + GeometryList[igeo]->GetName();
	      gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetName(HistName.Data());
              gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetTitle(HistTitle.Data());
              gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(mycolor);
              gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetLineColor(mycolor);
              gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetLineWidth(2);
	    
	      gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta] = new TGraph();
              HistName = TString("gr_effic_NoFakes_aveVsPhi_vs_mom_theta") + long(itheta+1) + TString("_geo") + long(igeo+1);
              HistTitle = TString("#phi average Effic_{NoFakes} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + 
	                  TString(", and geo ") + GeometryList[igeo]->GetName();
	      gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetName(HistName.Data());
              gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetTitle(HistTitle.Data());
              gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(mycolor);
              gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetLineColor(mycolor);
              gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetLineWidth(2);
	    
	      gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta] = new TGraph();
              HistName = TString("gr_effic_1Fake_aveVsPhi_vs_mom_theta") + long(itheta+1) + TString("_geo") + long(igeo+1);
              HistTitle = TString("#phi average Effic_{1-Fake} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + 
	                  TString(", and geo ") + GeometryList[igeo]->GetName();
	      gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetName(HistName.Data());
              gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetTitle(HistTitle.Data());
              gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(mycolor);
              gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetLineColor(mycolor);
              gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetLineWidth(2);
	    
	      gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta] = new TGraph();
              HistName = TString("gr_effic_2orMoreFakes_aveVsPhi_vs_mom_theta") + long(itheta+1) + TString("_geo") + long(igeo+1);
              HistTitle = TString("#phi average Effic_{#geq 2 Fakes} vs ") + MyMomentumVariable + TString(" for ") + TString(thetaTitle) + 
	                  TString(", and geo ") + GeometryList[igeo]->GetName();
	      gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetName(HistName.Data());
              gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetTitle(HistTitle.Data());
              gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(mycolor);
              gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetLineColor(mycolor);
              gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetLineWidth(2);
	    }
	  }
	  
	} // end loop over phi values
      } // end loop over theta values
      
      for(int ip=0;ip<MyNbins_mom_redu;ip++) { // begin loop over momentum values
	double pi;
	if(EfficAnalysisPars[0].UseAllMomVals) pi = momArr[ip];
	else                                   pi = momArr[0] + (ip+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
        pi /= global->GetUnit(MonUnits);
	
	sprintf(ytitle,"%.3f",pi);

	gr_PhiAveSigmaParVsTheta[igeo][ip][ipar] = new TGraphErrors();
	HistName = TString("gr_PhiAveSigmaParVsTheta_res_par") + long(ipar) + TString("_mom") + long(ip+1) + TString("_geo_") + long(igeo+1);
        HistTitle = TString("Effic & #phi average sigma-par") + long(ipar+1) + TString(" vs ") + MyPolarVariable +  
                    TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ") + GeometryList[igeo]->GetName();
	gr_PhiAveSigmaParVsTheta[igeo][ip][ipar]->SetName(HistName.Data());
	gr_PhiAveSigmaParVsTheta[igeo][ip][ipar]->SetTitle(HistTitle.Data());
	gr_PhiAveSigmaParVsTheta[igeo][ip][ipar]->SetMarkerColor(mycolor);
	gr_PhiAveSigmaParVsTheta[igeo][ip][ipar]->SetLineColor(mycolor);
	gr_PhiAveSigmaParVsTheta[igeo][ip][ipar]->SetLineWidth(2);
	
	if(ipar == 0) {
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip] = new TGraphErrors();
	  HistName = TString("gr_PhiAveEffic_tot_VsTheta_mom") + TString("_mom") + long(ip+1) + TString("_geo_") + long(igeo+1);
          HistTitle = TString("#phi average Effic_{tot} vs ") + MyPolarVariable +  
                      TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ") + GeometryList[igeo]->GetName();
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip]->SetName(HistName.Data());
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip]->SetTitle(HistTitle.Data());
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip]->SetMarkerColor(mycolor);
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip]->SetLineColor(mycolor);
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip]->SetLineWidth(2);
	  
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip] = new TGraphErrors();
	  HistName = TString("gr_PhiAveEffic_NoFakes_VsTheta_mom") + TString("_mom") + long(ip+1) + TString("_geo_") + long(igeo+1);
          HistTitle = TString("#phi average Effic_{NoFakes} vs ") + MyPolarVariable +  
                      TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ") + GeometryList[igeo]->GetName();
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip]->SetName(HistName.Data());
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip]->SetTitle(HistTitle.Data());
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip]->SetMarkerColor(mycolor);
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip]->SetLineColor(mycolor);
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip]->SetLineWidth(2);
	  
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip] = new TGraphErrors();
	  HistName = TString("gr_PhiAveEffic_1Fake_VsTheta_mom") + TString("_mom") + long(ip+1) + TString("_geo_") + long(igeo+1);
          HistTitle = TString("#phi average Effic_{1-Fake} vs ") + MyPolarVariable +  
                      TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ") + GeometryList[igeo]->GetName();
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip]->SetName(HistName.Data());
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip]->SetTitle(HistTitle.Data());
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip]->SetMarkerColor(mycolor);
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip]->SetLineColor(mycolor);
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip]->SetLineWidth(2);
	  
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip] = new TGraphErrors();
	  HistName = TString("gr_PhiAveEffic_2orMoreFakes_VsTheta_mom") + TString("_mom") + long(ip+1) + TString("_geo_") + long(igeo+1);
          HistTitle = TString("#phi average Effic_{#geq 2 Fakes} vs ") + MyPolarVariable +  
                      TString(" for ") + MyMomentumVariable + TString(" = ") + TString(ytitle) + TString(" ") + MonUnits + TString(" and geo ") + GeometryList[igeo]->GetName();
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip]->SetName(HistName.Data());
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip]->SetTitle(HistTitle.Data());
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip]->SetMarkerColor(mycolor);
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip]->SetLineColor(mycolor);
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip]->SetLineWidth(2);
	}

      } // end loop over momentum values
      
    } //end loop over parameters
    
  } //end loop over geometries
  

  //////////////////////////////////////////////////////////////////
  /// Looping over geometries and momentum, theta and phi values ///
  //////////////////////////////////////////////////////////////////
  
  fPrintFreq = 10;
  //fPrintFreq = 1;
  for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
    char thetaTitle[100];
    if(polarVariable == TString("theta"))         sprintf(thetaTitle,"theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
    else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    char thetaVal[100];
    if(polarVariable == TString("theta"))         sprintf(thetaVal,"%.1f deg",thetaArr[itheta]/global->GetAngleUnit("deg"));
    else if(polarVariable == TString("costheta")) sprintf(thetaVal,"%.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(thetaVal,"%.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
      for(int iP=0; iP<int(momArr.size()); iP++) { // loop on momentum

        double pi = momArr[iP];
	TVector3 momentum = global->GetMomentum(global->GetMomentumFromMomVar(pi,particle,momVariable),thetaArr[itheta],phiArr[iphi]);
        if(verbose) momentum.Print();

        for(int igeo=0;igeo<NGeometries;igeo++) { //loop over geometries
	  TrackerList[igeo]->SetTrajectoryInitConditions(ParticleOrigin,momentum);
	  
	  Efficiencies_t      Efficiencies;
	  TMatrixD            AveFitCovMatrix;
	  TrackerList[igeo]->GetFitTrackPseudoEfficiency(Efficiencies,AveFitCovMatrix);

	  gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetPoint(gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetN(),pi/global->GetUnit(MonUnits),Efficiencies.Effic_tot*100.0);
	  gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetPoint(gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetN(),pi/global->GetUnit(MonUnits),Efficiencies.Effic_NoFakes*100.0);
	  gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetPoint(gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetN(),pi/global->GetUnit(MonUnits),Efficiencies.Effic_1Fake*100.0);
	  gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetPoint(gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetN(),pi/global->GetUnit(MonUnits),Efficiencies.Effic_2orMoreFakes*100.0);

	  bool IsGoodFit = false;
	  if(Efficiencies.Effic_tot*100 > MinEffic) IsGoodFit = true;
	  
          bool IsNAN = false;
	  for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //loop over track parameters
	    if(IsGoodFit) {
	      if(std::isnan(AveFitCovMatrix(ipar,ipar))) {
		if(verbose) {
	          cout << "WARNNIN:: parameter " << ipar << " (" << titleY[ipar] << ") for " << MyMomentumVariable.Data() << " = " << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << ", "
	               << thetaTitle << ", and phi = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg, Is nan!!!!"
		       << endl;
		}
	        IsNAN = true;
	        break;
	      }
	    }
	  }

	  for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //loop over track parameters
	    //filling up the TGraph with the parameter resolution vs momentum
	    
	    if(IsGoodFit && !IsNAN) { //begin if good fit
	      double y = TrackerList[igeo]->GetTrajectory()->GetParamError(ipar,AveFitCovMatrix);
	      gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetPoint(gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN(),pi/global->GetUnit(MonUnits),y);
	    }
	    else {
	      gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetPoint(gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN(),pi/global->GetUnit(MonUnits),-1);
	    }
	  }  //end of loop over track parameters

          if(verbose) {
	    cout << "igeo = " << igeo+1 << ", tracking efficiency : ";
	    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) { //loop over track parameters
	      cout << titleY[ipar].Data() << " = " << TrackerList[igeo]->GetTrajectory()->GetParamError(ipar,AveFitCovMatrix);
	      if(ipar < TrackerList[igeo]->GetTrajectory()->GetNParameters() - 1) cout << ", ";
	      else                                                                cout << "  ";
              cout << "for " << MyMomentumVariable.Data() << " (bin = " << iP+1 << " of " << momArr.size() << ") = " << pi/global->GetUnit(MonUnits) << " " << MonUnits.Data() << " ";
	      cout << ", " << thetaTitle << "   ";
	      cout << "and phi value = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg   ";
	    }
            global->fWatch.Print();
            global->fWatch.Continue();
          }

          if(verbose || true) {
            if(!((iP+1)%fPrintFreq)) {
              cout << iP+1 << " " << MyMomentumVariable.Data() << " bins processed (out of " << NmomPoints << ") for polar value " << itheta+1 << " (out of " << thetaArr.size() << ")  " << thetaTitle
                   << " , phi value " << iphi+1 << " (out of " << phiArr.size() << ") phi = " << phiArr[iphi]/global->GetAngleUnit("deg") << " deg"
                   << " and geometry " << GeometryList[igeo]->GetName().Data() << " !!!  ";
              global->fWatch.Print();
              global->fWatch.Continue();
            }
          }

        } //end of loop over geometries
      } //end of loop on momentum
    } //end of loop on phi values
  } //end of loop on theta values

  cout << "End of loop in doTrackingPseudoEffic function  ";
  global->fWatch.Print();
  global->fWatch.Continue();
  cout << endl;

  // Check if there was at least one scanned momentum/theta/phi/geometry that 
  // produced a good track parameters covariance matrix
  bool NoGoodCovMatrix = true;
  bool GoodRangeParam[Nmax_params];
  for(int ipar=0;ipar<Nmax_params;ipar++) GoodRangeParam[ipar] = false;

  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
      for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	  for(int i=0;i<gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN();i++) {
	    double x,y;
            gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetPoint(i, x, y);
	    if(y >= 0.0) {
	      NoGoodCovMatrix      = false;
	      GoodRangeParam[ipar] = true;
	      break;
	    }
	  }
	}
      }
    }
  }
  if(NoGoodCovMatrix) {
    cout << endl;
    cout << "Guariguanchi::doTrackingPseudoEfficVsMomentum::WANNING  None of the geometries for the scanned values of momentum/theta/phi produced a good Track param cov matrix." << endl;
    cout << endl;
    return;
  }

  double REffic_tot[2];
  REffic_tot[0] = +1.0e+20;
  REffic_tot[0] = -1.0e+20;
  double REffic_NoFakes[2];
  REffic_NoFakes[0] = +1.0e+20;
  REffic_NoFakes[0] = -1.0e+20;
  double REffic_1Fake[2];
  REffic_1Fake[0] = +1.0e+20;
  REffic_1Fake[0] = -1.0e+20;
  double REffic_2orMoreFakes[2];
  REffic_2orMoreFakes[0] = +1.0e+20;
  REffic_2orMoreFakes[0] = -1.0e+20;
  double RX[2];
  RX[0] = +1.0e+20;
  RX[1] = -1.0e+20;
  double RY[Nmax_params][2];
  for(int ipar=0;ipar<Nmax_params;ipar++) {
    RY[ipar][0] = +1.0e+20;
    RY[ipar][1] = -1.0e+20;
  }

  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
      for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
	for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	  for(int iP=0; iP<gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN(); iP++) {
            double x,y;
            gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetPoint(iP, x, y);
	  
            if(RX[0] > x) RX[0] = x;
            if(RX[1] < x) RX[1] = x;

	    if(y > 0.0) {
              if(RY[ipar][0] > y) RY[ipar][0] = y;
              if(RY[ipar][1] < y) RY[ipar][1] = y;
	    }

	    if(ipar == 0) {
	      gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetPoint(iP, x, y);
	      if(y > MinEffic) {
	        if(REffic_tot[0] > y) REffic_tot[0] = y;
	        if(REffic_tot[1] < y) REffic_tot[1] = y;
	      }
	    
	      gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetPoint(iP, x, y);
	      if(y > MinEffic) {
	        if(REffic_NoFakes[0] > y) REffic_NoFakes[0] = y;
	        if(REffic_NoFakes[1] < y) REffic_NoFakes[1] = y;
	      }
	    
	      gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetPoint(iP, x, y);
	      if(y > MinEffic) {
	        if(REffic_1Fake[0] > y) REffic_1Fake[0] = y;
	        if(REffic_1Fake[1] < y) REffic_1Fake[1] = y;
	      }
	    
	      gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetPoint(iP, x, y);
	      if(y > MinEffic) {
	        if(REffic_2orMoreFakes[0] > y) REffic_2orMoreFakes[0] = y;
	        if(REffic_2orMoreFakes[1] < y) REffic_2orMoreFakes[1] = y;
	      }
	      
	    }
	    
          }
        }
        
      }
    }
  }

  double min_frac = 0.6;
  double porcent,delta,min;
  delta = RX[1] - RX[0];
  porcent = 0.1;
  if(TMath::Abs(delta*global->GetUnit(MonUnits)) < 1.0e-6*global->GetUnit(MonUnits)) {
    delta   = 1.0;
    porcent = 1.0;
  }
  RX[0] -= porcent*delta;
  RX[1] += porcent*delta;
  if(RX[0] < 0.0) RX[0] = 0.0;

  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    delta = RY[ipar][1] - RY[ipar][0];
    if(TMath::Abs(delta) < 1.0e-6) delta = RY[ipar][0];
    porcent = 0.1;
    TString ParErrorUnit = TrackerList[0]->GetTrajectory()->GetParameterErrorUnit(ipar);
    if(TMath::Abs(delta*global->GetUnit(ParErrorUnit)) < 1.0e-6*global->GetUnit(ParErrorUnit)) {
      delta   = 1.0;
      porcent = 1.0;
    }
    
    min = RY[ipar][0];
    RY[ipar][0] -= porcent*delta;
    RY[ipar][1] += porcent*delta;
    if(RY[ipar][0] < 0.0) RY[ipar][0] = min*min_frac;
  }

  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    if(!GoodRangeParam[ipar]) {
      RY[ipar][0] =  1.0e-6;
      RY[ipar][1] =  1.0;
    }
  }
  
  //Efficiency ranges
  delta = REffic_tot[1] - REffic_tot[0];
  porcent = 0.1;
  if(TMath::Abs(delta) < 1.0e-6) {
    delta   = 1.0;
    porcent = 1.0;
  }
  min = REffic_tot[0];
  REffic_tot[0] -= porcent*delta;
  REffic_tot[1] += porcent*delta;
  if(REffic_tot[0] < 0.0) REffic_tot[0] = min*min_frac;
  
  delta = REffic_NoFakes[1] - REffic_NoFakes[0];
  porcent = 0.1;
  if(TMath::Abs(delta) < 1.0e-6) {
    delta   = 1.0;
    porcent = 1.0;
  }
  min = REffic_NoFakes[0];
  REffic_NoFakes[0] -= porcent*delta;
  REffic_NoFakes[1] += porcent*delta;
  if(REffic_NoFakes[0] < 0.0) REffic_NoFakes[0] = min*min_frac;
  
  delta = REffic_1Fake[1] - REffic_1Fake[0];
  porcent = 0.1;
  if(TMath::Abs(delta) < 1.0e-6) {
    delta   = 1.0;
    porcent = 1.0;
  }
  min = REffic_1Fake[0];
  REffic_1Fake[0] -= porcent*delta;
  REffic_1Fake[1] += porcent*delta;
  if(REffic_1Fake[0] < 0.0) REffic_1Fake[0] = min*min_frac;
  
  delta = REffic_2orMoreFakes[1] - REffic_2orMoreFakes[0];
  porcent = 0.1;
  if(TMath::Abs(delta) < 1.0e-6) {
    delta   = 1.0;
    porcent = 1.0;
  }
  min = REffic_2orMoreFakes[0];
  REffic_2orMoreFakes[0] -= porcent*delta;
  REffic_2orMoreFakes[1] += porcent*delta;
  if(REffic_2orMoreFakes[0] < 0.0) REffic_2orMoreFakes[0] = min*min_frac;

  
  //Calculating the average of the resolution on track parameters vs phi
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
	std::vector<double> momList;
	std::vector<double> ParList;
	std::vector<double> Par2List;
	std::vector<int>    pointsList;
	momList.clear();
	ParList.clear();
	Par2List.clear();
	pointsList.clear();
	for(int iphi=0;iphi<NphiPoints;iphi++) {
	  for(int i=0;i<gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN();i++) {
	    double x,y;
	    gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetPoint(i,x,y);
	    //if(y < 0) continue;
	    
	    bool IsInlist = false;
	    for(int imom=0;imom<int(momList.size());imom++) {
	      if(TMath::Abs(x - momList[imom]) < 1.0e-6) {
		IsInlist = true;
		break;
	      }
	    }
	    if(!IsInlist) {
	      momList.push_back(x);
	      ParList.push_back(0);
	      Par2List.push_back(0);
	      pointsList.push_back(0);
	    }
	  }
	}
	
	int counter_mom = 0;
	for(int imom=0;imom<int(momList.size());imom++) {
	  for(int iphi=0;iphi<NphiPoints;iphi++) {
	    double x,y;
	    for(int i=0;i<gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN();i++) {
	      gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetPoint(i,x,y);
	      if(TMath::Abs(momList[imom] - x) < 1.0e-6) break;
	    }
	    if(x == -999.0) continue;
	    if(y < 0.0) continue;
	    ParList[imom]  += y;
	    Par2List[imom] += pow(y,2);
	    pointsList[imom]++;
	  }
	  if(pointsList[imom] > 0) {
	    ParList[imom]  /= pointsList[imom];
	    Par2List[imom] /= pointsList[imom];
	    Par2List[imom] -= pow(ParList[imom],2);
	    Par2List[imom]  = sqrt(Par2List[imom]);
	  }
	  else {
	    ParList[imom] = -1.0;
	    Par2List[imom] = 1.0e-20;
	  }
	  
	  gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetPoint(counter_mom,momList[imom],ParList[imom]);
	  //gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetPointError(counter_mom,1.0e-8,Par2List[imom]/sqrt(pointsList[imom]));

	  counter_mom++;
	}
	
      }
    }
  }
  
  //Calculating the average tracking efficiency vs phi
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      std::vector<double> momList;
      std::vector<double> Effic_tot_List;
      std::vector<double> Effic_tot2_List;
      std::vector<double> Effic_NoFakes_List;
      std::vector<double> Effic_NoFakes2_List;
      std::vector<double> Effic_1Fake_List;
      std::vector<double> Effic_1Fake2_List;
      std::vector<double> Effic_2orMoreFakes_List;
      std::vector<double> Effic_2orMoreFakes2_List;
      std::vector<int>    pointsList;
      momList.clear();
      Effic_tot_List.clear();
      Effic_tot2_List.clear();
      Effic_NoFakes_List.clear();
      Effic_NoFakes2_List.clear();
      Effic_1Fake_List.clear();
      Effic_1Fake2_List.clear();
      Effic_2orMoreFakes_List.clear();
      Effic_2orMoreFakes2_List.clear();
      pointsList.clear();
      for(int iphi=0;iphi<NphiPoints;iphi++) {
	for(int i=0;i<gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	  double x,y;
	  gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	  //if(y < 0) continue;
	  
	  bool IsInlist = false;
	  for(int imom=0;imom<int(momList.size());imom++) {
	    if(TMath::Abs(x - momList[imom]) < 1.0e-6) {
	      IsInlist = true;
	      break;
	    }
	  }
	  if(!IsInlist) {
	    momList.push_back(x);
	    Effic_tot_List.push_back(0);
	    Effic_tot2_List.push_back(0);
	    Effic_NoFakes_List.push_back(0);
	    Effic_NoFakes2_List.push_back(0);
	    Effic_1Fake_List.push_back(0);
	    Effic_1Fake2_List.push_back(0);
	    Effic_2orMoreFakes_List.push_back(0);
	    Effic_2orMoreFakes2_List.push_back(0);
	    pointsList.push_back(0);
	  }
	}
      }
      
      int counter_mom = 0;
      for(int imom=0;imom<int(momList.size());imom++) {
	for(int iphi=0;iphi<NphiPoints;iphi++) {
	  double x,y_tot,y_NoFakes,y_1Fake,y_2orMoreFakes;
	  for(int i=0;i<gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	    gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y_tot);
	    gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y_NoFakes);
	    gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y_1Fake);
	    gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y_2orMoreFakes);
	    if(TMath::Abs(momList[imom] - x) < 1.0e-6) break;
	  }
	  if(x == -999.0) continue;
	  if(y_tot < 0.0) continue;
	  
	  Effic_tot_List[imom]           += y_tot;
	  Effic_tot2_List[imom]          += pow(y_tot,2);
	  Effic_NoFakes_List[imom]       += y_NoFakes;
	  Effic_NoFakes2_List[imom]      += pow(y_NoFakes,2);
	  Effic_1Fake_List[imom]         += y_1Fake;
	  Effic_1Fake2_List[imom]        += pow(y_1Fake,2);
	  Effic_2orMoreFakes_List[imom]  += y_2orMoreFakes;
	  Effic_2orMoreFakes2_List[imom] += pow(y_2orMoreFakes,2);
	  
	  pointsList[imom]++;
	}
	if(pointsList[imom] > 0) {
	  Effic_tot_List[imom]           /= pointsList[imom];
	  Effic_tot2_List[imom]          /= pointsList[imom];
	  Effic_tot2_List[imom]          -= pow(Effic_tot_List[imom],2);
	  Effic_tot2_List[imom]           = sqrt(Effic_tot2_List[imom]);
	  
	  Effic_NoFakes_List[imom]       /= pointsList[imom];
	  Effic_NoFakes2_List[imom]      /= pointsList[imom];
	  Effic_NoFakes2_List[imom]      -= pow(Effic_NoFakes_List[imom],2);
	  Effic_NoFakes2_List[imom]       = sqrt(Effic_NoFakes2_List[imom]);
	  
	  Effic_1Fake_List[imom]         /= pointsList[imom];
	  Effic_1Fake2_List[imom]        /= pointsList[imom];
	  Effic_1Fake2_List[imom]        -= pow(Effic_1Fake_List[imom],2);
	  Effic_1Fake2_List[imom]         = sqrt(Effic_1Fake2_List[imom]);
	  
	  Effic_2orMoreFakes_List[imom]  /= pointsList[imom];
	  Effic_2orMoreFakes2_List[imom] /= pointsList[imom];
	  Effic_2orMoreFakes2_List[imom] -= pow(Effic_2orMoreFakes_List[imom],2);
	  Effic_2orMoreFakes2_List[imom]  = sqrt(Effic_2orMoreFakes2_List[imom]);
	}
	else {
	  Effic_tot_List[imom] = -1.0;
	  Effic_tot2_List[imom] = 1.0e-20;
	  
	  Effic_NoFakes_List[imom] = -1.0;
	  Effic_NoFakes2_List[imom] = 1.0e-20;
	  
	  Effic_1Fake_List[imom] = -1.0;
	  Effic_1Fake2_List[imom] = 1.0e-20;
	  
	  Effic_2orMoreFakes_List[imom] = -1.0;
	  Effic_2orMoreFakes2_List[imom] = 1.0e-20;
	}
	
	gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetPoint(counter_mom,momList[imom],Effic_tot_List[imom]);
	//gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetPointError(counter_mom,1.0e-8,Effic_tot2_List[imom]/sqrt(pointsList[imom]));
	
	gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetPoint(counter_mom,momList[imom],Effic_NoFakes_List[imom]);
	//gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetPointError(counter_mom,1.0e-8,Effic_NoFakes2_List[imom]/sqrt(pointsList[imom]));
	
	gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetPoint(counter_mom,momList[imom],Effic_1Fake_List[imom]);
	//gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetPointError(counter_mom,1.0e-8,Effic_1Fake2_List[imom]/sqrt(pointsList[imom]));
	
	gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetPoint(counter_mom,momList[imom],Effic_2orMoreFakes_List[imom]);
	//gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetPointError(counter_mom,1.0e-8,Effic_2orMoreFakes2_List[imom]/sqrt(pointsList[imom]));
	
	counter_mom++;
      }
       
    }
  }
  
  double Rtheta[2];
  double Rphi[2];
  Rtheta[0] = Rphi[0] = +1.0e+20;
  Rtheta[1] = Rphi[1] = -1.0e+20;
  for(int itheta=0;itheta<NthetaPoints;itheta++) { //begin loop over theta
    double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
    if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
    else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
    else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
    
    if(Rtheta[0] > thetaValue) Rtheta[0] = thetaValue;
    if(Rtheta[1] < thetaValue) Rtheta[1] = thetaValue;
    
    for(int iphi=0;iphi<NphiPoints;iphi++) { //loop on phi values
      if(Rphi[0] > phiArr[iphi]) Rphi[0] = phiArr[iphi];
      if(Rphi[1] < phiArr[iphi]) Rphi[1] = phiArr[iphi];
    } //end loop on phi values
  }
  
  
  //Get the track parameter resolution vs theta averaging over phi values
  double Rpar_vs_theta[Nmax_params][2];
  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters() ;ipar++) {
    Rpar_vs_theta[ipar][0] = +1.0e+20;
    Rpar_vs_theta[ipar][1] = -1.0e+20;
  }
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int ipar=0;ipar<TrackerList[igeo]->GetTrajectory()->GetNParameters();ipar++) {
      for(int iP=0; iP<MyNbins_mom_redu; iP++) {
	double pi;
	if(EfficAnalysisPars[0].UseAllMomVals) pi = momArr[iP];
	else                                   pi = momArr[0] + (iP+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
	pi /= global->GetUnit(MonUnits);
	int counter_theta = 0;
	for(int itheta=0;itheta<NthetaPoints;itheta++) {
	  double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
          if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
          else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
          else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
	  
	  //if(gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetN() == 0) continue;
	  int idx = -999;
	  double mom1,mom2,y1,y2;
	  for(int i=0;i<gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetN()-1;i++) {
	    gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetPoint(i,  mom1,y1);
	    gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetPoint(i+1,mom2,y2);
	    if(pi >= mom1 && pi <= mom2) {
	      idx = i;
	      break;
	    }
	  }
	  if(idx == -999) continue;
	  
	  double a = (y2-y1)/(mom2 - mom1);
	  double b = y2 - a*mom2;
	  
	  double value = a*pi + b;
	  gr_PhiAveSigmaParVsTheta[igeo][iP][ipar]->SetPoint(counter_theta,thetaValue,value);
	  if(Rpar_vs_theta[ipar][0] > value) Rpar_vs_theta[ipar][0] = value;
	  if(Rpar_vs_theta[ipar][1] < value) Rpar_vs_theta[ipar][1] = value;
	  
	  counter_theta++;
	}
      }
    }
  }
  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    delta = Rpar_vs_theta[ipar][1] - Rpar_vs_theta[ipar][0];
    if(TMath::Abs(delta) < 1.0e-6) delta = Rpar_vs_theta[ipar][0];
    porcent = 0.1;
    
    min = Rpar_vs_theta[ipar][0];
    Rpar_vs_theta[ipar][0] -= porcent*delta;
    Rpar_vs_theta[ipar][1] += porcent*delta;
    if(Rpar_vs_theta[ipar][0] < 0.0) Rpar_vs_theta[ipar][0] = min*min_frac;
  }
  
  //Get the tracking efficiency vs theta averaging over phi values
  double REffic_tot_vs_theta[2];
  REffic_tot_vs_theta[0] = +1.0e+20;
  REffic_tot_vs_theta[1] = -1.0e+20;
  double REffic_NoFakes_vs_theta[2];
  REffic_NoFakes_vs_theta[0] = +1.0e+20;
  REffic_NoFakes_vs_theta[1] = -1.0e+20;
  double REffic_1Fake_vs_theta[2];
  REffic_1Fake_vs_theta[0] = +1.0e+20;
  REffic_1Fake_vs_theta[1] = -1.0e+20;
  double REffic_2orMoreFakes_vs_theta[2];
  REffic_2orMoreFakes_vs_theta[0] = +1.0e+20;
  REffic_2orMoreFakes_vs_theta[1] = -1.0e+20;
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int iP=0; iP<MyNbins_mom_redu; iP++) {
      double pi;
      if(EfficAnalysisPars[0].UseAllMomVals) pi = momArr[iP];
      else                                   pi = momArr[0] + (iP+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
      pi /= global->GetUnit(MonUnits);
      int counter_theta = 0;
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
	double thetaValue = thetaArr[itheta]/global->GetUnit("deg"); 
	if(polarVariable == TString("theta"))         thetaValue = thetaArr[itheta]/global->GetUnit("deg");
	else if(polarVariable == TString("costheta")) thetaValue = global->FromThetaToCosTheta(thetaArr[itheta]);
	else if(polarVariable == TString("eta"))      thetaValue = global->FromThetaToEta(thetaArr[itheta]);
	
	//if(gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetN() == 0) continue;
	int idx = -999;
	double mom1,mom2,y1_tot,y2_tot,y1_NoFakes,y2_NoFakes,y1_1Fake,y2_1Fake,y1_2orMoreFakes,y2_2orMoreFakes;
	for(int i=0;i<gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetN()-1;i++) {
	  gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,  mom1,y1_tot);
	  gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i+1,mom2,y2_tot);
	  gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,  mom1,y1_NoFakes);
	  gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i+1,mom2,y2_NoFakes);
	  gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,  mom1,y1_1Fake);
	  gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i+1,mom2,y2_1Fake);
	  gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,  mom1,y1_2orMoreFakes);
	  gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i+1,mom2,y2_2orMoreFakes);
	  if(pi >= mom1 && pi <= mom2) {
	    idx = i;
	    break;
	  }
	}
	if(idx == -999) continue;
	
	double a,b,value;
	
	a = (y2_tot-y1_tot)/(mom2 - mom1);
	b = y2_tot - a*mom2;
	value = a*pi + b;
	gr_PhiAveEffic_tot_VsTheta[igeo][iP]->SetPoint(counter_theta,thetaValue,value);
	if(value > MinEffic) {
	  if(REffic_tot_vs_theta[0] > value) REffic_tot_vs_theta[0] = value;
          if(REffic_tot_vs_theta[1] < value) REffic_tot_vs_theta[1] = value;
	}
	
	a = (y2_NoFakes-y1_NoFakes)/(mom2 - mom1);
	b = y2_NoFakes - a*mom2;
	value = a*pi + b;
	gr_PhiAveEffic_NoFakes_VsTheta[igeo][iP]->SetPoint(counter_theta,thetaValue,value);
	if(value > MinEffic) {
	  if(REffic_NoFakes_vs_theta[0] > value) REffic_NoFakes_vs_theta[0] = value;
          if(REffic_NoFakes_vs_theta[1] < value) REffic_NoFakes_vs_theta[1] = value;
	}
	
	a = (y2_1Fake-y1_1Fake)/(mom2 - mom1);
	b = y2_1Fake - a*mom2;
	value = a*pi + b;
	gr_PhiAveEffic_1Fake_VsTheta[igeo][iP]->SetPoint(counter_theta,thetaValue,value);
	if(value > MinEffic) {
	  if(REffic_1Fake_vs_theta[0] > value) REffic_1Fake_vs_theta[0] = value;
          if(REffic_1Fake_vs_theta[1] < value) REffic_1Fake_vs_theta[1] = value;
	}
	
	a = (y2_2orMoreFakes-y1_2orMoreFakes)/(mom2 - mom1);
	b = y2_2orMoreFakes - a*mom2;
	value = a*pi + b;
	gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][iP]->SetPoint(counter_theta,thetaValue,value);
	if(value > MinEffic) {
	  if(REffic_2orMoreFakes_vs_theta[0] > value) REffic_2orMoreFakes_vs_theta[0] = value;
          if(REffic_2orMoreFakes_vs_theta[1] < value) REffic_2orMoreFakes_vs_theta[1] = value;
	}

	counter_theta++;
      }
    }
  }
  
  delta = REffic_tot_vs_theta[1] - REffic_tot_vs_theta[0];
  if(TMath::Abs(delta) < 1.0e-6) delta = REffic_tot_vs_theta[0];
  porcent = 0.1;  
  min = REffic_tot_vs_theta[0];
  REffic_tot_vs_theta[0] -= porcent*delta;
  REffic_tot_vs_theta[1] += porcent*delta;
  if(REffic_tot_vs_theta[0] < 0.0) REffic_tot_vs_theta[0] = min*min_frac;
  
  delta = REffic_NoFakes_vs_theta[1] - REffic_NoFakes_vs_theta[0];
  if(TMath::Abs(delta) < 1.0e-6) delta = REffic_NoFakes_vs_theta[0];
  porcent = 0.1;  
  min = REffic_NoFakes_vs_theta[0];
  REffic_NoFakes_vs_theta[0] -= porcent*delta;
  REffic_NoFakes_vs_theta[1] += porcent*delta;
  if(REffic_NoFakes_vs_theta[0] < 0.0) REffic_NoFakes_vs_theta[0] = min*min_frac;
  
  delta = REffic_1Fake_vs_theta[1] - REffic_1Fake_vs_theta[0];
  if(TMath::Abs(delta) < 1.0e-6) delta = REffic_1Fake_vs_theta[0];
  porcent = 0.1;  
  min = REffic_1Fake_vs_theta[0];
  REffic_1Fake_vs_theta[0] -= porcent*delta;
  REffic_1Fake_vs_theta[1] += porcent*delta;
  if(REffic_1Fake_vs_theta[0] < 0.0) REffic_1Fake_vs_theta[0] = min*min_frac;
  
  delta = REffic_2orMoreFakes_vs_theta[1] - REffic_2orMoreFakes_vs_theta[0];
  if(TMath::Abs(delta) < 1.0e-6) delta = REffic_2orMoreFakes_vs_theta[0];
  porcent = 0.1;  
  min = REffic_2orMoreFakes_vs_theta[0];
  REffic_2orMoreFakes_vs_theta[0] -= porcent*delta;
  REffic_2orMoreFakes_vs_theta[1] += porcent*delta;
  if(REffic_2orMoreFakes_vs_theta[0] < 0.0) REffic_2orMoreFakes_vs_theta[0] = min*min_frac;
  
  delta = Rtheta[1] - Rtheta[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Rtheta[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = 1.0;
  porcent = 0.1;
  Rtheta[0] -= porcent*delta;
  Rtheta[1] += porcent*delta;
  if(Rtheta[0] < 0.0) Rtheta[0] = 0.0;
  
  delta = Rphi[1] - Rphi[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = Rphi[0];
  if(TMath::Abs(delta) < 1.0e-8) delta = 1.0;
  porcent = 0.1;
  Rphi[0] -= porcent*delta;
  Rphi[1] += porcent*delta;
  //if(Rphi[0] < 0.0) Rphi[0] = 0.0;
  
  Rphi[0]      /= global->GetAngleUnit("deg");
  Rphi[1]      /= global->GetAngleUnit("deg");

  TString MyPolarVarName;
  TString MyPolarVarNameWithUnits;
  if(polarVariable == TString("theta")) {
    MyPolarVarName          = TString("#theta");
    MyPolarVarNameWithUnits = TString("#theta (") + ThetaUnits + TString(")");
  }
  else if(polarVariable == TString("costheta")) {
    MyPolarVarName          = TString("cos(#theta)");
    MyPolarVarNameWithUnits = TString("cos(#theta)");
  }
  else if(polarVariable == TString("eta")) {
    MyPolarVarName          = TString("#eta");
    MyPolarVarNameWithUnits = TString("#eta");
  }
  
  //Total Efficiency frame histos  
  HistName = TString("href_Effic_tot");
  TH1F* href_Effic_tot = new TH1F(HistName.Data(),
				  "Total Efficiency",
				  100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_Effic_tot->SetXTitle(HistName.Data());
  href_Effic_tot->GetXaxis()->CenterTitle(true);
  href_Effic_tot->SetYTitle("Efficiency (%)");
  href_Effic_tot->GetYaxis()->CenterTitle(true);
  href_Effic_tot->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_tot->SetLineColor(1);
  href_Effic_tot->SetLineWidth(1);
  href_Effic_tot->SetMinimum(REffic_tot[0]);
  href_Effic_tot->SetMaximum(REffic_tot[1]);
  href_Effic_tot->SetNdivisions(5 + 100*5,"X");
  href_Effic_tot->SetNdivisions(5 + 100*5,"Y");
  href_Effic_tot->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_tot->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_tot->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_tot->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_tot->SetStats(false);
    
  HistName  = TString("href_Effic_tot_vs_theta");
  HistTitle = TString("Total Efficiency vs ") + MyPolarVarName;
  TH1F* href_Effic_tot_vs_theta = new TH1F(HistName.Data(),
					   HistTitle.Data(),
					   100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_Effic_tot_vs_theta->SetXTitle(HistName.Data());
  HistName = TString("#phi Ave Efficiency (%)");
  href_Effic_tot_vs_theta->GetXaxis()->CenterTitle(true);
  href_Effic_tot_vs_theta->SetYTitle(HistName.Data());
  href_Effic_tot_vs_theta->GetYaxis()->CenterTitle(true);
  href_Effic_tot_vs_theta->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_tot_vs_theta->SetLineColor(1);
  href_Effic_tot_vs_theta->SetLineWidth(1);
  href_Effic_tot_vs_theta->SetMinimum(REffic_tot_vs_theta[0]);
  href_Effic_tot_vs_theta->SetMaximum(REffic_tot_vs_theta[1]);
  href_Effic_tot_vs_theta->SetNdivisions(5 + 100*5,"X");
  href_Effic_tot_vs_theta->SetNdivisions(5 + 100*5,"Y");
  href_Effic_tot_vs_theta->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_tot_vs_theta->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_tot_vs_theta->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_tot_vs_theta->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_tot_vs_theta->SetStats(false);
  
  //NoFakes Efficiency frame histos  
  HistName = TString("href_Effic_NoFakes");
  TH1F* href_Effic_NoFakes = new TH1F(HistName.Data(),
				      "No-Fakes Efficiency",
				      100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_Effic_NoFakes->SetXTitle(HistName.Data());
  href_Effic_NoFakes->GetXaxis()->CenterTitle(true);
  href_Effic_NoFakes->SetYTitle("Efficiency (%)");
  href_Effic_NoFakes->GetYaxis()->CenterTitle(true);
  href_Effic_NoFakes->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_NoFakes->SetLineColor(1);
  href_Effic_NoFakes->SetLineWidth(1);
  href_Effic_NoFakes->SetMinimum(REffic_NoFakes[0]);
  href_Effic_NoFakes->SetMaximum(REffic_NoFakes[1]);
  href_Effic_NoFakes->SetNdivisions(5 + 100*5,"X");
  href_Effic_NoFakes->SetNdivisions(5 + 100*5,"Y");
  href_Effic_NoFakes->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_NoFakes->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_NoFakes->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_NoFakes->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_NoFakes->SetStats(false);
    
  HistName  = TString("href_Effic_NoFakes_vs_theta");
  HistTitle = TString("No-Fakes Efficiency vs ") + MyPolarVarName;
  TH1F* href_Effic_NoFakes_vs_theta = new TH1F(HistName.Data(),
					       HistTitle.Data(),
					       100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_Effic_NoFakes_vs_theta->SetXTitle(HistName.Data());
  HistName = TString("#phi Ave Efficiency (%)");
  href_Effic_NoFakes_vs_theta->GetXaxis()->CenterTitle(true);
  href_Effic_NoFakes_vs_theta->SetYTitle(HistName.Data());
  href_Effic_NoFakes_vs_theta->GetYaxis()->CenterTitle(true);
  href_Effic_NoFakes_vs_theta->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_NoFakes_vs_theta->SetLineColor(1);
  href_Effic_NoFakes_vs_theta->SetLineWidth(1);
  href_Effic_NoFakes_vs_theta->SetMinimum(REffic_NoFakes_vs_theta[0]);
  href_Effic_NoFakes_vs_theta->SetMaximum(REffic_NoFakes_vs_theta[1]);
  href_Effic_NoFakes_vs_theta->SetNdivisions(5 + 100*5,"X");
  href_Effic_NoFakes_vs_theta->SetNdivisions(5 + 100*5,"Y");
  href_Effic_NoFakes_vs_theta->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_NoFakes_vs_theta->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_NoFakes_vs_theta->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_NoFakes_vs_theta->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_NoFakes_vs_theta->SetStats(false);
  
  //1Fake Efficiency frame histos
  HistName = TString("href_Effic_1Fake");
  TH1F* href_Effic_1Fake = new TH1F(HistName.Data(),
				    "1-Fake Efficiency",
				    100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_Effic_1Fake->SetXTitle(HistName.Data());
  href_Effic_1Fake->GetXaxis()->CenterTitle(true);
  href_Effic_1Fake->SetYTitle("Efficiency (%)");
  href_Effic_1Fake->GetYaxis()->CenterTitle(true);
  href_Effic_1Fake->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_1Fake->SetLineColor(1);
  href_Effic_1Fake->SetLineWidth(1);
  href_Effic_1Fake->SetMinimum(REffic_1Fake[0]);
  href_Effic_1Fake->SetMaximum(REffic_1Fake[1]);
  href_Effic_1Fake->SetNdivisions(5 + 100*5,"X");
  href_Effic_1Fake->SetNdivisions(5 + 100*5,"Y");
  href_Effic_1Fake->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_1Fake->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_1Fake->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_1Fake->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_1Fake->SetStats(false);
    
  HistName  = TString("href_Effic_1Fake_vs_theta");
  HistTitle = TString("1-Fake Efficiency vs ") + MyPolarVarName;
  TH1F* href_Effic_1Fake_vs_theta = new TH1F(HistName.Data(),
					     HistTitle.Data(),
					     100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_Effic_1Fake_vs_theta->SetXTitle(HistName.Data());
  HistName = TString("#phi Ave Efficiency (%)");
  href_Effic_1Fake_vs_theta->GetXaxis()->CenterTitle(true);
  href_Effic_1Fake_vs_theta->SetYTitle(HistName.Data());
  href_Effic_1Fake_vs_theta->GetYaxis()->CenterTitle(true);
  href_Effic_1Fake_vs_theta->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_1Fake_vs_theta->SetLineColor(1);
  href_Effic_1Fake_vs_theta->SetLineWidth(1);
  href_Effic_1Fake_vs_theta->SetMinimum(REffic_1Fake_vs_theta[0]);
  href_Effic_1Fake_vs_theta->SetMaximum(REffic_1Fake_vs_theta[1]);
  href_Effic_1Fake_vs_theta->SetNdivisions(5 + 100*5,"X");
  href_Effic_1Fake_vs_theta->SetNdivisions(5 + 100*5,"Y");
  href_Effic_1Fake_vs_theta->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_1Fake_vs_theta->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_1Fake_vs_theta->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_1Fake_vs_theta->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_1Fake_vs_theta->SetStats(false);
  
  //2 or more Fakes Efficiency frame histos
  HistName = TString("href_Effic_2orMoreFakes");
  TH1F* href_Effic_2orMoreFakes = new TH1F(HistName.Data(),
					   "2 #geq Fakes Efficiency",
					   100,RX[0],RX[1]);
  HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
  href_Effic_2orMoreFakes->SetXTitle(HistName.Data());
  href_Effic_2orMoreFakes->GetXaxis()->CenterTitle(true);
  href_Effic_2orMoreFakes->SetYTitle("Efficiency (%)");
  href_Effic_2orMoreFakes->GetYaxis()->CenterTitle(true);
  href_Effic_2orMoreFakes->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_2orMoreFakes->SetLineColor(1);
  href_Effic_2orMoreFakes->SetLineWidth(1);
  href_Effic_2orMoreFakes->SetMinimum(REffic_2orMoreFakes[0]);
  href_Effic_2orMoreFakes->SetMaximum(REffic_2orMoreFakes[1]);
  href_Effic_2orMoreFakes->SetNdivisions(5 + 100*5,"X");
  href_Effic_2orMoreFakes->SetNdivisions(5 + 100*5,"Y");
  href_Effic_2orMoreFakes->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_2orMoreFakes->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_2orMoreFakes->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_2orMoreFakes->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_2orMoreFakes->SetStats(false);
    
  HistName  = TString("href_Effic_2orMoreFakes_vs_theta");
  HistTitle = TString("2 $geq Fakes Efficiency vs ") + MyPolarVarName;
  TH1F* href_Effic_2orMoreFakes_vs_theta = new TH1F(HistName.Data(),
						    HistTitle.Data(),
						    100,Rtheta[0],Rtheta[1]);
  HistName = MyPolarVarNameWithUnits;
  href_Effic_2orMoreFakes_vs_theta->SetXTitle(HistName.Data());
  HistName = TString("#phi Ave Efficiency (%)");
  href_Effic_2orMoreFakes_vs_theta->GetXaxis()->CenterTitle(true);
  href_Effic_2orMoreFakes_vs_theta->SetYTitle(HistName.Data());
  href_Effic_2orMoreFakes_vs_theta->GetYaxis()->CenterTitle(true);
  href_Effic_2orMoreFakes_vs_theta->GetYaxis()->SetTitleOffset(TitleOffSet);
  href_Effic_2orMoreFakes_vs_theta->SetLineColor(1);
  href_Effic_2orMoreFakes_vs_theta->SetLineWidth(1);
  href_Effic_2orMoreFakes_vs_theta->SetMinimum(REffic_2orMoreFakes_vs_theta[0]);
  href_Effic_2orMoreFakes_vs_theta->SetMaximum(REffic_2orMoreFakes_vs_theta[1]);
  href_Effic_2orMoreFakes_vs_theta->SetNdivisions(5 + 100*5,"X");
  href_Effic_2orMoreFakes_vs_theta->SetNdivisions(5 + 100*5,"Y");
  href_Effic_2orMoreFakes_vs_theta->GetXaxis()->SetTitleSize(TheSize);
  href_Effic_2orMoreFakes_vs_theta->GetXaxis()->SetLabelSize(TheSize);
  href_Effic_2orMoreFakes_vs_theta->GetYaxis()->SetTitleSize(TheSize);
  href_Effic_2orMoreFakes_vs_theta->GetYaxis()->SetLabelSize(TheSize);
  href_Effic_2orMoreFakes_vs_theta->SetStats(false);
  
  TH1F* href[Nmax_params];
  TH1F* href_vs_theta[Nmax_params];
  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
    HistName = TString("href_par") + long(ipar);
    href[ipar] = new TH1F(HistName.Data(),title[ipar].Data(),
			  100,RX[0],RX[1]);
    HistName = MyMomentumVariable + TString(" (") + MonUnits + TString(")");
    href[ipar]->SetXTitle(HistName.Data());
    href[ipar]->GetXaxis()->CenterTitle(true);
    href[ipar]->SetYTitle(titleY[ipar].Data());
    href[ipar]->GetYaxis()->CenterTitle(true);
    href[ipar]->GetYaxis()->SetTitleOffset(TitleOffSet);
    href[ipar]->SetLineColor(1);
    href[ipar]->SetLineWidth(1);
    href[ipar]->SetMinimum(RY[ipar][0]);
    href[ipar]->SetMaximum(RY[ipar][1]);
    href[ipar]->SetNdivisions(5 + 100*5,"X");
    href[ipar]->SetNdivisions(5 + 100*5,"Y");
    href[ipar]->GetXaxis()->SetTitleSize(TheSize);
    href[ipar]->GetXaxis()->SetLabelSize(TheSize);
    href[ipar]->GetYaxis()->SetTitleSize(TheSize);
    href[ipar]->GetYaxis()->SetLabelSize(TheSize);
    href[ipar]->SetStats(false);
    
    HistName  = TString("href_par") + long(ipar) + TString("_vs_theta");
    HistTitle = titleY[ipar] + TString(" vs ") + MyPolarVarName;
    href_vs_theta[ipar] = new TH1F(HistName.Data(),
				   HistTitle.Data(),
				   100,Rtheta[0],Rtheta[1]);
    HistName = MyPolarVarNameWithUnits;
    href_vs_theta[ipar]->SetXTitle(HistName.Data());
    HistName = TString("#phi Ave ") + titleY[ipar];
    href_vs_theta[ipar]->GetXaxis()->CenterTitle(true);
    href_vs_theta[ipar]->SetYTitle(HistName.Data());
    href_vs_theta[ipar]->GetYaxis()->CenterTitle(true);
    href_vs_theta[ipar]->GetYaxis()->SetTitleOffset(TitleOffSet);
    href_vs_theta[ipar]->SetLineColor(1);
    href_vs_theta[ipar]->SetLineWidth(1);
    href_vs_theta[ipar]->SetMinimum(Rpar_vs_theta[ipar][0]);
    href_vs_theta[ipar]->SetMaximum(Rpar_vs_theta[ipar][1]);
    href_vs_theta[ipar]->SetNdivisions(5 + 100*5,"X");
    href_vs_theta[ipar]->SetNdivisions(5 + 100*5,"Y");
    href_vs_theta[ipar]->GetXaxis()->SetTitleSize(TheSize);
    href_vs_theta[ipar]->GetXaxis()->SetLabelSize(TheSize);
    href_vs_theta[ipar]->GetYaxis()->SetTitleSize(TheSize);
    href_vs_theta[ipar]->GetYaxis()->SetLabelSize(TheSize);
    href_vs_theta[ipar]->SetStats(false);
  }

  
  TLegend* leg = NULL;
  leg = new TLegend(0.10,0.10,0.90,0.80);
  leg->SetFillColor(10);
  
  TString EPSName  = TheOutputFile + TString("_") + long(GlobalFileCounter) + TString(".eps");
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");

  TCanvas* c1 = new TCanvas("c1","c1",LargeCanvasX,LargeCanvasY);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.20);
  c1->SetBottomMargin(0.15);
  //c1->SetRightMargin(0.10);

  c1->Print(EPSNameO.Data());
  
  bool AlsoEfficLogY = false;

  for(int igeo=0;igeo<NGeometries;igeo++) leg->AddEntry(gr_res_ave_vs_mom[igeo][0][0][0],GeometryList[igeo]->GetName().Data(),"lp");

  double R[2];
  for(int itheta=0;itheta<NthetaPoints;itheta++) { //loop on theta values
    char thetaTitle[300];
    if(polarVariable == TString("theta"))         sprintf(thetaTitle,"#theta = %.1f deg",thetaArr[itheta]/global->GetAngleUnit(ThetaUnits));
    else if(polarVariable == TString("costheta")) sprintf(thetaTitle,"cos(#theta) = %.3f",global->FromThetaToCosTheta(thetaArr[itheta]));
    else if(polarVariable == TString("eta"))      sprintf(thetaTitle,"#eta = %.3f",global->FromThetaToEta(thetaArr[itheta]));
    
    for(int iphi=0;iphi<NphiPoints;iphi++) {
      char phiTitle[300];
      sprintf(phiTitle,"#phi = %.1f %s",phiArr[iphi]/global->GetUnit(ThetaUnits),ThetaUnits.Data());
      
      TLine* l_ref_mom[Nmax_params];
      
      if(!EfficAnalysisPars[0].PlotOnlyPhiAveraged) {
	c1->Clear();
        c1->Divide(3,2);
	c1->cd(1);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
	R[0] = +1.0e+20;
	R[1] = -1.0e+20;
	if(!EfficAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < MinEffic) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
            if(TMath::Abs(delta) < 1.0e-4) delta = R[0];
            porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	    
	  href_Effic_tot->SetMinimum(R[0]);
	  href_Effic_tot->SetMaximum(R[1]);
	}
	
        href_Effic_tot->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_effic_tot_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_effic_tot_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_effic_tot_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
	}
	
        c1->cd(2);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
	R[0] = +1.0e+20;
	R[1] = -1.0e+20;
	if(!EfficAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < MinEffic) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
            if(TMath::Abs(delta) < 1.0e-4) delta = R[0];
            porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	    
	  href_Effic_NoFakes->SetMinimum(R[0]);
	  href_Effic_NoFakes->SetMaximum(R[1]);
	}
	
        href_Effic_NoFakes->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
	}
	
	c1->cd(4);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
	R[0] = +1.0e+20;
	R[1] = -1.0e+20;
	if(!EfficAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < MinEffic) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
            if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
            porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	    
	  href_Effic_1Fake->SetMinimum(R[0]);
	  href_Effic_1Fake->SetMaximum(R[1]);
	}
	
        href_Effic_1Fake->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
	}
	
	c1->cd(5);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
	R[0] = +1.0e+20;
	R[1] = -1.0e+20;
	if(!EfficAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetN();i++) {
	      double x,y;
	      gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetPoint(i,x,y);
	      if(y < MinEffic) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
            if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
            porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	    
	  href_Effic_2orMoreFakes->SetMinimum(R[0]);
	  href_Effic_2orMoreFakes->SetMaximum(R[1]);
	}
	
        href_Effic_2orMoreFakes->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
          if(gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetN() > 0) {
	    if(gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetN() == 1) {
	      gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetMarkerColor(gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->GetLineColor());
	      gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->Draw("PEL");
	  }
	}

        c1->cd(6);
        HistName  = TString(thetaTitle) + TString(", ") + TString(phiTitle);
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.10,0.88,HistName.Data());
        leg->Draw("same");

	c1->cd();
	PlotLogo(0.15,0.89,0.35,0.5);
        c1->Print(EPSName.Data());

        c1->Clear();
        c1->Divide(3,2);
        int counter_pad = 0;
        for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
          counter_pad++;
          if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++; 
          c1->cd(counter_pad);
          gPad->SetFillColor(10);
          gPad->SetFrameFillColor(10);
          gPad->SetTickx(1);
          gPad->SetTicky(1);
          gPad->SetLeftMargin(0.20);
          gPad->SetBottomMargin(0.15);
	  gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes);
	  R[0] = +1.0e+20;
	  R[1] = -1.0e+20;
	  if(!EfficAnalysisPars[0].SameRange) {
	    for(int igeo=0;igeo<NGeometries;igeo++) {
	      for(int i=0;i<gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN();i++) {
	        double x,y;
	        gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetPoint(i,x,y);
		if(y < 0.0) continue;
		
	        if(R[0] > y) R[0] = y;
	        if(R[1] < y) R[1] = y;
	      }
	    }
	    bool NoRange = false;
	    if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	    if(NoRange) {
	      R[0] = 1.0e-6;
	      R[1] = 1.0;
	    }
	    else {
	      delta = R[1] - R[0];
              if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
              porcent = 0.1;
	      min = R[0];
	      R[0] -= porcent*delta;
	      R[1] += porcent*delta;
	      if(R[0] < 0.0) R[0] = min*min_frac;
	    }
	    
	    href[ipar]->SetMinimum(R[0]);
	    href[ipar]->SetMaximum(R[1]);
	  }
	  l_ref_mom[ipar] = new TLine(momArr[0]/global->GetUnit(MonUnits),href[ipar]->GetMinimum(),
				      momArr[0]/global->GetUnit(MonUnits),href[ipar]->GetMaximum());
	  l_ref_mom[ipar]->SetLineColor(2);
	  l_ref_mom[ipar]->SetLineWidth(1);
	  l_ref_mom[ipar]->SetLineStyle(2);
	
          href[ipar]->Draw();
          for(int igeo=0;igeo<NGeometries;igeo++) {
            if(gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN() > 0) {
	      if(gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetN() == 1) {
		gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetMarkerColor(gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->GetLineColor());
		gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetMarkerStyle(SinglePointMarkerSytle);
		gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->SetMarkerSize(SinglePointMarkerSize);
	      }
	      gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->Draw("PEL");
	    }
          }
          l_ref_mom[ipar]->Draw();
        }
        c1->cd(6);
        HistName  = TString(thetaTitle) + TString(", ") + TString(phiTitle);
        latex->SetTextColor(kBlack);
        latex->DrawLatex(0.10,0.88,HistName.Data());
        leg->Draw("same");
	
	c1->cd();
	PlotLogo(0.15,0.89,0.35,0.5);
        c1->Print(EPSName.Data());
      }
    }
    
    if(EfficAnalysisPars[0].PlotOnlyPhiAveraged || EfficAnalysisPars[0].PlotPhiAveraged) {
      c1->Clear();
      c1->Divide(3,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      TString Titles_tot;
      TString YTitles_tot;

      Titles_tot = TString(href_Effic_tot->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_tot->GetTitle());
      href_Effic_tot->SetTitle(HistName.Data());
      
      YTitles_tot = TString(href_Effic_tot->GetYaxis()->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_tot->GetYaxis()->GetTitle());
      href_Effic_tot->SetYTitle(HistName.Data());

      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetN();i++) {
	    double x,y;
	    gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-4) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_tot->SetMinimum(R[0]);
	href_Effic_tot->SetMaximum(R[1]);
      }
     
      href_Effic_tot->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetN() > 0) {
	  if(gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetN() == 1) {
	    gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->GetLineColor());
	    gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerStyle(SinglePointMarkerSytle);
	    gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerSize(SinglePointMarkerSize);
	  }
	  gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->Draw("PEL");
	}
      }
      
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      TString Titles_NoFakes;
      TString YTitles_NoFakes;

      Titles_NoFakes = TString(href_Effic_NoFakes->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_NoFakes->GetTitle());
      href_Effic_NoFakes->SetTitle(HistName.Data());
      
      YTitles_NoFakes = TString(href_Effic_NoFakes->GetYaxis()->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_NoFakes->GetYaxis()->GetTitle());
      href_Effic_NoFakes->SetYTitle(HistName.Data());

      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetN();i++) {
	    double x,y;
	    gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-4) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_NoFakes->SetMinimum(R[0]);
	href_Effic_NoFakes->SetMaximum(R[1]);
      }
     
      href_Effic_NoFakes->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetN() > 0) {
	  if(gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetN() == 1) {
	    gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->GetLineColor());
	    gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerStyle(SinglePointMarkerSytle);
	    gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerSize(SinglePointMarkerSize);
	  }
	  gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->Draw("PEL");
	}
      }
      
      c1->cd(4);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      TString Titles_1Fake;
      TString YTitles_1Fake;

      Titles_1Fake = TString(href_Effic_1Fake->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_1Fake->GetTitle());
      href_Effic_1Fake->SetTitle(HistName.Data());
      
      YTitles_1Fake = TString(href_Effic_1Fake->GetYaxis()->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_1Fake->GetYaxis()->GetTitle());
      href_Effic_1Fake->SetYTitle(HistName.Data());

      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetN();i++) {
	    double x,y;
	    gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_1Fake->SetMinimum(R[0]);
	href_Effic_1Fake->SetMaximum(R[1]);
      }
     
      href_Effic_1Fake->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetN() > 0) {
	  if(gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetN() == 1) {
	    gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->GetLineColor());
	    gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerStyle(SinglePointMarkerSytle);
	    gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerSize(SinglePointMarkerSize);
	  }
	  gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->Draw("PEL");
	}
      }
      
      c1->cd(5);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      TString Titles_2orMoreFakes;
      TString YTitles_2orMoreFakes;

      Titles_2orMoreFakes = TString(href_Effic_2orMoreFakes->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_2orMoreFakes->GetTitle());
      href_Effic_2orMoreFakes->SetTitle(HistName.Data());
      
      YTitles_2orMoreFakes = TString(href_Effic_2orMoreFakes->GetYaxis()->GetTitle());
      HistName = TString("Ave ") + TString(href_Effic_2orMoreFakes->GetYaxis()->GetTitle());
      href_Effic_2orMoreFakes->SetYTitle(HistName.Data());

      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetN();i++) {
	    double x,y;
	    gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_2orMoreFakes->SetMinimum(R[0]);
	href_Effic_2orMoreFakes->SetMaximum(R[1]);
      }
     
      href_Effic_2orMoreFakes->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetN() > 0) {
	  if(gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetN() == 1) {
	    gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerColor(gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->GetLineColor());
	    gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerStyle(SinglePointMarkerSytle);
	    gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->SetMarkerSize(SinglePointMarkerSize);
	  }
	  gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->Draw("PEL");
	}
      }

      c1->cd(6);
      HistName  = TString(thetaTitle);
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      c1->cd();
      PlotLogo(0.15,0.89,0.35,0.5);
      c1->Print(EPSName.Data());
    
      href_Effic_tot->SetTitle(Titles_tot.Data());
      href_Effic_tot->SetYTitle(YTitles_tot.Data());
      href_Effic_NoFakes->SetTitle(Titles_NoFakes.Data());
      href_Effic_NoFakes->SetYTitle(YTitles_NoFakes.Data());
      href_Effic_1Fake->SetTitle(Titles_1Fake.Data());
      href_Effic_1Fake->SetYTitle(YTitles_1Fake.Data());
      href_Effic_2orMoreFakes->SetTitle(Titles_2orMoreFakes.Data());
      href_Effic_2orMoreFakes->SetYTitle(YTitles_2orMoreFakes.Data());

            
      c1->Clear();
      c1->Divide(3,2);
      int counter_pad = 0;
      TString Titles[Nmax_params];
      TString YTitles[Nmax_params];
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
        counter_pad++;
        if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++; 
        c1->cd(counter_pad);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes);
      
        Titles[ipar] = TString(href[ipar]->GetTitle());
        HistName = TString("Ave ") + TString(href[ipar]->GetTitle());
        href[ipar]->SetTitle(HistName.Data());
      
        YTitles[ipar] = TString(href[ipar]->GetYaxis()->GetTitle());
        HistName = TString("Ave ") + TString(href[ipar]->GetYaxis()->GetTitle());
        href[ipar]->SetYTitle(HistName.Data());
      
        R[0] = +1.0e+20;
        R[1] = -1.0e+20;
        if(!EfficAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetN();i++) {
	      double x,y;
	      gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href[ipar]->SetMinimum(R[0]);
	  href[ipar]->SetMaximum(R[1]);
        }
      
        href[ipar]->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
	  if(gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetN() > 0) {
	    if(gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetN() == 1) {
	      gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetMarkerColor(gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->GetLineColor());
	      gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetMarkerStyle(SinglePointMarkerSytle);
	      gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->SetMarkerSize(SinglePointMarkerSize);
	    }
	    gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->Draw("PEL");
	  }
        }
      }
      c1->cd(6);
      HistName  = TString(thetaTitle);
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      
      c1->cd();
      PlotLogo(0.15,0.89,0.35,0.5);
      c1->Print(EPSName.Data());
    
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
        href[ipar]->SetTitle(Titles[ipar].Data());
        href[ipar]->SetYTitle(YTitles[ipar].Data());
      }

    }
    
  }

  if(EfficAnalysisPars[0].PlotPerformancesVsTheta) {
    for(int iP=0; iP<MyNbins_mom_redu; iP++) {
      double mom;
      if(EfficAnalysisPars[0].UseAllMomVals) mom = momArr[iP];
      else                                   mom = momArr[0] + (iP+0.5) * (momArr[momArr.size()-1] - momArr[0])/MyNbins_mom_redu;
      mom /= global->GetUnit(MonUnits);
      
      c1->Clear();
      c1->Divide(3,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      //gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_PhiAveEffic_tot_VsTheta[igeo][iP]->GetN();i++) {
	    double x,y;
	    gr_PhiAveEffic_tot_VsTheta[igeo][iP]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-4) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_tot_vs_theta->SetMinimum(R[0]);
	href_Effic_tot_vs_theta->SetMaximum(R[1]);
      }
      
      href_Effic_tot_vs_theta->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_PhiAveEffic_tot_VsTheta[igeo][iP]->GetN() > 0) gr_PhiAveEffic_tot_VsTheta[igeo][iP]->Draw("PEL");
      }
      
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      //gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_PhiAveEffic_NoFakes_VsTheta[igeo][iP]->GetN();i++) {
	    double x,y;
	    gr_PhiAveEffic_NoFakes_VsTheta[igeo][iP]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-4) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_NoFakes_vs_theta->SetMinimum(R[0]);
	href_Effic_NoFakes_vs_theta->SetMaximum(R[1]);
      }
      
      href_Effic_NoFakes_vs_theta->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_PhiAveEffic_NoFakes_VsTheta[igeo][iP]->GetN() > 0) gr_PhiAveEffic_NoFakes_VsTheta[igeo][iP]->Draw("PEL");
      }
      
      c1->cd(4);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      //gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_PhiAveEffic_1Fake_VsTheta[igeo][iP]->GetN();i++) {
	    double x,y;
	    gr_PhiAveEffic_1Fake_VsTheta[igeo][iP]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_1Fake_vs_theta->SetMinimum(R[0]);
	href_Effic_1Fake_vs_theta->SetMaximum(R[1]);
      }
      
      href_Effic_1Fake_vs_theta->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_PhiAveEffic_1Fake_VsTheta[igeo][iP]->GetN() > 0) gr_PhiAveEffic_1Fake_VsTheta[igeo][iP]->Draw("PEL");
      }
      
      c1->cd(5);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.20);
      gPad->SetBottomMargin(0.15);
      //gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes && AlsoEfficLogY);
      R[0] = +1.0e+20;
      R[1] = -1.0e+20;
      if(!EfficAnalysisPars[0].SameRange) {
	for(int igeo=0;igeo<NGeometries;igeo++) {
	  for(int i=0;i<gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][iP]->GetN();i++) {
	    double x,y;
	    gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][iP]->GetPoint(i,x,y);
	    if(y < MinEffic) continue;
	      
	    if(R[0] > y) R[0] = y;
	    if(R[1] < y) R[1] = y;
	  }
	}
	bool NoRange = false;
	if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	if(NoRange) {
	  R[0] = 1.0e-6;
	  R[1] = 1.0;
	}
	else {
	  delta = R[1] - R[0];
	  if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	  porcent = 0.1;
	  min = R[0];
	  R[0] -= porcent*delta;
	  R[1] += porcent*delta;
	  if(R[0] < 0.0) R[0] = min*min_frac;
	}
	href_Effic_2orMoreFakes_vs_theta->SetMinimum(R[0]);
	href_Effic_2orMoreFakes_vs_theta->SetMaximum(R[1]);
      }
      
      href_Effic_2orMoreFakes_vs_theta->Draw();
      for(int igeo=0;igeo<NGeometries;igeo++) {
	if(gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][iP]->GetN() > 0) gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][iP]->Draw("PEL");
      }
      
      c1->cd(6);
      sprintf(ytitle,"%.3f",mom);
      HistName  = TString("p = ") + TString(ytitle) + TString(" ") + MonUnits;
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      
      c1->cd();
      PlotLogo(0.15,0.89,0.35,0.5);
      c1->Print(EPSName.Data());
      

      c1->Clear();
      c1->Divide(3,2);
      int counter_pad = 0;
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) {
        counter_pad++;
        if(TrackerList[0]->GetTrajectory()->GetNParameters() == 4 && ipar+1 == TrackerList[0]->GetTrajectory()->GetNParameters()-1) counter_pad++;
        c1->cd(counter_pad);
        gPad->SetFillColor(10);
        gPad->SetFrameFillColor(10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.20);
        gPad->SetBottomMargin(0.15);
	//gPad->SetLogy(EfficAnalysisPars[0].UseLogYAxes);
        R[0] = +1.0e+20;
        R[1] = -1.0e+20;
        if(!EfficAnalysisPars[0].SameRange) {
	  for(int igeo=0;igeo<NGeometries;igeo++) {
	    for(int i=0;i<gr_PhiAveSigmaParVsTheta[igeo][iP][ipar]->GetN();i++) {
	      double x,y;
	      gr_PhiAveSigmaParVsTheta[igeo][iP][ipar]->GetPoint(i,x,y);
	      if(y < 0.0) continue;
	      
	      if(R[0] > y) R[0] = y;
	      if(R[1] < y) R[1] = y;
	    }
	  }
	  bool NoRange = false;
	  if(R[0] == +1.0e+20 && R[1] == -1.0e+20) NoRange = true;
	  if(NoRange) {
	    R[0] = 1.0e-6;
	    R[1] = 1.0;
	  }
	  else {
	    delta = R[1] - R[0];
	    if(TMath::Abs(delta) < 1.0e-8) delta = R[0];
	    porcent = 0.1;
	    min = R[0];
	    R[0] -= porcent*delta;
	    R[1] += porcent*delta;
	    if(R[0] < 0.0) R[0] = min*min_frac;
	  }
	  href_vs_theta[ipar]->SetMinimum(R[0]);
	  href_vs_theta[ipar]->SetMaximum(R[1]);
        }
      
        href_vs_theta[ipar]->Draw();
        for(int igeo=0;igeo<NGeometries;igeo++) {
	  if(gr_PhiAveSigmaParVsTheta[igeo][iP][ipar]->GetN() > 0) gr_PhiAveSigmaParVsTheta[igeo][iP][ipar]->Draw("PEL");
        }
      }
      c1->cd(6);
      sprintf(ytitle,"%.3f",mom);
      HistName  = TString("p = ") + TString(ytitle) + TString(" ") + MonUnits;
      latex->SetTextColor(kBlack);
      latex->DrawLatex(0.10,0.88,HistName.Data());
      leg->Draw("same");
      
      c1->cd();
      PlotLogo(0.15,0.89,0.35,0.5);
      c1->Print(EPSName.Data());
    }
  }

  c1->Print(EPSNameC.Data());

  if(SavePlots) {
    cout << endl;
    TString ROOTName = TheOutputFile + TString(".root");
    cout << "Saving track parameters resolution vs momentum to " << ROOTName.Data() << " file" << endl;
    TFile file(ROOTName.Data(),"UPDATE");

    //Saving the TGraph objects
    for(int igeo=0;igeo<NGeometries;igeo++) {
      for(int itheta=0;itheta<NthetaPoints;itheta++) {
	for(int iphi=0;iphi<NphiPoints;iphi++) {
	  if(!EfficAnalysisPars[0].PlotOnlyPhiAveraged) {
	    gr_effic_tot_vs_mom[igeo][itheta][iphi]->Write();
	    gr_effic_NoFakes_vs_mom[igeo][itheta][iphi]->Write();
	    gr_effic_1Fake_vs_mom[igeo][itheta][iphi]->Write();
	    gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi]->Write();
	    
	    for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) gr_res_ave_vs_mom[igeo][itheta][iphi][ipar]->Write();
	  }
        }
        
        if(EfficAnalysisPars[0].PlotOnlyPhiAveraged || EfficAnalysisPars[0].PlotPhiAveraged) {
	  gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta]->Write();
	  gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta]->Write();
	  gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta]->Write();
	  gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta]->Write();
	  
          for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar]->Write();
	}
	
      }
      
      if(EfficAnalysisPars[0].PlotPerformancesVsTheta) {
        for(int ip=0;ip<MyNbins_mom_redu;ip++) {
	  gr_PhiAveEffic_tot_VsTheta[igeo][ip]->Write();
	  gr_PhiAveEffic_NoFakes_VsTheta[igeo][ip]->Write();
	  gr_PhiAveEffic_1Fake_VsTheta[igeo][ip]->Write();
	  gr_PhiAveEffic_2orMoreFakes_VsTheta[igeo][ip]->Write();
	  
	  for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) gr_PhiAveSigmaParVsTheta[igeo][ip][ipar]->Write();
	}
      }
      
    }

    file.Close();
  }
  
  //Freeing the memory
  for(int igeo=0;igeo<NGeometries;igeo++) {
    for(int itheta=0;itheta<NthetaPoints;itheta++) {
      for(int iphi=0;iphi<NphiPoints;iphi++) {
	
	delete gr_effic_tot_vs_mom[igeo][itheta][iphi];
	delete gr_effic_NoFakes_vs_mom[igeo][itheta][iphi];
	delete gr_effic_1Fake_vs_mom[igeo][itheta][iphi];
	delete gr_effic_2orMoreFakes_vs_mom[igeo][itheta][iphi];
	
        for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) delete gr_res_ave_vs_mom[igeo][itheta][iphi][ipar];
      }
      
      delete gr_effic_tot_aveVsPhi_vs_mom[igeo][itheta];
      delete gr_effic_NoFakes_aveVsPhi_vs_mom[igeo][itheta];
      delete gr_effic_1Fake_aveVsPhi_vs_mom[igeo][itheta];
      delete gr_effic_2orMoreFakes_aveVsPhi_vs_mom[igeo][itheta];
      
      for(int ipar=0;ipar<TrackerList[0]->GetTrajectory()->GetNParameters();ipar++) delete gr_res_ave_aveVsPhi_vs_mom[igeo][itheta][ipar];
    }
  }

  return;
  
}
//====================================================================



