#define Subtract_cxx
#include "Subtract.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TF1.h>


#define DRAW_TWOPHOTON 1
//#define DRAW_CASCADE 1



// Linear background function                                                                                          
double background(double *x, double *par)
{
  return par[0]*x[0] + par[1];
}


// Gaussian Peak function                                                                                              
double gaussianPeak(double *x, double *par)
{
  return (par[0] * exp( -pow(x[0]-par[1],2) / (2*pow(par[2],2)) ));
}


// fit function                                                                                                       
double fitFunction(double *x, double *par)
{
  return background(x,&par[3]) + gaussianPeak(x,par);
}


const int noAngleGr = 2;
const int noTransitions = 2;

void Subtract::FillHisto()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   const int noRegions = 2;
   const int binning = 260;
   const int maxEnergy = 2000;
   double errors[binning][2][2];
   const int regions[noRegions+1] = {0,150,350};
   const int energyLow = 100;
   const double timeUp = 0.3;
   const double timeDown = -0.3;
   const int randomGate1[2] = {-20,-5};
   const int randomGate2[2] = {5,20};
   const double scaleFactor = (timeUp-timeDown)/double(fabs(randomGate1[0]-randomGate1[1]) + fabs(randomGate2[0]-randomGate2[1])); 
   const int timeEnergyGate[2] = {650,674};

   TH1D* timeSpectrum[noAngleGr][noRegions]; 
   TH1D* sumEnergyVeto[noAngleGr][noRegions]; 
   TH1D* sumEnergyVetoRandom[noAngleGr][noRegions]; 
   TH1D* sumEnergyVetoSubtr[noAngleGr][noRegions]; 
   TH1D* casTimeSpectrum[noAngleGr][noTransitions]; 
   TH1D* casVeto[noAngleGr][noTransitions]; 
   TH1D* casVetoRandom[noAngleGr][noTransitions]; 
   TH1D* casVetoSubtr[noAngleGr][noTransitions]; 


   for(int ang=0;ang<noAngleGr;ang++)
   {
      for(int reg=0;reg<noRegions;reg++)
      {
         timeSpectrum[ang][reg] = new TH1D(Form("timeSpectrum, angle %d deg, region %d",(ang+1)*72,reg),
					   Form("timeSpectrum, angle %d deg, region %d",(ang+1)*72,reg),200,-21,21);
         sumEnergyVeto[ang][reg] = new TH1D(Form("sumEnergyVeto, angle %d deg, region %d",(ang+1)*72,reg),
					   Form("sumEnergyVeto, angle %d deg, region %d",(ang+1)*72,reg),binning,0,maxEnergy);
         sumEnergyVetoRandom[ang][reg] = new TH1D(Form("sumEnergyVetoRandom, angle %d deg, region %d",(ang+1)*72,reg),
					   Form("sumEnergyVetoRandom, angle %d deg, region %d",(ang+1)*72,reg),binning,0,maxEnergy);
         sumEnergyVetoSubtr[ang][reg] = new TH1D(Form("sumEnergyVetoSubtracted, angle %d deg, region %d",(ang+1)*72,reg),
					   Form("sumEnergyVetoSubtracted, angle %d deg, region %d",(ang+1)*72,reg),binning,0,maxEnergy);
      }
      for(int trans=0;trans<noTransitions;trans++)
      {
         casTimeSpectrum[ang][trans] = new TH1D(Form("casTimeSpectrum, angle %d deg, transition %d",(ang+1)*72,trans),
					   Form("casTimeSpectrum, angle %d deg, transition %d",(ang+1)*72,trans),200,-21,21);
         casVeto[ang][trans] = new TH1D(Form("casVeto, angle %d deg, transition %d",(ang+1)*72,trans),
					   Form("casVeto, angle %d deg, transition %d",(ang+1)*72,trans),binning,0,maxEnergy);
         casVetoRandom[ang][trans] = new TH1D(Form("casVetoRandom, angle %d deg, transition %d",(ang+1)*72,trans),
					   Form("casVetoRandom, angle %d deg, transition %d",(ang+1)*72,trans),binning,0,maxEnergy);
         casVetoSubtr[ang][trans] = new TH1D(Form("casVetoSubtracted, angle %d deg, transition %d",(ang+1)*72,trans),
					   Form("casVetoSubtracted, angle %d deg, transition %d",(ang+1)*72,trans),binning,0,maxEnergy);
      }
   }


   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);  

      int angleGroup = -1;
      if(fChain->GetTreeNumber()<3)  // runs 206, 207, 208 have a problem with pairNo==4
      {
	if( pairNo==0 || pairNo==3 || pairNo==7 || pairNo==9)   
            angleGroup = 0;
      }
      else
      {
	if( pairNo==0 || pairNo==3 || pairNo==7 || pairNo==9 || pairNo==4)   
            angleGroup = 0;
      }
      if(pairNo==1 || pairNo==2 || pairNo==5 || pairNo==6 || pairNo==8)
         angleGroup = 1;
      if(angleGroup==-1 && pairNo!=4)
      {
         std::cout << "Error: Something wrong with pairNo." << std::endl;
	 exit(-1);
      }
      if(angleGroup==-1 && pairNo==4)
         continue;

      // two-Photons histograms
      for(int reg=0;reg<noRegions;reg++)
      {
	if( e1>energyLow && e2>energyLow && veto==0 && angleGroup!=-1 && fabs(e1-e2)>regions[reg] && fabs(e1-e2)<regions[reg+1] )
	  //	    && checkSequVeto(angleGroup,e1) && checkSequVeto(angleGroup,e2)  )
	  //	  && checkCascade(e1) && checkCascade(e2) )
         {
	    if( tdiff>timeDown && tdiff<timeUp )
	       sumEnergyVeto[angleGroup][reg]->Fill(e1+e2);
            if( (tdiff>randomGate1[0] && tdiff<randomGate1[1]) || (tdiff>randomGate2[0] && tdiff<randomGate2[1]) )  
	       sumEnergyVetoRandom[angleGroup][reg]->Fill(e1+e2);
	    if( e1+e2>timeEnergyGate[0] && e1+e2<timeEnergyGate[1] )
	       timeSpectrum[angleGroup][reg]->Fill(tdiff);
         }
      }

      // cascade histograms (283, 379)
      for(int trans=0;trans<noTransitions;trans++)
      {
	 if( e1>energyLow && e2>energyLow && veto==0 && angleGroup!=-1 && e1+e2>timeEnergyGate[0] && e1+e2<timeEnergyGate[1] )
	 {
	    if( tdiff>timeDown && tdiff<timeUp )
	    {
	       casVeto[angleGroup][trans]->Fill(e1);
	       casVeto[angleGroup][trans]->Fill(e2);
	    }
            if( (tdiff>randomGate1[0] && tdiff<randomGate1[1]) || (tdiff>randomGate2[0] && tdiff<randomGate2[1]) )  
	    {
	       casVetoRandom[angleGroup][trans]->Fill(e1);
	       casVetoRandom[angleGroup][trans]->Fill(e2);
	    }
	 }
      }


      if(jentry%100000==0 && jentry!=0)
      {
         std::cout << "Processed Events:  " << jentry << "\n\033[1A";
      }
      if(jentry==nentries-1)
	std::cout << std::endl;      
   }

   TCanvas* c[noAngleGr];
   for(int ang=0;ang<noAngleGr;ang++)
   {
      c[ang] = new TCanvas(Form("c%d",ang+1),Form("c%d",ang+1),700,3000);
      c[ang]->Divide(1,4);
   }


#ifdef DRAW_CASCADE
   for(int ang=0;ang<noAngleGr;ang++)
   {
      c[ang]->cd(1);
      timeSpectrum[ang][0]->SetLineColor(kBlue);
      timeSpectrum[ang][0]->Draw();
      timeSpectrum[ang][1]->SetLineColor(kGreen);
      timeSpectrum[ang][1]->Draw("same");
      c[ang]->cd(2);
      casVeto[ang][0]->DrawCopy();
      casVetoRandom[ang][0]->Scale(scaleFactor);
      casVetoRandom[ang][1]->Scale(scaleFactor);
      casVetoRandom[ang][0]->SetLineColor(kRed);
      casVetoRandom[ang][0]->DrawCopy("same");
      for(int j=0;j<binning;j++)
      {
         errors[j][ang][0] = casVeto[ang][0]->GetBinError(j);
         errors[j][ang][1] = casVeto[ang][1]->GetBinError(j);
      }
      casVeto[ang][0]->Add(casVetoRandom[ang][0],-1);
      casVeto[ang][1]->Add(casVetoRandom[ang][1],-1);

      casVetoSubtr[ang][0] = casVeto[ang][0];
      casVetoSubtr[ang][1] = casVeto[ang][1];
      double err1[binning] = {0};
      double err2[binning] = {0};
      for(int j=0;j<binning;j++)
      {
         err1[j] = errors[j][ang][0];
	 err2[j] = errors[j][ang][1];
      }
      casVetoSubtr[ang][0]->SetError(err1);
      casVetoSubtr[ang][1]->SetError(err2);
      casVetoSubtr[ang][0]->SetMarkerSize(0.5);
      casVetoSubtr[ang][1]->SetMarkerSize(0.5);
      casVetoSubtr[ang][0]->SetMarkerStyle(20);
      casVetoSubtr[ang][1]->SetMarkerStyle(20);
      casVetoSubtr[ang][0]->SetLineColor(kGreen);
      casVetoSubtr[ang][1]->SetLineColor(kGreen);      

      c[ang]->cd(3);
      gStyle->SetErrorX(0.000001);
      TAxis *axisX1 = casVetoSubtr[ang][0]->GetXaxis();
      TAxis *axisY1 = casVetoSubtr[ang][0]->GetYaxis();
      TAxis *axisX2 = casVetoSubtr[ang][1]->GetXaxis();
      TAxis *axisY2 = casVetoSubtr[ang][1]->GetYaxis();
      axisX1->SetLabelSize(0.07);
      axisY1->SetLabelSize(0.07);
      axisX2->SetLabelSize(0.07);
      axisY2->SetLabelSize(0.07);
      axisX1->SetRange(0,binning/3);
      axisX2->SetRange(0,binning/3);
      casVetoSubtr[ang][0]->Draw();

      c[ang]->cd(4);
      gStyle->SetErrorX(0.000001);
      casVetoSubtr[ang][1]->Draw();
   }
#endif


#ifdef DRAW_TWOPHOTON
   for(int ang=0;ang<noAngleGr;ang++)
   {
      c[ang]->cd(1);
      timeSpectrum[ang][0]->SetLineColor(kBlue);
      timeSpectrum[ang][0]->Draw();
      timeSpectrum[ang][1]->SetLineColor(kGreen);
      timeSpectrum[ang][1]->Draw("same");
      c[ang]->cd(2);
      sumEnergyVeto[ang][0]->DrawCopy();
      sumEnergyVetoRandom[ang][0]->Scale(scaleFactor);
      sumEnergyVetoRandom[ang][1]->Scale(scaleFactor);
      sumEnergyVetoRandom[ang][0]->SetLineColor(kRed);
      sumEnergyVetoRandom[ang][0]->DrawCopy("same");

      for(int j=0;j<binning;j++)
      {
         errors[j][ang][0] = sumEnergyVeto[ang][0]->GetBinError(j);
         errors[j][ang][1] = sumEnergyVeto[ang][1]->GetBinError(j);
      }
      sumEnergyVeto[ang][0]->Add(sumEnergyVetoRandom[ang][0],-1);
      sumEnergyVeto[ang][1]->Add(sumEnergyVetoRandom[ang][1],-1);

      sumEnergyVetoSubtr[ang][0] = sumEnergyVeto[ang][0];
      sumEnergyVetoSubtr[ang][1] = sumEnergyVeto[ang][1];
      double err1[binning] = {0};
      double err2[binning] = {0};
      for(int j=0;j<binning;j++)
      {
         err1[j] = errors[j][ang][0];
	 err2[j] = errors[j][ang][1];
      }
      sumEnergyVetoSubtr[ang][0]->SetError(err1);
      sumEnergyVetoSubtr[ang][1]->SetError(err2);
      sumEnergyVetoSubtr[ang][0]->SetMarkerSize(0.5);
      sumEnergyVetoSubtr[ang][1]->SetMarkerSize(0.5);
      sumEnergyVetoSubtr[ang][0]->SetMarkerStyle(20);
      sumEnergyVetoSubtr[ang][1]->SetMarkerStyle(20);
      sumEnergyVetoSubtr[ang][0]->SetLineColor(kGreen);
      sumEnergyVetoSubtr[ang][1]->SetLineColor(kGreen);      
      c[ang]->cd(3);
      gStyle->SetErrorX(0.000001);
      TAxis *axisX1 = sumEnergyVetoSubtr[ang][0]->GetXaxis();
      TAxis *axisY1 = sumEnergyVetoSubtr[ang][0]->GetYaxis();
      TAxis *axisX2 = sumEnergyVetoSubtr[ang][1]->GetXaxis();
      TAxis *axisY2 = sumEnergyVetoSubtr[ang][1]->GetYaxis();
      axisX1->SetLabelSize(0.07);
      axisY1->SetLabelSize(0.07);
      axisX2->SetLabelSize(0.07);
      axisY2->SetLabelSize(0.07);
      axisX1->SetRange(binning/4+5,binning/2);
      axisX2->SetRange(binning/4+5,binning/2);
      sumEnergyVetoSubtr[ang][0]->Draw();
      c[ang]->cd(4);
      gStyle->SetErrorX(0.000001);
      sumEnergyVetoSubtr[ang][1]->Draw();
   }
#endif
}



Bool_t Subtract::checkSequVeto(const int angleGr, const double e)
{
   const int width = 20;
   const int energyVeto[noAngleGr][4] = {{125-width,125+width,478-width,478+width},{172-width,172+width,478-width,478+width}};
   if( (e>energyVeto[angleGr][0] && e<energyVeto[angleGr][1]) || (e>energyVeto[angleGr][2] && e<energyVeto[angleGr][3]) )
   {
      return false;
   }
   return true;
}


Bool_t Subtract::checkCascade(const double e)
{
   const int width = 1;
   const int cascadeEnergy[4] = {283-width,283+width,379-width,379+width};
   if( (e>cascadeEnergy[0] && e<cascadeEnergy[1]) || (e>cascadeEnergy[2] && e<cascadeEnergy[3]) )
   {
      return false;
   }
   return true;
}







   /*
   std::ofstream myfile;
   myfile.open ("timeSpectrum.dat");
   for(int i=0;i<binning;i++)
   {
      myfile << i/4.0 << "   " << timeSpectrum->GetBinContent(i) << std::endl;
   }
   myfile.close();
   */

   /*
   TF1* fitFcn = new TF1("fitFcn",fitFunction,startFit,stopFit,5);
   fitFcn->SetParameters(50,662,20,-0.3,250);
   fitFcn->FixParameter(2,25/2.35);
   sumEnergyVetoSubtracted->Fit("fitFcn","RQ");
   std::cout << "Chi-square:         "<< fitFcn->GetChisquare() << std::endl;
   std::cout << "Degrees of freedom: " << fitFcn->GetNDF() << std::endl;
   std::cout << "Integral:           "
             << fitFcn->GetParameter(0)*fitFcn->GetParameter(2)*sqrt(2*3.14)*((double)binning)/maxEnergy
             << std::endl;*/

  

/*
   const int startFit = 600;                                                                                           
   const int stopFit = 750;
*/
