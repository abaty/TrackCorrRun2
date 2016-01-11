#ifndef GETTRKCORR
#define GETTRKCORR

#include "TMath.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TrkSettings.h"
#include <iostream>
#include <vector>

class TrkCorr{
  public:
    //void UpdateEventInfo(float *pt, float *eta, float *phi, int nTrk);
    double getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin, int correction=0);
    TrkCorr(std::string inputDirectory = "trkCorrections/");
    ~TrkCorr();    

  private:
    const static int nFiles = 20;
    const static int nSteps = 4;
    TH1D *       eff[nFiles][nSteps];
    TH1D *      fake[nFiles][nSteps];
    TH2D *      eff2[nFiles][nSteps];
    TH2D *     fake2[nFiles][nSteps];
    TH1D * secondary[nFiles];
    TH1D *  multiple[nFiles];

    TrkSettings * s;
};

TrkCorr::TrkCorr(std::string inputDirectory)
{
  std::cout << "Initializing tracking correction files..." << std::endl;

  s = new TrkSettings(Form("%sTrkCorrInputFile.txt",inputDirectory.c_str()));
  TFile * f[nFiles];
  for(int i = 0; i<nFiles; i++)
  {
    f[i] = TFile::Open(Form("%scorrHists_job%d.root",inputDirectory.c_str(),i),"read");
    for(int j = 0; j<nSteps; j++)
    {
      if(j!=2)
      {
        eff[i][j] = (TH1D*)f[i]->Get(Form("finalEff_type%d",j));
        eff[i][j]->SetDirectory(0);
        fake[i][j]= (TH1D*)f[i]->Get(Form("finalFake_type%d",j));
        fake[i][j]->SetDirectory(0);
      }
      else 
      {
        eff2[i][j] = (TH2D*)f[i]->Get(Form("finalEff_type%d",j));
        eff2[i][j]->SetDirectory(0);
        fake2[i][j]= (TH2D*)f[i]->Get(Form("finalFake_type%d",j));
        fake2[i][j]->SetDirectory(0);
      }
    }
    secondary[i] = (TH1D*)f[i]->Get("SecondaryRate");
    secondary[i]->SetDirectory(0);
    multiple[i] = (TH1D*)f[i]->Get("MultipleRecoRate");
    multiple[i]->SetDirectory(0);
    f[i]->Close();
  }
  std::cout << "Initialization complete." << std::endl; 
}

//correction=0 is total, 1 is eff, 2 is fake, 3 is second, 4 is mult
double TrkCorr::getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin, int correction)
{
  if(pt<0.5 || pt>=300){  std::cout << "\nPt of " << pt << " less than 500 MeV or > 300 GeV, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(eta<-2.4 || eta>2.4){  std::cout << "\nEta outside of |2.4|, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(hiBin<0 || hiBin>199){  std::cout << "\nhiBin not within 0 to 200, please place a cut to prevent this.  Returning a correction of 1" << std::endl; return 1;}
  
  float netEff = 1;
  float netFake = 1;
  float netSec = 0;
  float netMult = 0; 
 
  //calculating what file to take corrections out of 
  int coarseBin = 0;
  float cent = hiBin/2;
  if(cent>=10 && cent<20) coarseBin = coarseBin+1;
  else if(cent>=20 && cent<50) coarseBin = coarseBin+2;
  else if(cent>=50 && cent<100) coarseBin = coarseBin+3;
  if(pt>=1 && pt<3) coarseBin = coarseBin+4;
  else if(pt>=3 && pt<10) coarseBin = coarseBin+8;
  else if(pt>=10 && pt<30) coarseBin = coarseBin+12;
  else if(pt>=30 && pt<300) coarseBin = coarseBin+16;
  //end bin calculation
 
  netMult = multiple[coarseBin]->GetBinContent(multiple[coarseBin]->FindBin(pt));

  netSec  = secondary[coarseBin]->GetBinContent(secondary[coarseBin]->FindBin(pt));
  netEff *= eff[coarseBin][0]->GetBinContent(eff[coarseBin][0]->FindBin(pt));
  netEff *= eff[coarseBin][1]->GetBinContent(eff[coarseBin][1]->FindBin(cent));
  netEff *= eff2[coarseBin][2]->GetBinContent(eff2[coarseBin][2]->GetXaxis()->FindBin(eta),eff2[coarseBin][2]->GetYaxis()->FindBin(phi));
  netEff *= eff[coarseBin][3]->GetBinContent(eff[coarseBin][3]->FindBin(rmin));
  
  netFake *= fake[coarseBin][0]->GetBinContent(fake[coarseBin][0]->FindBin(pt));
  netFake *= fake[coarseBin][1]->GetBinContent(fake[coarseBin][1]->FindBin(cent));
  netFake *= fake2[coarseBin][2]->GetBinContent(fake2[coarseBin][2]->GetXaxis()->FindBin(eta),fake2[coarseBin][2]->GetYaxis()->FindBin(phi));
  netFake *= fake[coarseBin][3]->GetBinContent(fake[coarseBin][3]->FindBin(rmin));

  if(netFake<1) netFake = 1;
  if(netEff>1)  netEff = 1;

  std::cout << coarseBin << std::endl;
  std::cout << "pt: " << pt << " cent: " << cent << " eta: " << eta << " phi: " << phi   << std::endl;
  std::cout << "Efficiency: " << netEff << " (pt: " << eff[coarseBin][0]->GetBinContent(eff[coarseBin][0]->FindBin(pt)) << " cent: " <<  eff[coarseBin][1]->GetBinContent(eff[coarseBin][1]->FindBin(cent)) << " eta/phi: " << eff2[coarseBin][2]->GetBinContent(eff2[coarseBin][2]->GetXaxis()->FindBin(eta),eff2[coarseBin][2]->GetYaxis()->FindBin(phi)) << " rmin: " << eff[coarseBin][3]->GetBinContent(eff[coarseBin][3]->FindBin(rmin))<< std::endl;
  std::cout << "Fake Rate: " << (1-1./netFake) << " (pt: " << 1-1./fake[coarseBin][0]->GetBinContent(fake[coarseBin][0]->FindBin(pt)) << " cent: " << 1-1./fake[coarseBin][1]->GetBinContent(fake[coarseBin][1]->FindBin(cent)) << " eta/phi: " << 1-1./fake2[coarseBin][2]->GetBinContent(fake2[coarseBin][2]->GetXaxis()->FindBin(eta),fake2[coarseBin][2]->GetYaxis()->FindBin(phi)) << fake[coarseBin][3]->GetBinContent(fake[coarseBin][3]->FindBin(rmin)) << std::endl;
  //std::cout << "Secondary Rate: " <<  netSec << std::endl;
  //std::cout << "Multiple Reco Rate: " << netMult << "\nTotal Correction: " << (1.0-netSec)/(netEff*netFake*(1+netMult)) << std::endl;


  if(1/netEff>1000) std::cout << "problem here!" << netEff <<  " " <<pt << " " << eta << " " << phi << " " << " " << coarseBin << std::endl;

  if(correction==1) return 1/(netEff);
  else if(correction==2) return 1/(netFake);
  else if(correction==3) return 1-netSec;
  else if(correction==4) return 1/(1+netMult);
  else return 1.0/(netEff*netFake*(1+netMult));
//  else return 1.0/(netEff*netFake*(1+netMult)); 
}

TrkCorr::~TrkCorr()
{
  delete s;
}
#endif
