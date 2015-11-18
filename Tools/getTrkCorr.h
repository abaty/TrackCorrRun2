#ifndef GETTRKCORR
#define GETTRKCORR

#include "TMath.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

//written for implementing a specific PbPb mode; nFiles, nSteps, and the names of historgams might need to be changed if another training is used...

class TrkCorr{
  public:
    void UpdateEventInfo(float *pt, float *eta, float *phi, int nTrk);
    double getTrkCorr(float pt, float eta, float phi, int correction=0);
    TrkCorr();
    ~TrkCorr();    

  private:
    const static int nFiles = 5;
    const static int nSteps = 3;
    
    double getArea(double eta1, double R);

    TH2D * localDensity;
    TH1D *       eff[nFiles][nSteps];
    TH1D *      fake[nFiles][nSteps];
    TH2D *      eff2[nFiles][nSteps];
    TH2D *     fake2[nFiles][nSteps];
    TH1D * secondary[nFiles];
    TH1D *  multiple[nFiles];

    float dMapR; 
    int nEtaBin;
    int nPhiBin;
};

TrkCorr::TrkCorr()
{
  std::cout << "Initializing tracking correction files..." << std::endl;
  dMapR = 0.1;
  nEtaBin = 192;
  nPhiBin = 251;
  localDensity = new TH2D("densityMap","densityMap:eta:phi",nEtaBin,-2.4,2.4,nPhiBin,-TMath::Pi(),TMath::Pi());

  TFile * f[nFiles];
  for(int i = 0; i<nFiles; i++)
  {
    f[i] = TFile::Open(Form("trkCorrections/corrHists_job%d.root",i),"read");
    for(int j = 0; j<nSteps; j++)
    {
      if(j!=1)
      {
        eff[i][j] = (TH1D*)f[i]->Get(Form("finalEff_step%d",j));
        eff[i][j]->SetDirectory(0);
        fake[i][j]= (TH1D*)f[i]->Get(Form("finalFake_step%d",j));
        fake[i][j]->SetDirectory(0);
      }
      else 
      {
        eff2[i][j] = (TH2D*)f[i]->Get(Form("finalEff_step%d",j));
        eff2[i][j]->SetDirectory(0);
        fake2[i][j]= (TH2D*)f[i]->Get(Form("finalFake_step%d",j));
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

double TrkCorr::getArea(double eta1, double R)
{  
  if(TMath::Abs(eta1)<(2.4-R)) return TMath::Pi()*R*R;
  else
  {
    double theta = 2*TMath::ACos((2.4-TMath::Abs(eta1))/R);
    double area = R*R*(TMath::Pi()-(theta-TMath::Sin(theta))/2.0);
    return area;
  }
}

//updating the event by event properties (centrality, local density, jets, etc)
void TrkCorr::UpdateEventInfo(float pt[], float eta[], float phi[], int nTrk)
{
  localDensity->Reset();

  //Filling density histogram (tracks with >3 GeV)
  for(int j = 0; j<nTrk; j++)
  {
    if(TMath::Abs(eta[j])>2.4 || pt[j]<=3) continue;
    //loop over strip in phi (have to be careful about the -pi to pi wrap around...)
    //for case where we don't have to worry about wrap around
    if(TMath::Pi()-TMath::Abs(phi[j])>dMapR)
    {
      for(int phi_iter = localDensity->GetYaxis()->FindBin(phi[j]-dMapR); phi_iter<=localDensity->GetYaxis()->FindBin(phi[j]+dMapR); phi_iter++)
      {
        float dEtaMax = TMath::Power(dMapR*dMapR-TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi_iter))),2),0.5);
        //loop over the eta bins needed
        for(int eta_iter = localDensity->GetXaxis()->FindBin(eta[j]-dEtaMax); eta_iter<=localDensity->GetXaxis()->FindBin(eta[j]+dEtaMax); eta_iter++)
        {
          if(TMath::Power(eta[j]-localDensity->GetXaxis()->GetBinCenter(eta_iter),2)+TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi_iter))),2)<dMapR*dMapR ) localDensity->SetBinContent(eta_iter,phi_iter,localDensity->GetBinContent(eta_iter,phi_iter)+1); 
        }
      }
    }
    else
    //for case with -pi and pi wrap around 
    for(int phi_iter = 1; phi_iter<=localDensity->GetYaxis()->FindBin(phi[j]+dMapR-(phi[j]>0?2*TMath::Pi():0)); phi_iter++) 
    {
      //loop over the eta bins needed
      float dEtaMax = TMath::Power(dMapR*dMapR-TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi_iter))),2),0.5);
      for(int eta_iter = localDensity->GetXaxis()->FindBin(eta[j]-dEtaMax); eta_iter<=localDensity->GetXaxis()->FindBin(eta[j]+dEtaMax); eta_iter++)
      {
        if(TMath::Power(eta[j]-localDensity->GetXaxis()->GetBinCenter(eta_iter),2)+TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi_iter))),2)<dMapR*dMapR ) localDensity->SetBinContent(eta_iter,phi_iter,localDensity->GetBinContent(eta_iter,phi_iter)+1); 
      }
    }
    for(int phi_iter = localDensity->GetYaxis()->FindBin(phi[j]-dMapR+(phi[j]<0?2*TMath::Pi():0)); phi_iter<=nPhiBin; phi_iter++) 
    {  
      //loop over the eta bins needed
      float dEtaMax = TMath::Power(dMapR*dMapR-TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi_iter))),2),0.5);
      for(int eta_iter = localDensity->GetXaxis()->FindBin(eta[j]-dEtaMax); eta_iter<=localDensity->GetXaxis()->FindBin(eta[j]+dEtaMax); eta_iter++)
      {
        if(TMath::Power(eta[j]-localDensity->GetXaxis()->GetBinCenter(eta_iter),2)+TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi_iter))),2)<dMapR*dMapR ) localDensity->SetBinContent(eta_iter,phi_iter,localDensity->GetBinContent(eta_iter,phi_iter)+1); 
      }
    }
  }
}


//correction=0 is total, 1 is eff, 2 is fake, 3 is second, 4 is mult
double TrkCorr::getTrkCorr(float pt, float eta, float phi, int correction)
{
  if(pt<0.5 || pt>=300){  std::cout << "\nPt of " << pt << " less than 500 MeV or > 300 GeV, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(eta<-2.4 || eta>2.4){  std::cout << "\nEta outside of |2.4|, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  //if(hiBin<0 || hiBin>199){  std::cout << "\nhiBin not within 0 to 200, please place a cut to prevent this.  Returning a correction of 1" << std::endl; return 1;}
  
  float netEff = 1;
  float netFake = 1;
  float netSec = 0;
  float netMult = 0; 

  float density = localDensity->GetBinContent(localDensity->GetXaxis()->FindBin(eta), localDensity->GetYaxis()->FindBin(phi))/getArea(eta,dMapR);
 
  //calculating what file to take corrections out of 
  int coarseBin = 0;
  //int cent = hiBin/2;
  //if(cent>=10 && cent<20) coarseBin = coarseBin+1;
  //else if(cent>=20 && cent<50) coarseBin = coarseBin+2;
  //else if(cent>=50 && cent<100) coarseBin = coarseBin+3;
  if(pt>=1 && pt<3) coarseBin = coarseBin+1;
  else if(pt>=3 && pt<10) coarseBin = coarseBin+2;
  else if(pt>=10 && pt<30) coarseBin = coarseBin+3;
  else if(pt>=30 && pt<300) coarseBin = coarseBin+4;
  //end bin calculation
 
  netMult = multiple[coarseBin]->GetBinContent(multiple[coarseBin]->FindBin(pt));

  netSec  = secondary[coarseBin]->GetBinContent(secondary[coarseBin]->FindBin(pt));
  netEff *= eff[coarseBin][0]->GetBinContent(eff[coarseBin][0]->FindBin(pt));
  //netEff *= eff[coarseBin][1]->GetBinContent(eff[coarseBin][1]->FindBin(cent));
  netEff *= eff2[coarseBin][1]->GetBinContent(eff2[coarseBin][1]->GetXaxis()->FindBin(eta),eff2[coarseBin][1]->GetYaxis()->FindBin(phi));
  netEff *= eff[coarseBin][2]->GetBinContent(eff[coarseBin][2]->FindBin(density));
  
  netFake *= fake[coarseBin][0]->GetBinContent(fake[coarseBin][0]->FindBin(pt));
  //netFake *= fake[coarseBin][1]->GetBinContent(fake[coarseBin][1]->FindBin(cent));
  netFake *= fake2[coarseBin][1]->GetBinContent(fake2[coarseBin][1]->GetXaxis()->FindBin(eta),fake2[coarseBin][1]->GetYaxis()->FindBin(phi));
  netFake *= fake[coarseBin][2]->GetBinContent(fake[coarseBin][2]->FindBin(density));

  if(netFake<1) netFake = 1;
  if(netEff>1)  netEff = 1;
/*
  std::cout << "pt: " << pt << " cent: " << cent << " eta: " << eta << " phi: " << phi << " density: " << density << std::endl;
  std::cout << "Efficiency: " << netEff << " (pt: " << eff[coarseBin][0]->GetBinContent(eff[coarseBin][0]->FindBin(pt)) << " cent: " <<  eff[coarseBin][1]->GetBinContent(eff[coarseBin][1]->FindBin(cent)) << " eta/phi: " << eff2[coarseBin][2]->GetBinContent(eff2[coarseBin][2]->GetXaxis()->FindBin(eta),eff2[coarseBin][2]->GetYaxis()->FindBin(phi)) << " density: " << eff[coarseBin][3]->GetBinContent(eff[coarseBin][3]->FindBin(density)) <<")"<< std::endl;
  std::cout << "Fake Rate: " << (1-1./netFake) << " (pt: " << 1-1./fake[coarseBin][0]->GetBinContent(fake[coarseBin][0]->FindBin(pt)) << " cent: " << 1-1./fake[coarseBin][1]->GetBinContent(fake[coarseBin][1]->FindBin(cent)) << " eta/phi: " << 1-1./fake2[coarseBin][2]->GetBinContent(fake2[coarseBin][2]->GetXaxis()->FindBin(eta),fake2[coarseBin][2]->GetYaxis()->FindBin(phi)) << " density: " << fake[coarseBin][3]->GetBinContent(fake[coarseBin][3]->FindBin(density)) <<")"<< std::endl;
//  std::cout << "Secondary Rate: " <<  netSec << std::endl;
  std::cout << "Multiple Reco Rate: " << netMult << "\nTotal Correction: " << (1.0-netSec)/(netEff*netFake*(1+netMult)) << std::endl;
*/

  if(1/netEff>1000) std::cout << "problem here!" << netEff <<  " " <<pt << " " << eta << " " << phi << " " << density << " " << coarseBin << std::endl;

  if(correction==1) return 1/(netEff);
  else if(correction==2) return 1/(netFake);
  else if(correction==3) return 1-netSec;
  else if(correction==4) return 1/(1+netMult);
  else return (1.0-netSec)/(netEff*netFake*(1+netMult));
//  else return 1.0/(netEff*netFake*(1+netMult)); 
}

TrkCorr::~TrkCorr()
{
}
#endif
