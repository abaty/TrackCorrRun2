#ifndef GETTRKCORR
#define GETTRKCORR

#include "TMath.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

class TrkCorr{
  public:
    void UpdateEventInfo();
    double getTrkCorr();
    const int nFiles = 20;
    const int nSteps = 20;

  private:
    double getArea();

    TH2D * localDensity;
    TH1D *       eff[nFiles][nSteps];
    TH1D *      fake[nFiles][nSteps];
    TH2D *      eff2[nFiles][nSteps];
    TH2D *     fake2[nFiles][nSteps];
    TH1D * secondary[nFiles];
    TH1D *  multiple[nFiles];

    TTree * trkTree_;
    bool hasTree;
    int nTrk_;
    float pt_[100000];
    float eta_[100000];
    float phi_[100000];
    bool  highPurity_[100000];

    float dMapR; 
    int nEtaBin;
    int nPhiBin;
};

void TrkCorr::TrkCorr()
{
  std::cout << "Initializing tracking correction files..." << std::endl;
  bool hasTree = false;

  float dMapR = 0.1;
  int nEtaBin = 192;
  int nPhiBin = 251;
  localDensity = new TH2D("densityMap","densityMap:eta:phi",nEtaBin,-2.4,2.4,nPhiBin,-TMath::Pi(),TMath::Pi());

  TFile * f[nFiles];
  for(int i = 0; i<nFiles; i++)
  {
    f[i] = TFile::Open(Form("trkCorrections/corrHists_job%d",i),"read");
    for(int j = 0; j<nSteps; j++)
    {
      eff[i][j] = (TH1D*)f[i]->Get(Form(,j));
      eff[i][j]->SetDirectory(0);
      fake[i][j]= (TH1D*)f[i]->Get(Form(,j));
      eff[i][j]->SetDirectory(0);
    }
    secondary[i] = (TH1D*)f[i]->Get();
    secondary[i]->SetDirectory(0);
    multiple[i] = (TH1D*)f[i]->Get();
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

//optional way to calling with a Tree instead of w/ the arrays themselves (probably not preferred)
void TrkCorr::UpdateEventInfo(TTree* trkTree, int evtNumber, bool resetTree)
{
  if(resetTree || hasTree==0)
  {
    trkTree = trkTree;
    trkTree_->SetBranchAddress("nTrk",&nTrk_);
    trkTree_->SetBranchAddress("trkPt",&pt_);
    trkTree_->SetBranchAddress("trkEta",&eta_);
    trkTree_->SetBranchAddress("trkPhi",&phi_);
    trkTree_->SetBranchAddress("highPurity",&highPurity_);
    hasTree = 1;
  }
  _trkTree->GetEntry(evtNumber);
  TrkCorr::UpdateEventInfo(&pt_,&eta_,&phi_,&highPurity_,nTrk_);
}

//updating the event by event properties (centrality, local density, jets, etc)
void TrkCorr::UpdateEventInfo(float *pt, float *eta, float *phi, bool *highPurity, int nTrk)
{
  localDensity->Reset();

  //Filling density histogram (tracks with >3 GeV and highPurity)
  for(int j = 0; j<nTrk; j++)
  {
    if(TMath::Abs(eta[j])>2.4 || pt[j]<3 || highPurity[j]!=1) continue;
    //loop over strip in phi (have to be careful about the -pi to pi wrap around...)
    //for case where we don't have to worry about wrap around
    if(TMath::Pi()-TMath::Abs(phi[j])>dMapR)
    {
      for(int phi = localDensity->GetYaxis()->FindBin(phi[j]-dMapR); phi<=localDensity->GetYaxis()->FindBin(phi[j]+dMapR); phi++)
      {
        float dEtaMax = TMath::Power(dMapR*dMapR-TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi))),2),0.5);
        //loop over the eta bins needed
        for(int eta = localDensity->GetXaxis()->FindBin(eta[j]-dEtaMax); eta<=localDensity->GetXaxis()->FindBin(eta[j]+dEtaMax); eta++)
        {
          if(TMath::Power(eta[j]-localDensity->GetXaxis()->GetBinCenter(eta),2)+TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi))),2)<dMapR*dMapR ) localDensity->SetBinContent(eta,phi,localDensity->GetBinContent(eta,phi)+1); 
        }
      }
    }
    else
    //for case with -pi and pi wrap around 
    for(int phi = 1; phi<=localDensity->GetYaxis()->FindBin(phi[j]+dMapR-(phi[j]>0?2*TMath::Pi():0)); phi++) 
    {
      //loop over the eta bins needed
      float dEtaMax = TMath::Power(dMapR*dMapR-TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi))),2),0.5);
      for(int eta = localDensity->GetXaxis()->FindBin(eta[j]-dEtaMax); eta<=localDensity->GetXaxis()->FindBin(eta[j]+dEtaMax); eta++)
      {
        if(TMath::Power(eta[j]-localDensity->GetXaxis()->GetBinCenter(eta),2)+TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi))),2)<dMapR*dMapR ) localDensity->SetBinContent(eta,phi,localDensity->GetBinContent(eta,phi)+1); 
      }
    }
    for(int phi = localDensity->GetYaxis()->FindBin(phi[j]-dMapR+(phi[j]<0?2*TMath::Pi():0)); phi<=nPhiBin; phi++) 
    {  
      //loop over the eta bins needed
      float dEtaMax = TMath::Power(dMapR*dMapR-TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi))),2),0.5);
      for(int eta = localDensity->GetXaxis()->FindBin(eta[j]-dEtaMax); eta<=localDensity->GetXaxis()->FindBin(eta[j]+dEtaMax); eta++)
      {
        if(TMath::Power(eta[j]-localDensity->GetXaxis()->GetBinCenter(eta),2)+TMath::Power(TMath::ACos(TMath::Cos(phi[j]-localDensity->GetYaxis()->GetBinCenter(phi))),2)<dMapR*dMapR ) localDensity->SetBinContent(eta,phi,localDensity->GetBinContent(eta,phi)+1); 
      }
    }
  }
}

double TrkCorr::getTrkCorr(float pt, float eta, float phi, int hiBinPU)
{

}

void TrkCorr::~TrkCorr()
{
  delete localDensity;
}
#endif
