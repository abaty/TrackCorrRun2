#include "../src/Settings.h"
#include "../src/getWeights.C"
#include "getTrkCorr.C"
#include "TMath.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TAxis.h"
#include "TFile.h"
#include "TChain.h"
#include <cstring>
#include <vector>
#include <iostream>

void getClosure()
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  Settings s("trkCorrections/TrkCorrInputFile.txt");
  closureTest(s);
  return;
}

//calculates the area falling outside the acceptance (circle overlapping a rectagle)
double getArea(double eta1, double R)
{  
  if(TMath::Abs(eta1)<(2.4-R)) return TMath::Pi()*R*R;
  else
  {
    double theta = 2*TMath::ACos((2.4-TMath::Abs(eta1))/R);
    double area = R*R*(TMath::Pi()-(theta-TMath::Sin(theta))/2.0);
    return area;
  }
}

void closureTest(Settings s)
{
//Setup variables for skim
  TChain * trkCh;
  TChain * centCh;
  TChain * evtCh;
  TChain * jet;

  //track
  int nTrk;
  float trkPt[100000];
  float trkEta[100000];
  float trkPhi[100000];
  float trkStatus[100000]; //for trkStatus, -999 = fake, -99 = secondary, 1 & 2 are matched tracks
  bool highPurity[100000];
  int nVtx;

  //gen parameters
  int nParticle;
  float genPt[100000];
  float mtrkPt[100000];
  float genEta[100000];
  float genPhi[100000];
  int   mtrkQual[100000];
  float pNRec[100000];

  //other parameters
  float localTrackDensity = 0;
  
  //event parameters
  int hiBin;
  float vz;
  float pthat;
  int nref;
  float jtpt[100];
  float jtphi[100];
  float jteta[100];
  float weight = 1;
  int pcoll;

  //Setup input trees  
  //track tree     
  trkCh = new TChain("anaTrack/trackTree");
  for(int i = 0; i<s.nMC; i++)  trkCh->Add(s.MCFiles.at(i).c_str()); 
  trkCh->SetBranchAddress("nTrk",&nTrk); 
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkStatus",&trkStatus);
  if(s.doCentPU && s.nPb==0) trkCh->SetBranchAddress("nVtx",&nVtx);
  
  trkCh->SetBranchAddress("nParticle",&nParticle);
  trkCh->SetBranchAddress("pPt",&genPt);
  trkCh->SetBranchAddress("pEta",&genEta);
  trkCh->SetBranchAddress("pPhi",&genPhi);
  trkCh->SetBranchAddress("pNRec",&pNRec);
  trkCh->SetBranchAddress("mtrkPt",&mtrkPt);
//  trkCh->SetBranchAddress("mtrkQual",&mtrkQual); //for 2.76 samples
  trkCh->SetBranchAddress("mhighPurity",&mtrkQual);  //for 5.02 samples

  //centrality and vz
  centCh = new TChain("hiEvtAnalyzer/HiTree");
  for(int i = 0; i<s.nMC; i++)  centCh->Add(s.MCFiles.at(i).c_str());  
  centCh->SetBranchAddress("vz",&vz);
  if(s.doCentPU && s.nPb==2) centCh->SetBranchAddress("hiBin",&hiBin);
  trkCh->AddFriend(centCh);  
  
  //pthat and jets
  jet = new TChain(Form("%sJetAnalyzer/t",s.jetDefinition.c_str()));
  for(int i = 0; i<s.nMC; i++)  jet->Add(s.MCFiles.at(i).c_str());  
  jet->SetBranchAddress("pthat", &pthat);
  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("jtpt",&jtpt);
  jet->SetBranchAddress("jteta",&jteta);
  jet->SetBranchAddress("jtphi",&jtphi);
  trkCh->AddFriend(jet);
  
  evtCh = new TChain("skimanalysis/HltTree");
  for(int i = 0; i<s.nMC; i++)  evtCh->Add(s.MCFiles.at(i).c_str());
  evtCh->SetBranchAddress("pcollisionEventSelection",&pcoll);
  trkCh->AddFriend(evtCh);

  //event loop
  std::cout << "starting event loop" << std::endl; 
  for(int i = 0; i<10000; i++)//trkCh->GetEntries(); i++)
  {
    if(i%2000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    if(s.nPb==2)  centCh->GetEntry(i);
    else trkCh->GetEntry(i);
   
    if(TMath::Abs(vz)>s.vz_window) continue;
    if(s.nPb==2) trkCh->GetEntry(i);
    if(pcoll==0 || pthat>800) continue;
  
    //getting weight parameters
    int centPU;
    if(s.nPb==2)
    {
      weight = getWeight(s,pthat,vz,hiBin);    
      centPU = hiBin/2.0;
    }
    if(s.nPb==0) 
    {
      weight = getWeight(s,pthat,vz,nVtx);
      centPU = nVtx;
    }

    float maxJetPt = -999;
    for(int k = 0; k<nref; k++)
    {
      if(TMath::Abs(jteta[k])>2) continue;
      if(jtpt[k]>maxJetPt) maxJetPt=jtpt[k];
    }

    //track loop  
    for(int j = 0; j<nTrk; j++)
    {
      if(TMath::Abs(trkEta[j])>2.4) continue;
      if(highPurity[j]!=1) continue;
      //TODO: Calo matching here
      //other cut here as well maybe?
      //trkStauts cut here?
      if(trkPt[j]<0.5 || trkPt[j]>300) continue;

      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        if(jtpt[k]<50) break;
        if(TMath::Abs(jteta[k])>2) continue;
        float R = TMath::Power(jteta[k]-trkEta[j],2) + TMath::Power(jtphi[k]-trkPhi[j],2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }

      //fill histograms here    
  
    }
 
    //gen 
    for(int j = 0; j<nParticle; j++)
    {
      if(TMath::Abs(genEta[j])>2.4) continue;
      if(genPt[j]<0.5 || genPt[j]>300) continue;
    
      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        if(jtpt[k]<50) break;
        if(TMath::Abs(jteta[k])>2) continue;
        float R = TMath::Power(jteta[k]-genEta[j],2) + TMath::Power(jtphi[k]-genPhi[j],2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }
       
      //fill histograms
    }
  }
  //save files
}
