#include "Settings.h"
#include "getWeights.C"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include <cstring>
#include <vector>
#include <iostream>

void makeSkim(Settings s, int job, const char * skimType)
{
//figuring out which subset of data to skim based on the job number
  int nSkip = 0;
  double ptMin = 0.5;
  double ptMax = 1;
  double centPUMin = 0;
  double centPUMax = 1;

  int jobTemp = job;
  for(int i = 0; i<s.nPtBinCoarse; i++)
  {
    if(jobTemp>=s.nCentPUBinCoarse.at(i)) jobTemp = jobTemp-s.nCentPUBinCoarse.at(i);
    else
    { 
      nSkip = s.eventSkip.at(i).at(jobTemp);
      ptMin = s.ptBinCoarse.at(i);
      ptMax = s.ptBinCoarse.at(i+1);
      centPUMin = s.centPUBinCoarse.at(i).at(jobTemp);
      centPUMax = s.centPUBinCoarse.at(i).at(jobTemp+1);
      break;
    }
  }

  std::cout << "\nJob number: " << job << "\nCorresponds to the following parameters\nnSkip: " << nSkip
  << "\nptMin: " << ptMin << "\nptMax: " << ptMax << "\ncentMin: " << centPUMin << "\ncentMax " << centPUMax << std::endl;

  std::cout << "\nCreating initial skim of important information" << std::endl;

//Setup variables for skim
  TChain * trkCh;
  TChain * centCh;
  TChain * jet;
  int nTrk;
  float trkPt[25000];
  float trkEta[25000];
  float trkPhi;
  float highPurity;
  float genPt;
  float genEta;
  float genPhi;
  float mtrkPt;
  int nVtx;
  float localtrackDensity;
  float eff = 1;
  float fake = 0;
  
  int hiBin;
  float vz;
  float pthat;
  float weight = 1;

  //Setup input trees  
  //track tree     
  trkCh = new TChain("anaTrack/trackTree");
  for(int i = 0; i<s.nMC; i++)  trkCh->Add(s.MCFiles.at(i).c_str()); 
  trkCh->SetBranchAddress("nTrk",&nTrk); 
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  //trkCh->SetBranchAddress("trkPhi",&trkPhi);
  //trkCh->SetBranchAddress("highPurity",&highPurity);
  //trkCh->SetBranchAddress("nVtx",&nVtx);
  
  //centrality and vz
  centCh = new TChain("hiEvtAnalyzer/HiTree");
  for(int i = 0; i<s.nMC; i++)  centCh->Add(s.MCFiles.at(i).c_str());  
  centCh->SetBranchAddress("vz",&vz);
  if(s.doCentPU && s.nPb==2) centCh->SetBranchAddress("hiBin",&hiBin);
  trkCh->AddFriend(centCh);  
  
  //pthat
  jet = new TChain(Form("%sJetAnalyzer/t",s.jetDefinition.c_str()));
  for(int i = 0; i<s.nMC; i++)  jet->Add(s.MCFiles.at(i).c_str());  
  jet->SetBranchAddress("pthat", &pthat);
  trkCh->AddFriend(jet);
  
  //Setup output Ntuples
  std::string trackVars;
  std::string particleVars;
  if(strcmp(skimType,"Eff")==0)
  {
    //particleVars="genPt:matchedpt:genEta:genPhi:genDensity:centPU:eff:weight";
    //trackVars=   "trlPt:trkEta:trkPhi:trkDensity:trackselect:trackstatus:centPU:eff:weight";
    trackVars = "trkPt:trkEta";
  }
  else if(strcmp(skimType,"Fake")==0)
  {
    particleVars="pt:matchedpt:eta:phi:density:trackselect:centPU:fake:weight";
    trackVars=   "pt:eta:phi:density:trackselect:trackstatus:centPU:trkfake:fake:weight";
  }

  TFile * skimOut = TFile::Open("trackSkim_step0.root","recreate");
  TNtuple * gen  = new TNtuple("Gen","",particleVars.data()); 
  TNtuple * reco = new TNtuple("Reco","",trackVars.data());

  std::cout << "starting skim loop" << std::endl;
  //Actual skimming
  int processed = 0;
  for(int i = 0; i<1000;i++)//trkCh->GetEntries(); i++)
  {
    std::cout << i << std::endl;
    if(s.nPb==2)  centCh->GetEntry(i);
    else trkCh->GetEntry(i);
  
    if((s.nPb==2) && ((hiBin/2 < centPUMin) || (hiBin/2 >= centPUMax))) continue;
    else if((s.nPb==0) && ((nVtx < centPUMin) || (nVtx >= centPUMax))) continue;
    else if(TMath::Abs(vz)>s.vz_window) continue;
    else if(processed%nSkip !=0)
    {
      processed++;
      //std::cout << "skipping event: " << i << std::endl;
      continue;
    }
    if(s.nPb==2) trkCh->GetEntry(i);
  
    //make sure to fix this
    //if(s.nPb==2) weight = getWeight(s,pthat,vz,hiBin);    
    //if(s.nPb==0) weight = getWeight(s,pthat,vz,nVtx);
  
    float entry[] = {trkPt[0],trkEta[0]};
    reco->Fill(entry); 
    processed++;
  }
  skimOut->Write();
  skimOut->Close();   
}
