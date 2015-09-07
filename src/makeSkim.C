#include "Settings.h"
#include "getWeights.C"
#include "TMath.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TAxis.h"
#include "TFile.h"
#include "TChain.h"
#include <cstring>
#include <vector>
#include <iostream>

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

  //track
  int nTrk;
  float trkPt[100000];
  float trkEta[100000];
  float trkPhi[100000];
  float trkFake[100000];
  bool highPurity[100000];
  int nVtx;

  //gen parameters
  int nParticle;
  float genPt[100000];
  float genEta[100000];
  float genPhi[100000];

  //other parameters
  float localTrackDensity = 0;
  
  //event parameters
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
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  if(strcmp(skimType,"Fake")==0) trkCh->SetBranchAddress("trkFake",&trkFake);
  if(s.doCentPU && s.nPb==0) trkCh->SetBranchAddress("nVtx",&nVtx);
  
  trkCh->SetBranchAddress("nParticle",&nParticle);
  trkCh->SetBranchAddress("pPt",&genPt);
  trkCh->SetBranchAddress("pEta",&genEta);
  trkCh->SetBranchAddress("pPhi",&genPhi);
  
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
  trkCh->AddFriend(jet);
  
  //Setup output Ntuples
  std::string trackVars;
  std::string particleVars;
  if(strcmp(skimType,"Eff")==0)
  {
    particleVars="genPt:genEta:genPhi:genDensity:weight";
    trackVars=   "trkPt:trkEta:trkPhi:trkDensity:weight";
  }
  else if(strcmp(skimType,"Fake")==0)
  {
    particleVars="genPt:genEta:genPhi:genDensity:weight";
    trackVars=   "trkPt:trkEta:trkPhi:trkDensity:trkFake:weight";
  }

  TFile * skimOut = TFile::Open("trackSkim_step0.root","recreate");
  TNtuple * gen  = new TNtuple("Gen","",particleVars.data()); 
  TNtuple * reco = new TNtuple("Reco","",trackVars.data());

  std::cout << "starting skim loop" << std::endl;
  //Actual skimming
  int processed = 0;
  //cone size for filling density map
  float dMapR = 0.2;
  int nEtaBin = 192;
  int nPhiBin = 251;
  //grid resolution is 0.025x0.02503 in eta x phi space
  TH2D * densityMap = new TH2D("densityMap","densityMap:eta:phi",nEtaBin,-2.4,2.4,nPhiBin,-TMath::Pi(),TMath::Pi());
  
  for(int i = 0; i<10;i++)//trkCh->GetEntries(); i++)
  {
    if(s.nPb==2)  centCh->GetEntry(i);
    else trkCh->GetEntry(i);
  
    if((s.nPb==2) && ((hiBin/2 < centPUMin) || (hiBin/2 >= centPUMax))) continue;
    else if((s.nPb==0) && ((nVtx < centPUMin) || (nVtx >= centPUMax))) continue;
    else if(TMath::Abs(vz)>s.vz_window) continue;
    else if(processed%nSkip !=0)
    {
      processed++;
      continue;
    }
    if(s.nPb==2) trkCh->GetEntry(i);
  
    int centPU;
    if(s.nPb==2)
    {
      weight = getWeight(s,pthat,vz,hiBin);    
      centPU = hiBin;
    }
    if(s.nPb==0) 
    {
      weight = getWeight(s,pthat,vz,nVtx);
      centPU = nVtx;
    }

    //Filling density histogram
    //The eta phi strips are used simply to speed up program by only looping over the relevant 2*dMapR x 2*dMapR 'square' of bins
    //Can speed up by reversing loops and narrowing the Eta range a bit
    for(int j = 0; j<nTrk; j++)
    {
     if(TMath::Abs(trkEta[j])>2.4) continue;
     //loop over strip in eta around the track
     for(int eta = densityMap->GetXaxis()->FindBin(trkEta[j]-dMapR); eta<=densityMap->GetXaxis()->FindBin(trkEta[j]+dMapR); eta++)
     {
       //loop over strip in phi (have to be careful about the -pi to pi wrap around...)
       //for case where we don't have to worry about wrap around
       if(TMath::Pi()-TMath::Abs(trkPhi[j])>dMapR)
       {
         for(int phi = densityMap->GetYaxis()->FindBin(trkPhi[j]-dMapR); phi<=densityMap->GetYaxis()->FindBin(trkPhi[j]+dMapR); phi++)
         {
           if(TMath::Power(trkEta[j]-densityMap->GetXaxis()->GetBinCenter(eta),2)+TMath::Power(TMath::ACos(TMath::Cos(trkPhi[j]-densityMap->GetYaxis()->GetBinCenter(phi))),2)<dMapR*dMapR ) densityMap->SetBinContent(eta,phi,densityMap->GetBinContent(eta,phi)+1); 
         }
       }
       else
       //for case with -pi and pi wrap around 
       {
         for(int phi = 1; phi<=densityMap->GetYaxis()->FindBin(trkPhi[j]+dMapR-(trkPhi[j]>0?2*TMath::Pi():0)); phi++) 
         { 
           if(TMath::Power(trkEta[j]-densityMap->GetXaxis()->GetBinCenter(eta),2)+TMath::Power(TMath::ACos(TMath::Cos(trkPhi[j]-densityMap->GetYaxis()->GetBinCenter(phi))),2)<dMapR*dMapR ) densityMap->SetBinContent(eta,phi,densityMap->GetBinContent(eta,phi)+1); 
         }
         for(int phi = densityMap->GetYaxis()->FindBin(trkPhi[j]-dMapR+(trkPhi[j]<0?2*TMath::Pi():0)); phi<=nPhiBin; phi++) 
         { 
           if(TMath::Power(trkEta[j]-densityMap->GetXaxis()->GetBinCenter(eta),2)+TMath::Power(TMath::ACos(TMath::Cos(trkPhi[j]-densityMap->GetYaxis()->GetBinCenter(phi))),2)<dMapR*dMapR ) densityMap->SetBinContent(eta,phi,densityMap->GetBinContent(eta,phi)+1); 
         }
       }
     }
    }//end density map fill
    /*TCanvas * c1 = new TCanvas("c1","c1",800,800);
    densityMap->Draw("colz");
    c1->SaveAs("Density.png");
    delete c1;*/
//need to divide by area to get density still
   
 for(int j = 0; j<nTrk; j++)
    {
      if(TMath::Abs(trkEta[j])>2.4) continue;
      if(highPurity[j]!=1) continue;
      //TODO: Calo matching here
      //other cut here as well maybe?
      //trkStauts cut here?
      if(trkPt[j]<ptMin || trkPt[j]>ptMax) continue;

      localTrackDensity = (float)densityMap->GetBinContent(densityMap->GetXaxis()->FindBin(trkEta[j]),densityMap->GetYaxis()->FindBin(trkPhi[j]))/getArea(trkEta[j],dMapR);
      float trkEntry[] = {trkPt[j],trkEta[j],trkPhi[j],localTrackDensity,weight};
      reco->Fill(trkEntry);
    }
    for(int j = 0; j<nParticle; j++)
    {
      if(TMath::Abs(genEta[j])>2.4) continue;
      if(genPt[j]<ptMin || genPt[j]>ptMax) continue;

      localTrackDensity = (float)densityMap->GetBinContent(densityMap->GetXaxis()->FindBin(genEta[j]),densityMap->GetYaxis()->FindBin(genPhi[j]))/getArea(genEta[j],dMapR);
      float genEntry[] = {genPt[j],genEta[j],genPhi[j],localTrackDensity,weight};
      gen->Fill(genEntry); 
    }
    densityMap->Reset();
    processed++;
  }
  delete densityMap;
  skimOut->Write();
  skimOut->Close();   
}
