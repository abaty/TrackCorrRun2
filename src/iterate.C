#include "Settings.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include <iostream>

TH1D * makeTH1(Settings s, int stepType, const char * titlePrefix)
{
  TH1D * hist;
  if(stepType ==0)  hist = new TH1D(Form("%s_ptStep",titlePrefix),";p_{T};",s.ptBinFine,s.ptMin,s.ptMax);
  if(stepType ==6)  hist = new TH1D(Form("%s_densityStep",titlePrefix),";trkDensity;",30,0,600);
  return hist;
}

TH2D * makeTH2(Settings s, int stepType, const char * titlePrefix)
{
  TH2D * hist;
  if(stepType ==1)  hist = new TH2D(Form("%s_accStep",titlePrefix),";#eta;#phi;",s.etaBinFine,-2.4,2.4,s.phiBinFine,TMath::Pi(),-TMath::Pi());
  return hist;
}


void iterate(Settings s,int iter, int stepType, const char * effOrFake)
{

  std::cout << "Loading appropriate information for denominator (gen for efficiency, reco for fake)..." << std::endl;
  TFile * genFile = TFile::Open("iterHists.root","update");
  TH1D * genHist;
  TH2D * genHist2;

  if(strcmp(effOrFake,"Eff")==0)
  {
    //make needed gen skim if it hasn't been done yet
    if(iter<s.nStep)
    {
      std::cout << "Denominator info not found, calculating and saving it..." << std::endl;
      if(stepType == 0 || stepType == 6) genHist = makeTH1(s,stepType,"gen");
      if(stepType == 1) genHist2 = makeTH2(s,stepType,"gen");
      
      TFile * skim = TFile::Open("trackSkim.root","read");
      TNtuple * gen = (TNtuple*)  skim->Get("Gen");
      float pt, eta, phi, density, weight; 
      gen->SetBranchAddress("genPt",&pt);
      gen->SetBranchAddress("genEta",&eta); 
      gen->SetBranchAddress("genPhi",&phi);
      gen->SetBranchAddress("genDensity",&density);
      gen->SetBranchAddress("weight",&weight);
   
      for(int i = 0; i<gen->GetEntries(); i++)
      {
        gen->GetEntry(i);
        if(stepType==0) genHist->Fill(pt,weight);
        if(stepType==1) genHist2->Fill(eta,phi,weight);  
        if(stepType==6) genHist->Fill(density,weight);
      }
      skim->Close();
      genFile->Write(); 
    }
   
    if(stepType==0) genHist = (TH1D*)genFile->get("gen_ptStep");
    if(stepType==1) genHist2 = (TH2D*)genFile->get("gen_accStep");  
    if(stepType==6) genHist = (TH1D*)genFile->get("gen_densityStep");
  }
   

 
    
  return;
}
