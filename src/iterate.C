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
  if(stepType ==0)  hist = new TH1D(Form("%s_pt",titlePrefix),";p_{T};",s.ptBinFine,s.ptMin,s.ptMax);
  if(stepType ==6)  hist = new TH1D(Form("%s_density",titlePrefix),";trkDensity;",30,0,600);
  return hist;
}

TH2D * makeTH2(Settings s, int stepType, const char * titlePrefix)
{
  TH2D * hist;
  if(stepType ==1)  hist = new TH2D(Form("%s_accept",titlePrefix),";#eta;#phi;",s.etaBinFine,-2.4,2.4,s.phiBinFine,TMath::Pi(),-TMath::Pi());
  return hist;
}


void iterate(Settings s,int iter, int stepType, const char * effOrFake)
{
  float pt, eta, phi, density, weight, highPurity,trkFake; 

  std::cout << "Loading appropriate information for denominator (gen for efficiency, reco for fake)..." << std::endl;
  TFile * histFile;
  if(iter==0) histFile = TFile::Open("iterHists.root","recreate");
  else        histFile = TFile::Open("iterHists.root","update");
  TH1D * genHist, *recoHist, *divHist;
  TH2D * genHist2,*recoHist2,*divHist2;

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
      histFile->Write(); 
    }
  
    //redundant for first step, but needed if the gen file was made and saved previously 
    if(stepType==0) genHist = (TH1D*)histFile->Get("gen_pt");
    if(stepType==1) genHist2 = (TH2D*)histFile->Get("gen_accept");  
    if(stepType==6) genHist = (TH1D*)histFile->Get("gen_density");
    std::cout << "Denominator histogram available now." << std::endl;
  }
   
  //********************************************************************************
  std::cout << "Calculating numerator for efficiency calculation..." << std::endl;
  if(stepType == 0 || stepType == 6) recoHist = makeTH1(s,stepType,Form("reco_step%d",iter));
  if(stepType == 1) recoHist2 = makeTH2(s,stepType,Form("reco_step%d",iter));

  TFile * skim = TFile::Open("trackSkim.root","read");
  TNtuple * reco = (TNtuple*)  skim->Get("Reco"); 
  reco->SetBranchAddress("trkPt",&pt);
  reco->SetBranchAddress("trkEta",&eta);
  reco->SetBranchAddress("trkPhi",&phi);
  reco->SetBranchAddress("trkDensity",&density);
  reco->SetBranchAddress("weight",&weight);
  
  for(int i = 0; i<reco->GetEntries(); i++)
  {
    float previousEffCorr = 1; 
    reco->GetEntry(i);
    if(strcmp(effOrFake,"Eff")==0 && iter!=0)
    {
    //put function here to get all previous steps eff corrections 
    }
    if(stepType==0) recoHist->Fill(pt,weight*previousEffCorr);
    if(stepType==1) recoHist2->Fill(eta,phi,weight*previousEffCorr);  
    if(stepType==6) recoHist->Fill(density,weight*previousEffCorr);
  }
  skim->Close();
  histFile->cd();    
  if(stepType==0 || stepType==6) recoHist->Write();
  else                           recoHist2->Write();
  histFile->Close(); 
    
  return;
}
