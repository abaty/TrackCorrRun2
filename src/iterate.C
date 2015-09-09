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
  //set log spacing 
  if(stepType ==0)
  {
    const int ptBins = s.ptBinFine+1;
    double ptAxis[ptBins];
    for(int x = 0; x<ptBins;x++) ptAxis[x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(s.ptBinFine))) + TMath::Log10(s.ptMin));
    hist = new TH1D(Form("%s_pt",titlePrefix),";p_{T};",ptBins,s.ptMin,s.ptMax);
  }

  if(stepType ==2) 
  {
    if(s.nPb==2)  hist = new TH1D(Form("%s_centPU",titlePrefix),";hiBin;",s.centPUBinFine,s.centPUMin,s.centPUMax); 
    if(s.nPb==0)  hist = new TH1D(Form("%s_centPU",titlePrefix),";nVtx;",s.centPUBinFine,s.centPUMin,s.centPUMax); 
  }
  
  if(stepType ==4) hist = new TH1D(Form("%s_eta",titlePrefix),";eta;",s.etaBinFine,-2.4,2.4);
 
  if(stepType ==6)
  {
    const int densityBins = 25;
    double densityAxis[densityBins+1]={0.0001,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,440,480,520,600,5000};
    hist = new TH1D(Form("%s_density",titlePrefix),";trkDensity;",densityBins,densityAxis);
  }
  return hist;
}

TH2D * makeTH2(Settings s, int stepType, const char * titlePrefix)
{
  TH2D * hist;
  if(stepType ==1)  hist = new TH2D(Form("%s_accept",titlePrefix),";#eta;#phi;",s.etaBinFine,-2.4,2.4,s.phiBinFine,-TMath::Pi(),TMath::Pi());
  return hist;
}


void iterate(Settings s,int iter, int stepType, const char * effOrFake)
{
  float pt, eta, phi, density, weight, highPurity,trkFake, centPU; 

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
      if(stepType == 0 || stepType==2 || stepType == 4 || stepType == 6) genHist = makeTH1(s,stepType,"gen");
      if(stepType == 1) genHist2 = makeTH2(s,stepType,"gen");
      
      TFile * skim = TFile::Open("trackSkim.root","read");
      TNtuple * gen = (TNtuple*)  skim->Get("Gen");
      gen->SetBranchAddress("genPt",&pt);
      gen->SetBranchAddress("genEta",&eta); 
      gen->SetBranchAddress("genPhi",&phi);
      gen->SetBranchAddress("genDensity",&density);
      gen->SetBranchAddress("weight",&weight);
      gen->SetBranchAddress("centPU",&centPU);
   
      for(int i = 0; i<gen->GetEntries(); i++)
      {
        gen->GetEntry(i);
        if(stepType==0) genHist->Fill(pt,weight);
        if(stepType==1) genHist2->Fill(eta,phi,weight); 
        if(stepType==2) genHist->Fill(centPU,weight);
        if(stepType==4) genHist->Fill(eta,weight); 
        if(stepType==6) genHist->Fill(density,weight);
      }
      skim->Close();
      histFile->Write(); 
    }
  
    //redundant for first step, but needed if the gen file was made and saved previously 
    if(stepType==0) genHist = (TH1D*)histFile->Get("gen_pt");
    if(stepType==1) genHist2 = (TH2D*)histFile->Get("gen_accept"); 
    if(stepType==2) genHist = (TH1D*)histFile->Get("gen_centPU");
    if(stepType==4) genHist = (TH1D*)histFile->Get("gen_eta"); 
    if(stepType==6) genHist = (TH1D*)histFile->Get("gen_density");
    std::cout << "Denominator histogram available now." << std::endl;
  }
   
  //********************************************************************************
  std::cout << "Calculating numerator for efficiency calculation..." << std::endl;
  //getting old eff histograms to calculate the updated efficiency (the number 30 is arbitrary, increase if more are needed)
  TH1D * previousEff[30];
  TH2D * previousEff2[30];
  for(int i=0; i<iter; i++)
  {
    int type = s.stepOrder.at(i%s.nStep); 
    if(type==0 || type==2 || type==4 || type == 6) previousEff[i] = (TH1D*)histFile->Get(Form("eff_step%d",i)); 
    if(type==1)              previousEff2[i] = (TH2D*)histFile->Get(Form("eff_step%d",i)); 
  }

  //setting up stuff for reading out of skim
  if(stepType == 0 || stepType==2 || stepType==4 || stepType == 6) recoHist = makeTH1(s,stepType,Form("reco_step%d",iter));
  if(stepType == 1) recoHist2 = makeTH2(s,stepType,Form("reco_step%d",iter));

  TFile * skim = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/trackSkim_job%d.root",s.job),"read");
  TNtuple * reco = (TNtuple*)  skim->Get("Reco"); 
  reco->SetBranchAddress("trkPt",&pt);
  reco->SetBranchAddress("trkEta",&eta);
  reco->SetBranchAddress("trkPhi",&phi);
  reco->SetBranchAddress("trkDensity",&density);
  reco->SetBranchAddress("weight",&weight);
  reco->SetBranchAddress("centPU",&centPU);

  //reading out of skim 
  for(int i = 0; i<reco->GetEntries(); i++)
  {
    //applying efficiencies from all previous steps
    float previousEffCorr = 1; 
    reco->GetEntry(i);
    if(strcmp(effOrFake,"Eff")==0 && iter!=0)
    {
      for(int n = 0; n<iter; n++)
      {
        int type = s.stepOrder.at(n%s.nStep);
        if(type==0) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(pt));
        if(type==1) previousEffCorr *= previousEff2[n]->GetBinContent(previousEff2[n]->GetXaxis()->FindBin(eta),previousEff2[n]->GetYaxis()->FindBin(phi));
        if(type==2) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(centPU));
        if(type==4) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(eta));
        if(type==6) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(density)); 
      } 
    }
    if(previousEffCorr==0) std::cout <<  "\nWarning!!! A correction is going to infinity.  This usually indicates an empty bin somewhere, try using a coarser binning or more events! \n" << std::endl;
   
    //filling histograms
    if(stepType==0) recoHist->Fill(pt,weight/previousEffCorr);
    if(stepType==1) recoHist2->Fill(eta,phi,weight/previousEffCorr); 
    if(stepType==2) recoHist->Fill(centPU,weight/previousEffCorr);
    if(stepType==4) recoHist->Fill(eta,weight/previousEffCorr); 
    if(stepType==6) recoHist->Fill(density,weight/previousEffCorr);
  }
  skim->Close();
 
  //saving reco and efficiencies 
  std::cout << "Calculating updated Efficiency and saving histograms" << std::endl;
  histFile->cd();    
  if(stepType==0 || stepType == 2 || stepType==4 || stepType==6)
  {
    divHist = (TH1D*)recoHist->Clone(Form("eff_step%d",iter));
    divHist->Divide(genHist);
    recoHist->Write();
    divHist->Write();
  }
  if(stepType==1)
  {
    divHist2 = (TH2D*)recoHist2->Clone(Form("eff_step%d",iter));
    divHist2->Divide(genHist2);
    recoHist2->Write();
    divHist2->Write();
  }
  histFile->Close();    
 
  std::cout << "Finished with iteration \n" << std::endl;
  return;
}
