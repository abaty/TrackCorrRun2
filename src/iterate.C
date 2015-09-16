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


void iterate(Settings s,int iter, int stepType)
{
  float pt, eta, phi, density, weight, centPU, rmin, maxJetPt,trkStatus; 

  std::cout << "Loading appropriate information for denominator (gen for efficiency, reco for fake)..." << std::endl;
  TFile * histFile;
  if(iter==0) histFile = TFile::Open(Form("corrHists_job%d.root",s.job),"recreate");
  else        histFile = TFile::Open(Form("corrHists_job%d.root",s.job),"update");
  TH1D * genHist, *recoHist;
  TH2D * genHist2, *recoHist2;

  //make needed gen skim if it hasn't been done yet
  if(iter<s.nStep)
  {
    std::cout << "Denominator info not found, calculating and saving it..." << std::endl;
    if(stepType == 0 || stepType==2 || stepType == 4 || stepType == 6) 
    {
      genHist = makeTH1(s,stepType,"gen");
      recoHist= makeTH1(s,stepType,"mreco");
    }
    if(stepType == 1)
    {
      genHist2 = makeTH2(s,stepType,"gen");
      recoHist2= makeTH2(s,stepType,"mreco");
    }
      
    TFile * skim;
    //if(s.reuseSkim) skim = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/tracking_Efficiencies/ntuples/trackSkim_job%d.root",s.job),"read");
    //else skim = TFile::Open(Form("trackSkim_job%d.root",s.job),"read");
    skim = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/trackSkim_job%d.root",s.job),"read");
    //for efficiency
    TNtuple * gen = (TNtuple*)  skim->Get("Gen");
    gen->SetBranchAddress("genPt",&pt);
    gen->SetBranchAddress("genEta",&eta); 
    gen->SetBranchAddress("genPhi",&phi);
    gen->SetBranchAddress("genDensity",&density);
    gen->SetBranchAddress("weight",&weight);
    gen->SetBranchAddress("centPU",&centPU);
    gen->SetBranchAddress("rmin",&rmin);
    gen->SetBranchAddress("jtpt",&maxJetPt); 
	   
    for(int i = 0; i<gen->GetEntries(); i++)
    {
      gen->GetEntry(i);
      if(stepType==0) genHist->Fill(pt,weight);
      if(stepType==1) genHist2->Fill(eta,phi,weight); 
      if(stepType==2) genHist->Fill(centPU,weight);
      if(stepType==4) genHist->Fill(eta,weight); 
      if(stepType==6) genHist->Fill(density,weight);
    }

    //for fake
    TNtuple * reco = (TNtuple*)  skim->Get("Reco"); 
    reco->SetBranchAddress("trkPt",&pt);
    reco->SetBranchAddress("trkEta",&eta);
    reco->SetBranchAddress("trkPhi",&phi);
    reco->SetBranchAddress("trkDensity",&density);
    reco->SetBranchAddress("weight",&weight);
    reco->SetBranchAddress("centPU",&centPU);
    reco->SetBranchAddress("rmin",&rmin);
    reco->SetBranchAddress("jtpt",&maxJetPt);
    reco->SetBranchAddress("trkStatus",&trkStatus);
    for(int i = 0; i<reco->GetEntries(); i++)
    {
      reco->GetEntry(i);
      if(trkStatus<0) continue;
      if(stepType==0) recoHist->Fill(pt,weight);
      if(stepType==1) recoHist2->Fill(eta,phi,weight); 
      if(stepType==2) recoHist->Fill(centPU,weight);
      if(stepType==4) recoHist->Fill(eta,weight); 
      if(stepType==6) recoHist->Fill(density,weight);
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
  std::cout << "Efficiency denominator histogram available now." << std::endl;
  if(stepType==0) recoHist = (TH1D*)histFile->Get("mreco_pt");
  if(stepType==1) recoHist2 = (TH2D*)histFile->Get("mreco_accept"); 
  if(stepType==2) recoHist = (TH1D*)histFile->Get("mreco_centPU");
  if(stepType==4) recoHist = (TH1D*)histFile->Get("mreco_eta"); 
  if(stepType==6) recoHist = (TH1D*)histFile->Get("mreco_density");
  std::cout << "Fake denominator histogram available now." << std::endl;
	   
  //************************************************************************************************************
  std::cout << "Calculating numerator for efficiency/fake calculation..." << std::endl;
  TH1D *mrecoHist, *divHist, *frecoHist, *fdivHist;
  TH2D *mrecoHist2,*divHist2, *frecoHist2, *fdivHist2;
  //getting old eff histograms to calculate the updated efficiency (the number 30 is arbitrary, increase if more are needed)
  TH1D * previousEff[30], *previousFake[30];
  TH2D * previousEff2[30], *previousFake2[30];
  for(int i=0; i<iter; i++)
  {
    int type = s.stepOrder.at(i%s.nStep); 
    if(type==0 || type==2 || type==4 || type == 6)
    {
      previousEff[i] = (TH1D*)histFile->Get(Form("eff_step%d",i)); 
      previousFake[i] = (TH1D*)histFile->Get(Form("fake_step%d",i));
    }
    if(type==1)
    {  
      previousEff2[i] = (TH2D*)histFile->Get(Form("eff_step%d",i)); 
      previousFake2[i] = (TH2D*)histFile->Get(Form("fake_step%d",i));
    }
  }

  //setting up stuff for reading out of skim
  if(stepType == 0 || stepType==2 || stepType==4 || stepType == 6)
  {
    mrecoHist = makeTH1(s,stepType,Form("mreco_eff_step%d",iter));
    frecoHist = makeTH1(s,stepType,Form("reco_fake_step%d",iter));
  }
  if(stepType == 1)
  {
    mrecoHist2 = makeTH2(s,stepType,Form("mreco_eff_step%d",iter));
    frecoHist2 = makeTH2(s,stepType,Form("reco_fake_step%d",iter));
  }

  TFile * skim;
  //if(s.reuseSkim) skim = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/tracking_Efficiencies/ntuples/trackSkim_job%d.root",s.job),"read");
  //else skim = TFile::Open(Form("trackSkim_job%d.root",s.job),"read");
  skim = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/trackSkim_job%d.root",s.job),"read");
  TNtuple * reco = (TNtuple*)  skim->Get("Reco"); 
  reco->SetBranchAddress("trkPt",&pt);
  reco->SetBranchAddress("trkEta",&eta);
  reco->SetBranchAddress("trkPhi",&phi);
  reco->SetBranchAddress("trkDensity",&density);
  reco->SetBranchAddress("weight",&weight);
  reco->SetBranchAddress("centPU",&centPU);
  reco->SetBranchAddress("rmin",&rmin);
  reco->SetBranchAddress("jtpt",&maxJetPt);
  reco->SetBranchAddress("trkStatus",&trkStatus);

  //reading out of skim 
  for(int i = 0; i<reco->GetEntries(); i++)
  {
    //applying efficiencies from all previous steps
    float previousEffCorr = 1;
    float previousFakeCorr = 1; 
    reco->GetEntry(i);
    
    //fake part
    if(trkStatus==-99) continue;
    if(iter!=0)
    {
      for(int n = 0; n<iter; n++)
      {
        int type = s.stepOrder.at(n%s.nStep);
        if(type==0) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(pt));
        if(type==1) previousFakeCorr *= previousFake2[n]->GetBinContent(previousFake2[n]->GetXaxis()->FindBin(eta),previousFake2[n]->GetYaxis()->FindBin(phi));
        if(type==2) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(centPU));
        if(type==4) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(eta));
        if(type==6) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(density));
      } 
    }
    if(previousFakeCorr==0) std::cout <<  "\nWarning!!! A correction is going to infinity.  This usually indicates an empty bin somewhere, try using a coarser binning or more events! \n" << std::endl;
    //filling histograms
    if(stepType==0) frecoHist->Fill(pt,weight/previousFakeCorr);
    if(stepType==1) frecoHist2->Fill(eta,phi,weight/previousFakeCorr); 
    if(stepType==2) frecoHist->Fill(centPU,weight/previousFakeCorr);
    if(stepType==4) frecoHist->Fill(eta,weight/previousFakeCorr); 
    if(stepType==6) frecoHist->Fill(density,weight/previousFakeCorr);

    //eff part
    if(trkStatus<0) continue;
    if(iter!=0)
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
    if(stepType==0) mrecoHist->Fill(pt,weight/previousEffCorr);
    if(stepType==1) mrecoHist2->Fill(eta,phi,weight/previousEffCorr); 
    if(stepType==2) mrecoHist->Fill(centPU,weight/previousEffCorr);
    if(stepType==4) mrecoHist->Fill(eta,weight/previousEffCorr); 
    if(stepType==6) mrecoHist->Fill(density,weight/previousEffCorr);
  }
  skim->Close();
 
  //saving reco and efficiencies 
  std::cout << "Calculating updated Efficiency/Fake Rate and saving histograms" << std::endl;
  histFile->cd();    
  if(stepType==0 || stepType == 2 || stepType==4 || stepType==6)
  {
    divHist = (TH1D*)mrecoHist->Clone(Form("eff_step%d",iter));
    divHist->Divide(genHist);
    mrecoHist->Write();
    divHist->Write();
    
    fdivHist = (TH1D*)frecoHist->Clone(Form("fake_step%d",iter));
    fdivHist->Divide(recoHist);
    frecoHist->Write();
    fdivHist->Write(); 
  }
  if(stepType==1)
  {
    divHist2 = (TH2D*)mrecoHist2->Clone(Form("eff_step%d",iter));
    divHist2->Divide(genHist2);
    mrecoHist2->Write();
    divHist2->Write();
    
    fdivHist2 = (TH2D*)frecoHist2->Clone(Form("fake_step%d",iter));
    fdivHist2->Divide(recoHist2);
    frecoHist2->Write();
    fdivHist2->Write();
  }
  
  //*********************************************************************************************
  //writing final correction tables (consolidating multiple steps of the same variable)
  //once again 5 is an arbitrary number, increase if needed...
  if(iter>=s.nStep*s.fullIterations && s.terminateStep==stepType)
  {
    std::cout << "Consolidating into final histograms by multiplying out efficiencies/fake rates per variable" << std::endl;
    TH1D * finalEff[5], *finalFake[5];
    TH2D * finalEff2[5], *finalFake2[5];
    for(int i=0; i<s.nStep; i++)
    {
      int type = s.stepOrder.at(i%s.nStep); 
      if(type == 0 || type==2 || type==4 || type == 6)
      {
        finalEff[i] = (TH1D*)previousEff[i]->Clone(Form("finalEff_step%d",i));
        finalFake[i] = (TH1D*)previousFake[i]->Clone(Form("finalFake_step%d",i));
      }
      if(type == 1)
      {
        finalEff2[i] = (TH2D*)previousEff2[i]->Clone(Form("finalEff_step%d",i));
        finalFake2[i] = (TH2D*)previousFake2[i]->Clone(Form("finalFake_step%d",i));
      }
    }
    for(int n = s.nStep; n<iter; n++)
    {
      int type2 = s.stepOrder.at(n%s.nStep);
      if(type2==0)
      {
        finalEff[n%s.nStep]->Multiply(previousEff[n]);
        finalFake[n%s.nStep]->Multiply(previousFake[n]); 
      }
      if(type2==1)
      {
        finalEff2[n%s.nStep]->Multiply(previousEff2[n]);
        finalFake2[n%s.nStep]->Multiply(previousFake2[n]); 
      }
      if(type2==2)
      {
        finalEff[n%s.nStep]->Multiply(previousEff[n]);
        finalFake[n%s.nStep]->Multiply(previousFake[n]); 
      }
      if(type2==4)
      {
        finalEff[n%s.nStep]->Multiply(previousEff[n]);
        finalFake[n%s.nStep]->Multiply(previousFake[n]); 
      }
      if(type2==6)
      {
        finalEff[n%s.nStep]->Multiply(previousEff[n]); 
        finalFake[n%s.nStep]->Multiply(previousFake[n]); 
      }
    }  
    if(stepType == 0 || stepType==2 || stepType==4 || stepType == 6)
    {
      finalEff[iter%s.nStep]->Multiply((TH1D*)histFile->Get(Form("eff_step%d",iter)));
      finalFake[iter%s.nStep]->Multiply((TH1D*)histFile->Get(Form("fake_step%d",iter)));
    }
    if(stepType == 1)
    {
      finalEff2[iter%s.nStep]->Multiply((TH2D*)histFile->Get(Form("eff_step%d",iter)));
      finalFake2[iter%s.nStep]->Multiply((TH2D*)histFile->Get(Form("fake_step%d",iter)));
    }
    for(int i=0; i<s.nStep; i++)
    {
      int type = s.stepOrder.at(i%s.nStep); 
      if(type == 0 || type==2 || type==4 || type == 6){ finalEff[i]->Write(); finalFake[i]->Write();}
      if(type == 1){ finalEff2[i]->Write(); finalFake2[i]->Write();}
    }
  }
  histFile->Close();    
 
  std::cout << "Finished with iteration \n" << std::endl;
  return;
}
