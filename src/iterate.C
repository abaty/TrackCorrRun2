#include "TrkSettings.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include <iostream>

//settings for the histograms used
TH1D * makeTH1(TrkSettings s, int stepType, const char * titlePrefix)
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
  
  if(stepType ==3) hist = new TH1D(Form("%s_maxJetPt",titlePrefix),";jtpt;",30,0,300); 
  if(stepType ==4) hist = new TH1D(Form("%s_eta",titlePrefix),";eta;",s.etaBinFine,-2.4,2.4);

  const int rminBins = 16;
  double rminBinning[rminBins+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,2,3,10};
  if(stepType ==5) hist = new TH1D(Form("%s_rmin",titlePrefix),";rmin;",rminBins,rminBinning);
  /*if(stepType ==6)
  {
    const int densityBins = 10;
    double R = 0.1;
    double densityAxis[densityBins+1]={0};
    densityAxis[0]=0; densityAxis[1]=1.0/(R*R*TMath::Pi()*2)-0.0001; densityAxis[densityBins]=100000;
    for(int i=2;i<6;i++)  densityAxis[i]=(i-1)/(R*R*TMath::Pi())+1.0/(R*R*TMath::Pi()*2)-0.0001;
    for(int i=6;i<densityBins;i++)  densityAxis[i]=(2*i-6)/(R*R*TMath::Pi())+1.0/(R*R*TMath::Pi()*2)-0.0001;
    hist = new TH1D(Form("%s_density",titlePrefix),";trkDensity;",densityBins,densityAxis);
  }*/
  return hist;
}

TH2D * makeTH2(TrkSettings s, int stepType, const char * titlePrefix)
{
  TH2D * hist;
  if(s.ptMin>=10){s.etaBinFine = s.etaBinFine/2; s.phiBinFine = s.phiBinFine/2;}
  if(stepType ==1)  hist = new TH2D(Form("%s_accept",titlePrefix),";#eta;#phi;",s.etaBinFine,-2.4,2.4,s.phiBinFine,-TMath::Pi(),TMath::Pi());
  if(stepType ==7)
  {
    const int ptBins = s.ptBinFine+1;
    double ptAxis[ptBins];
    for(int x = 0; x<ptBins;x++) ptAxis[x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(s.ptBinFine))) + TMath::Log10(s.ptMin));
    hist = new TH2D(Form("%s_etaPt",titlePrefix),";#eta;#pt;",s.etaBinFine,-2.4,2.4,ptBins,s.ptMin,s.ptMax);
  }
  return hist;
}


//iteration code
void iterate(TrkSettings s,int iter, int stepType, bool doCondor)
{
  float pt, eta, phi, weight, centPU, rmin, maxJetPt,trkStatus,pNRec,mpt,mtrkQual; 
  
  TFile * histFile;
  if(iter==0) histFile = TFile::Open(Form("corrHists_job%d.root",s.job),"recreate");
  else        histFile = TFile::Open(Form("corrHists_job%d.root",s.job),"update");

  //make needed gen skim if it hasn't been done yet
  if(iter==0)
  {
    TH1D *genPre[20], *mrecoPre[20];
    TH2D *genPre2[20], *mrecoPre2[20];
    std::cout << "Denominator info not yet calculated; calculating and saving it..." << std::endl;
    for(int i = 0; i<8; i++)
    {
      if(i==6) continue;
      if(i != 1 && i!=7)
      {
        genPre[i] = makeTH1(s,i,"gen");
        mrecoPre[i] = makeTH1(s,i,"mreco");
      }
      else
      {
        genPre2[i] = makeTH2(s,i,"gen");
        mrecoPre2[i]= makeTH2(s,i,"mreco");
      }
    }
   
    TFile * skim;
    if(doCondor)
    {
      if(s.reuseSkim) skim = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/tracking_Efficiencies/ntuples/trackSkim_job%d.root",s.job),"read");
      else skim = TFile::Open(Form("trackSkim_job%d.root",s.job),"read");
    }
    else skim = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/trackSkim_job%d.root",s.job),"read");
    //for efficiency
    std::cout << "Doing Efficiency denominator" << std::endl;   
    TNtuple * gen = (TNtuple*)  skim->Get("Gen");
    gen->SetBranchAddress("genPt",&pt);
    gen->SetBranchAddress("genEta",&eta); 
    gen->SetBranchAddress("genPhi",&phi);
    //gen->SetBranchAddress("genDensity",&density);
    gen->SetBranchAddress("weight",&weight);
    gen->SetBranchAddress("centPU",&centPU);
    gen->SetBranchAddress("rmin",&rmin);
    gen->SetBranchAddress("jtpt",&maxJetPt);
    gen->SetBranchAddress("pNRec",&pNRec); 
	   
    for(int i = 0; i<gen->GetEntries(); i++)
    {
      gen->GetEntry(i);
      genPre[0]->Fill(pt,weight);
      genPre2[1]->Fill(eta,phi,weight); 
      genPre[2]->Fill(centPU,weight);
      genPre[3]->Fill(maxJetPt,weight);
      genPre[4]->Fill(eta,weight); 
      genPre[5]->Fill(rmin,weight);
      //genPre[6]->Fill(density,weight);
      genPre2[7]->Fill(eta,pt,weight);
    }

    //for fake
    std::cout << "Doing Fake denominator" << std::endl;   
    TNtuple * reco = (TNtuple*)  skim->Get("Reco"); 
    reco->SetBranchAddress("trkPt",&pt);
    reco->SetBranchAddress("trkEta",&eta);
    reco->SetBranchAddress("trkPhi",&phi);
    //reco->SetBranchAddress("trkDensity",&density);
    reco->SetBranchAddress("weight",&weight);
    reco->SetBranchAddress("centPU",&centPU);
    reco->SetBranchAddress("rmin",&rmin);
    reco->SetBranchAddress("jtpt",&maxJetPt);
    reco->SetBranchAddress("trkStatus",&trkStatus);
    for(int i = 0; i<reco->GetEntries(); i++)
    {
      reco->GetEntry(i);
      if(trkStatus<-100) continue;
      mrecoPre[0]->Fill(pt,weight);
      mrecoPre2[1]->Fill(eta,phi,weight); 
      mrecoPre[2]->Fill(centPU,weight);
      mrecoPre[3]->Fill(maxJetPt,weight);
      mrecoPre[4]->Fill(eta,weight); 
      mrecoPre[5]->Fill(rmin,weight);
      //mrecoPre[6]->Fill(density,weight);
      mrecoPre2[7]->Fill(eta,pt,weight);
    }
  
    //Secondary calculation (no iterations)
    std::cout << "Quickly calculating the Secondary Rate from the reco tree (No further iterations needed)" << std::endl;
    TH1D * Secondary_Matched = new TH1D("Secondary_Matched",";pt;",s.multiRecoBins.at(s.job/s.nPtBinCoarse),s.ptMin,s.ptMax); 
    TH1D * Secondary_Secondaries = new TH1D("Secondary_Secondaries",";pt;",s.multiRecoBins.at(s.job/s.nPtBinCoarse),s.ptMin,s.ptMax); 
    for(int i = 0; i<reco->GetEntries(); i++)
    {
      reco->GetEntry(i);
      if(trkStatus>-100)
      {
        Secondary_Matched->Fill(pt,weight);
        if(trkStatus==-99) Secondary_Secondaries->Fill(pt,weight);
      }
    }
    TH1D * Secondary = (TH1D*)Secondary_Secondaries->Clone("SecondaryRate");
    Secondary->Divide(Secondary_Matched);
    Secondary->SetDirectory(histFile);
    Secondary_Matched->SetDirectory(histFile);
    Secondary_Secondaries->SetDirectory(histFile);
    //end Secondary Reco calculation

    //multiReco calculation (no iterations)
    std::cout << "Quickly calculating the Multiple Reconstruction Rate from the gen tree (No further iterations needed)" << std::endl;
    TH1D * MultiGen = new TH1D("MultiGen",";pt;",s.multiRecoBins.at(s.job/s.nPtBinCoarse),s.ptMin,s.ptMax); 
    TH1D * MultiReco = new TH1D("MultiMatchedReco",";pt;",s.multiRecoBins.at(s.job/s.nPtBinCoarse),s.ptMin,s.ptMax); 
    for(int i = 0; i<gen->GetEntries(); i++)
    {
      gen->GetEntry(i);
      if(pNRec>-1)
      {
        MultiGen->Fill(pt,weight);
        if(pNRec>1) MultiReco->Fill(pt,(pNRec-1)*weight);
      }
    }
    TH1D * Multi = (TH1D*)MultiReco->Clone("MultipleRecoRate");
    Multi->Divide(MultiGen);
    Multi->SetDirectory(histFile);
    MultiReco->SetDirectory(histFile);
    MultiGen->SetDirectory(histFile);
    //end Multiple Reco calculation

    skim->Close();
    histFile->Write(); 
  }
	  
  //redundant for first step, but needed if the gen file was made and saved previously 
  std::cout << "Loading appropriate information for denominator (gen for efficiency, reco for fake)..." << std::endl;
  TH1D * genHist[20], *recoHist[20];
  TH2D * genHist2[20], *recoHist2[20];
  genHist[0] = (TH1D*)histFile->Get("gen_pt");
  genHist2[1] = (TH2D*)histFile->Get("gen_accept"); 
  genHist[2] = (TH1D*)histFile->Get("gen_centPU");
  genHist[3] = (TH1D*)histFile->Get("gen_maxJetPt");
  genHist[4] = (TH1D*)histFile->Get("gen_eta"); 
  genHist[5] = (TH1D*)histFile->Get("gen_rmin"); 
  //genHist[6] = (TH1D*)histFile->Get("gen_density");
  genHist2[7] = (TH2D*)histFile->Get("gen_etaPt");
  std::cout << "Efficiency denominator histogram available now." << std::endl;
  recoHist[0] = (TH1D*)histFile->Get("mreco_pt");
  recoHist2[1] = (TH2D*)histFile->Get("mreco_accept"); 
  recoHist[2] = (TH1D*)histFile->Get("mreco_centPU");
  recoHist[3] = (TH1D*)histFile->Get("mreco_maxJetPt");
  recoHist[4] = (TH1D*)histFile->Get("mreco_eta"); 
  recoHist[5] = (TH1D*)histFile->Get("mreco_rmin");
  //recoHist[6] = (TH1D*)histFile->Get("mreco_density");
  recoHist2[7] = (TH2D*)histFile->Get("mreco_etaPt");
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
    if(type==0 || type==2 || type==3 || type==4 || type==5)
    {
      previousEff[i] = (TH1D*)histFile->Get(Form("eff_step%d",i)); 
      previousFake[i] = (TH1D*)histFile->Get(Form("fake_step%d",i));
    }
    if(type==1 || type==7)
    {  
      previousEff2[i] = (TH2D*)histFile->Get(Form("eff_step%d",i)); 
      previousFake2[i] = (TH2D*)histFile->Get(Form("fake_step%d",i));
    }
  }

  //setting up stuff for reading out of skim
  if(stepType == 0 || stepType==2 || stepType==3 || stepType==4 || stepType==5)
  {
    mrecoHist = makeTH1(s,stepType,Form("mreco_eff_step%d",iter));
    frecoHist = makeTH1(s,stepType,Form("reco_fake_step%d",iter));
  }
  if(stepType == 1 || stepType == 7)
  {
    mrecoHist2 = makeTH2(s,stepType,Form("mreco_eff_step%d",iter));
    frecoHist2 = makeTH2(s,stepType,Form("reco_fake_step%d",iter));
  }

  TFile * skim;
  if(doCondor)
  {
    if(s.reuseSkim) skim = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/tracking_Efficiencies/ntuples/trackSkim_job%d.root",s.job),"read");
    else skim = TFile::Open(Form("trackSkim_job%d.root",s.job),"read");
  }
  else skim = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/trackSkim_job%d.root",s.job),"read");
  TNtuple * reco = (TNtuple*)  skim->Get("Reco"); 
  reco->SetBranchAddress("trkPt",&pt);
  reco->SetBranchAddress("trkEta",&eta);
  reco->SetBranchAddress("trkPhi",&phi);
//  reco->SetBranchAddress("trkDensity",&density);
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
    if(iter!=0)
    {
      for(int n = 0; n<iter; n++)
      {
        int type = s.stepOrder.at(n%s.nStep);
        if(type==0) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(pt));
        if(type==1) previousFakeCorr *= previousFake2[n]->GetBinContent(previousFake2[n]->GetXaxis()->FindBin(eta),previousFake2[n]->GetYaxis()->FindBin(phi));
        if(type==2) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(centPU));
        if(type==3) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(maxJetPt));
        if(type==4) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(eta));
        if(type==5) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(rmin));
//        if(type==6) previousFakeCorr *= previousFake[n]->GetBinContent(previousFake[n]->FindBin(density));
        if(type==7) previousFakeCorr *= previousFake2[n]->GetBinContent(previousFake2[n]->GetXaxis()->FindBin(eta),previousFake2[n]->GetYaxis()->FindBin(pt));
      } 
    }
    if(previousFakeCorr==0) std::cout <<  "\nWarning!!! A correction is going to infinity.  This usually indicates an empty bin somewhere, try using a coarser binning or more events! \n" << std::endl;
    //filling histograms
    if(stepType==0) frecoHist->Fill(pt,weight/previousFakeCorr);
    if(stepType==1) frecoHist2->Fill(eta,phi,weight/previousFakeCorr); 
    if(stepType==2) frecoHist->Fill(centPU,weight/previousFakeCorr);
    if(stepType==3) frecoHist->Fill(maxJetPt,weight/previousFakeCorr);
    if(stepType==4) frecoHist->Fill(eta,weight/previousFakeCorr); 
    if(stepType==5) frecoHist->Fill(rmin,weight/previousFakeCorr); 
//    if(stepType==6) frecoHist->Fill(density,weight/previousFakeCorr);
    if(stepType==7) frecoHist2->Fill(eta,pt,weight/previousFakeCorr);
  }
  
  //eff part
  TNtuple * gen = (TNtuple*)  skim->Get("Gen");
  gen->SetBranchAddress("genPt",&pt);
  gen->SetBranchAddress("genEta",&eta); 
  gen->SetBranchAddress("genPhi",&phi);
//  gen->SetBranchAddress("genDensity",&density);
  gen->SetBranchAddress("weight",&weight);
  gen->SetBranchAddress("centPU",&centPU);
  gen->SetBranchAddress("rmin",&rmin);
  gen->SetBranchAddress("jtpt",&maxJetPt);
  gen->SetBranchAddress("pNRec",&pNRec); 
  gen->SetBranchAddress("mtrkPt",&mpt);
  gen->SetBranchAddress("mtrkQual",&mtrkQual); 

  //reading out of skim 
  for(int i = 0; i<gen->GetEntries(); i++)
  {
    //applying efficiencies from all previous steps
    float previousEffCorr = 1;
    float previousFakeCorr = 1; 
    gen->GetEntry(i);
    if(mtrkQual<1||  mpt<=0) continue;
    if(iter!=0)
    {
      for(int n = 0; n<iter; n++)
      {
        int type = s.stepOrder.at(n%s.nStep);
        if(type==0) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(pt));
        if(type==1) previousEffCorr *= previousEff2[n]->GetBinContent(previousEff2[n]->GetXaxis()->FindBin(eta),previousEff2[n]->GetYaxis()->FindBin(phi));
        if(type==2) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(centPU));
        if(type==3) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(maxJetPt));
        if(type==4) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(eta));
        if(type==5) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(rmin));
//        if(type==6) previousEffCorr *= previousEff[n]->GetBinContent(previousEff[n]->FindBin(density));
        if(type==7) previousEffCorr *= previousEff2[n]->GetBinContent(previousEff2[n]->GetXaxis()->FindBin(eta),previousEff2[n]->GetYaxis()->FindBin(pt));
      } 
    }
    if(previousEffCorr==0) std::cout <<  "\nWarning!!! A correction is going to infinity.  This usually indicates an empty bin somewhere, try using a coarser binning or more events! \n" << std::endl; 
    //filling histograms
    if(stepType==0) mrecoHist->Fill(pt,weight/previousEffCorr);
    if(stepType==1) mrecoHist2->Fill(eta,phi,weight/previousEffCorr); 
    if(stepType==2) mrecoHist->Fill(centPU,weight/previousEffCorr);
    if(stepType==3) mrecoHist->Fill(maxJetPt,weight/previousEffCorr);
    if(stepType==4) mrecoHist->Fill(eta,weight/previousEffCorr); 
    if(stepType==5) mrecoHist->Fill(rmin,weight/previousEffCorr); 
//    if(stepType==6) mrecoHist->Fill(density,weight/previousEffCorr);
    if(stepType==7) mrecoHist2->Fill(eta,pt,weight/previousEffCorr);
  }
  skim->Close();
 
  //saving reco and efficiencies/fake rates 
  std::cout << "Calculating updated Efficiency/Fake Rate and saving histograms" << std::endl;
  histFile->cd();    
  if(stepType==0 || stepType == 2 || stepType == 3 || stepType==4 || stepType==5)
  {
    divHist = (TH1D*)mrecoHist->Clone(Form("eff_step%d",iter));
    divHist->Divide(genHist[stepType]);
    mrecoHist->Write();
    divHist->Write();
    
    fdivHist = (TH1D*)frecoHist->Clone(Form("fake_step%d",iter));
    fdivHist->Divide(recoHist[stepType]);
    frecoHist->Write();
    fdivHist->Write(); 
  }
  if(stepType==1 || stepType==7)
  {
    divHist2 = (TH2D*)mrecoHist2->Clone(Form("eff_step%d",iter));
    divHist2->Divide(genHist2[stepType]);
    mrecoHist2->Write();
    divHist2->Write();
    
    fdivHist2 = (TH2D*)frecoHist2->Clone(Form("fake_step%d",iter));
    fdivHist2->Divide(recoHist2[stepType]);
    frecoHist2->Write();
    fdivHist2->Write();
  }
  
  //*********************************************************************************************
  //*********************************************************************************************
  //Only executed on last step before exiting program
  if(iter>=s.nStep*s.fullIterations && s.terminateStep==stepType)
  {
    //writing final correction tables (consolidating multiple steps of the same variable)
    //once again 10 is an arbitrary number, increase if needed...
    std::cout << "Consolidating into final histograms by multiplying out efficiencies/fake rates per variable" << std::endl;
    TH1D * finalEff[10], *finalFake[10];
    TH2D * finalEff2[10], *finalFake2[10];
    for(int i=0; i<s.nStep; i++)
    {
      int type = s.stepOrder.at(i%s.nStep); 
      if(type == 0 || type==2 || type==3 || type==4 || type==5)
      {
        finalEff[i] = (TH1D*)previousEff[i]->Clone(Form("finalEff_type%d",i));
        finalFake[i] = (TH1D*)previousFake[i]->Clone(Form("finalFake_type%d",i));
      }
      if(type == 1 ||  type==7)
      {
        finalEff2[i] = (TH2D*)previousEff2[i]->Clone(Form("finalEff_type%d",i));
        finalFake2[i] = (TH2D*)previousFake2[i]->Clone(Form("finalFake_type%d",i));
      }
    }
    for(int n = s.nStep; n<iter; n++)
    {
      int type2 = s.stepOrder.at(n%s.nStep);
      if(type2==0 || type2==2 || type2==3 || type2==4 || type2==5)
      {
        finalEff[n%s.nStep]->Multiply(previousEff[n]);
        finalFake[n%s.nStep]->Multiply(previousFake[n]); 
      }
      if(type2==1 || type2==7)
      {
        finalEff2[n%s.nStep]->Multiply(previousEff2[n]);
        finalFake2[n%s.nStep]->Multiply(previousFake2[n]); 
      }
    }  
    if(stepType == 0 || stepType==2 || stepType==3 || stepType==4 || stepType==5)
    {
      finalEff[iter%s.nStep]->Multiply((TH1D*)histFile->Get(Form("eff_step%d",iter)));
      finalFake[iter%s.nStep]->Multiply((TH1D*)histFile->Get(Form("fake_step%d",iter)));
    }
    if(stepType == 1 || stepType == 7)
    {
      finalEff2[iter%s.nStep]->Multiply((TH2D*)histFile->Get(Form("eff_step%d",iter)));
      finalFake2[iter%s.nStep]->Multiply((TH2D*)histFile->Get(Form("fake_step%d",iter)));
    }
    for(int i=0; i<s.nStep; i++)
    {
      int type = s.stepOrder.at(i%s.nStep); 
      if(type == 0 || type==2 || type==3 || type==4 || type == 5){ finalEff[i]->Write(); finalFake[i]->Write();}
      if(type == 1 || type == 7){ finalEff2[i]->Write(); finalFake2[i]->Write();}
    }
   
    //******************************************************************************************************* 
    //writing calculating final closures in each variable checked after applying corrections
    std::cout << "Calculating Final Closures..." << std::endl;
    TH1D * finalEffClosure[10], *finalFakeClosure[10];
    TH2D * finalEffClosure2[10], *finalFakeClosure2[10];
    for(int i=0; i<8; i++)
    {
      int type = i;
      if(type==0 || type==2 || type==3 || type==4 || type==5)
      {
        finalEffClosure[i] = (TH1D*)genHist[i]->Clone(Form("finalEffClosure_var%d",i));
        finalFakeClosure[i] = (TH1D*)recoHist[i]->Clone(Form("finalFakeClosure_var%d",i));
        finalEffClosure[i]->Reset();
        finalFakeClosure[i]->Reset();
      }
      if(type==1 || type==7)
      {  
        finalEffClosure2[i] = (TH2D*)genHist2[i]->Clone(Form("finalEffClosure_var%d",i));
        finalFakeClosure2[i] = (TH2D*)recoHist2[i]->Clone(Form("finalFakeClosure_var%d",i));
        finalEffClosure2[i]->Reset();
        finalFakeClosure2[i]->Reset();
      }
    }
   
    if(doCondor)
    { 
      if(s.reuseSkim) skim = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/tracking_Efficiencies/ntuples/trackSkim_job%d.root",s.job),"read");
      else skim = TFile::Open(Form("trackSkim_job%d.root",s.job),"read");
    }
    else skim = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/trackSkim_job%d.root",s.job),"read");
    reco = (TNtuple*)  skim->Get("Reco"); 
    reco->SetBranchAddress("trkPt",&pt);
    reco->SetBranchAddress("trkEta",&eta);
    reco->SetBranchAddress("trkPhi",&phi);
    //reco->SetBranchAddress("trkDensity",&density);
    reco->SetBranchAddress("weight",&weight);
    reco->SetBranchAddress("centPU",&centPU);
    reco->SetBranchAddress("rmin",&rmin);
    reco->SetBranchAddress("jtpt",&maxJetPt);
    reco->SetBranchAddress("trkStatus",&trkStatus);

    //reading out of skim 
    for(int i = 0; i<reco->GetEntries(); i++)
    {
      //applying efficiencies from all previous steps
      float previousFakeCorr = 1; 
      reco->GetEntry(i);
      for(int n=0; n<s.nStep; n++)//getting correction
      {
        int type = s.stepOrder.at(n%s.nStep);
        if(type==0) previousFakeCorr *= finalFake[n]->GetBinContent(finalFake[n]->FindBin(pt));
        if(type==1) previousFakeCorr *= finalFake2[n]->GetBinContent(finalFake2[n]->GetXaxis()->FindBin(eta),finalFake2[n]->GetYaxis()->FindBin(phi));
        if(type==2) previousFakeCorr *= finalFake[n]->GetBinContent(finalFake[n]->FindBin(centPU));
        if(type==3) previousFakeCorr *= finalFake[n]->GetBinContent(finalFake[n]->FindBin(maxJetPt));
        if(type==4) previousFakeCorr *= finalFake[n]->GetBinContent(finalFake[n]->FindBin(eta));
        if(type==5) previousFakeCorr *= finalFake[n]->GetBinContent(finalFake[n]->FindBin(rmin));
      //  if(type==6) previousFakeCorr *= finalFake[n]->GetBinContent(finalFake[n]->FindBin(density));
        if(type==7) previousFakeCorr *= finalFake2[n]->GetBinContent(finalFake2[n]->GetXaxis()->FindBin(eta),finalFake2[n]->GetYaxis()->FindBin(pt));
      }
      if(previousFakeCorr<1) previousFakeCorr==1;
      finalFakeClosure[0]->Fill(pt,weight/previousFakeCorr);
      finalFakeClosure2[1]->Fill(eta,phi,weight/previousFakeCorr); 
      finalFakeClosure[2]->Fill(centPU,weight/previousFakeCorr);
      finalFakeClosure[3]->Fill(maxJetPt,weight/previousFakeCorr);
      finalFakeClosure[4]->Fill(eta,weight/previousFakeCorr); 
      finalFakeClosure[5]->Fill(rmin,weight/previousFakeCorr); 
     // finalFakeClosure[6]->Fill(density,weight/previousFakeCorr);  
      finalFakeClosure2[7]->Fill(eta,pt,weight/previousFakeCorr);  
    }
  
    gen = (TNtuple*)  skim->Get("Gen");
    gen->SetBranchAddress("genPt",&pt);
    gen->SetBranchAddress("genEta",&eta); 
    gen->SetBranchAddress("genPhi",&phi);
    //gen->SetBranchAddress("genDensity",&density);
    gen->SetBranchAddress("weight",&weight);
    gen->SetBranchAddress("centPU",&centPU);
    gen->SetBranchAddress("rmin",&rmin);
    gen->SetBranchAddress("jtpt",&maxJetPt);
    gen->SetBranchAddress("pNRec",&pNRec); 
    gen->SetBranchAddress("mtrkPt",&mpt);
    gen->SetBranchAddress("mtrkQual",&mtrkQual); 
    //reading out of skim 
    for(int i = 0; i<gen->GetEntries(); i++)
    {
      //applying efficiencies from all previous steps
      float previousEffCorr = 1;
      gen->GetEntry(i);
      if(mtrkQual<1 || mpt<=0) continue;
      for(int n=0; n<s.nStep; n++)//getting correction
      {
        int type = s.stepOrder.at(n%s.nStep);
        if(type==0) previousEffCorr *= finalEff[n]->GetBinContent(finalEff[n]->FindBin(pt));
        if(type==1) previousEffCorr *= finalEff2[n]->GetBinContent(finalEff2[n]->GetXaxis()->FindBin(eta),finalEff2[n]->GetYaxis()->FindBin(phi));
        if(type==2) previousEffCorr *= finalEff[n]->GetBinContent(finalEff[n]->FindBin(centPU));
        if(type==3) previousEffCorr *= finalEff[n]->GetBinContent(finalEff[n]->FindBin(maxJetPt));
        if(type==4) previousEffCorr *= finalEff[n]->GetBinContent(finalEff[n]->FindBin(eta));
        if(type==5) previousEffCorr *= finalEff[n]->GetBinContent(finalEff[n]->FindBin(rmin));
      //  if(type==6) previousEffCorr *= finalEff[n]->GetBinContent(finalEff[n]->FindBin(density));  
        if(type==7) previousEffCorr *= finalEff2[n]->GetBinContent(finalEff2[n]->GetXaxis()->FindBin(eta), finalEff2[n]->GetYaxis()->FindBin(pt));  
      }
      if(previousEffCorr>1) previousEffCorr==1;
      finalEffClosure[0]->Fill(pt,weight/previousEffCorr);
      finalEffClosure2[1]->Fill(eta,phi,weight/previousEffCorr); 
      finalEffClosure[2]->Fill(centPU,weight/previousEffCorr);
      finalEffClosure[3]->Fill(maxJetPt,weight/previousEffCorr);
      finalEffClosure[4]->Fill(eta,weight/previousEffCorr); 
      finalEffClosure[5]->Fill(rmin,weight/previousEffCorr); 
      //finalEffClosure[6]->Fill(density,weight/previousEffCorr);
      finalEffClosure2[7]->Fill(eta,pt,weight/previousEffCorr);
    }

    histFile->cd();
    for(int i=0; i<8; i++)
    {
      if(i==6) continue;
      if(i!=1 && i!=7)
      {
        finalFakeClosure[i]->Divide(recoHist[i]);
        finalEffClosure[i]->Divide(genHist[i]);
        finalFakeClosure[i]->Write();
        finalEffClosure[i]->Write();
      }
      else
      {
        finalFakeClosure2[i]->Divide(recoHist2[i]);
        finalEffClosure2[i]->Divide(genHist2[i]);
        finalFakeClosure2[i]->Write();
        finalEffClosure2[i]->Write();
      }
    }
    std::cout << "Calculating Final Closures..." << std::endl;
  }
  histFile->Close();    
 
  std::cout << "Finished with iteration \n" << std::endl;
  return;
}
