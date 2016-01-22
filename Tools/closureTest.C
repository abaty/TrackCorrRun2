#include "../src/TrkSettings.h"
#include "../src/getWeights.C"
#include "getTrkCorr.h"
#include "TMath.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TAxis.h"
#include "TFile.h"
#include "TChain.h"
#include <cstring>
#include <vector>
#include <iostream>

//settings for the histograms used
TH1D * makeTH1(TrkSettings s, int stepType, const char * titlePrefix)
{
  TH1D * hist;
  //set log spacing 
  if(stepType ==0)
  {
    const int ptBins = 31;
    double ptAxis[ptBins];
    for(int x = 0; x<ptBins;x++) ptAxis[x] = TMath::Power(10,(x*(TMath::Log10(300)-TMath::Log10(0.5))/((float)(ptBins-1))) + TMath::Log10(0.5));
    hist = new TH1D(Form("%s_pt",titlePrefix),";p_{T};",ptBins-1,ptAxis);
  }

  if(stepType ==2) 
  {
    if(s.nPb==2)  hist = new TH1D(Form("%s_centPU",titlePrefix),";hiBin;",20,0,100); 
    if(s.nPb==0)  hist = new TH1D(Form("%s_centPU",titlePrefix),";nVtx;",20,0,100); 
  }
  
  if(stepType ==3) hist = new TH1D(Form("%s_maxJetPt",titlePrefix),";jtpt;",30,0,300); 
  if(stepType ==4) hist = new TH1D(Form("%s_eta",titlePrefix),";eta;",s.etaBinFine,-2.4,2.4);

  const int rminBins = 16;
  double rminBinning[rminBins+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,2,3,10};
  if(stepType ==5) hist = new TH1D(Form("%s_rmin",titlePrefix),";rmin;",rminBins,rminBinning);
  if(stepType ==6)
  {
    const int densityBins = 10;
    double R = 0.1;
    double densityAxis[densityBins+1]={0};
    densityAxis[0]=0; densityAxis[1]=1.0/(R*R*TMath::Pi()*2)-0.0001; densityAxis[densityBins]=100000;
    for(int i=2;i<6;i++)  densityAxis[i]=(i-1)/(R*R*TMath::Pi())+1.0/(R*R*TMath::Pi()*2)-0.0001;
    for(int i=6;i<densityBins;i++)  densityAxis[i]=(2*i-6)/(R*R*TMath::Pi())+1.0/(R*R*TMath::Pi()*2)-0.0001;
    hist = new TH1D(Form("%s_density",titlePrefix),";trkDensity;",densityBins,densityAxis);
  }
  return hist;
}

TH2D * makeTH2(TrkSettings s, int stepType, const char * titlePrefix)
{
  TH2D * hist;
  //if(s.ptMin>=10){s.etaBinFine = s.etaBinFine/2; s.phiBinFine = s.phiBinFine/2;}
  if(stepType ==1)  hist = new TH2D(Form("%s_accept",titlePrefix),";#eta;#phi;",s.etaBinFine,-2.4,2.4,s.phiBinFine,-TMath::Pi(),TMath::Pi());
  if(stepType ==7)
  {
    const int ptBins = 31;
    double ptAxis[ptBins];
    for(int x = 0; x<ptBins;x++) ptAxis[x] = TMath::Power(10,(x*(TMath::Log10(300)-TMath::Log10(0.5))/((float)(ptBins-1))) + TMath::Log10(0.5));
    hist = new TH2D(Form("%s_etaPt",titlePrefix),";#eta;#pt;",s.etaBinFine,-2.4,2.4,ptBins-1,ptAxis);
  }
  return hist;
}


void closureTest(TrkSettings s)
{
  TrkCorr* trkCorr = new TrkCorr("Jan18_PbPb_Iterative_CrosscheckDifferentOrder/");
//Setup variables for skim
  TChain * trkCh;
  TChain * centCh;
  TChain * evtCh;
  TChain * jet;

  //track
  int nTrk;
  float trkPt[75000];
  float trkPtError[50000];
  float trkEta[75000];
  float trkPhi[75000];
  float trkStatus[75000]; //for trkStatus, -999 = fake, -99 = secondary, 1 & 2 are matched tracks
  bool highPurity[75000];
  float trkMVA[75000];
  float pfHcal[75000];
  float pfEcal[75000];
  float trkDxy1[75000];
  float trkDxyError1[75000];
  float trkDz1[75000];
  float trkDzError1[75000];
  int nVtx;
  float zVtx[100];

  //gen parameters
  int nParticle;
  float genPt[75000];
  float mtrkPt[75000];
  float genEta[75000];
  float genPhi[75000];
  //int   mtrkQual[75000];
  bool   mtrkQual[75000];//for 5.02 samples
  float mtrkPtError[50000];
  float mtrkMVA[75000];
  float mtrkDxy1[100000];
  float mtrkDxyError1[100000];
  float mtrkDz1[100000];
  float mtrkDzError1[100000];
  float mtrkPfHcal[100000];
  float mtrkPfEcal[100000];
  float pNRec[75000];
  
  //event parameters
  int hiBin;
  float vz;
  float pthat;
  int nref;
  float jtpt[100];
  float jtphi[100];
  float jteta[100];
  float rawpt[100];
  float chargedSum[100];
  float weight = 1;
  int pcoll;

  //Setup input trees  
  //track tree     
  trkCh = new TChain(Form("%s/trackTree",s.trackTreeName.c_str()));
  for(int i = 0; i<s.nMC; i++)  trkCh->Add(s.MCFiles.at(i).c_str()); 
  trkCh->SetBranchAddress("nTrk",&nTrk); 
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkPtError",&trkPtError);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkMVA",&trkMVA);
  trkCh->SetBranchAddress("trkStatus",&trkStatus);
  trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
  trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
  trkCh->SetBranchAddress("trkDz1",&trkDz1);
  trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
  trkCh->SetBranchAddress("pfHcal",&pfHcal); 
  trkCh->SetBranchAddress("pfEcal",&pfEcal); 
  
  trkCh->SetBranchAddress("nParticle",&nParticle);
  trkCh->SetBranchAddress("pPt",&genPt);
  trkCh->SetBranchAddress("pEta",&genEta);
  trkCh->SetBranchAddress("pPhi",&genPhi);
  trkCh->SetBranchAddress("pNRec",&pNRec);
  trkCh->SetBranchAddress("mtrkPt",&mtrkPt);
  trkCh->SetBranchAddress("mtrkPtError",&mtrkPtError);
  //trkCh->SetBranchAddress("mtrkQual",&mtrkQual); //for 2.76 samples
  trkCh->SetBranchAddress("mhighPurity",&mtrkQual);  //for 5.02 samples
  trkCh->SetBranchAddress("mtrkMVA",&mtrkMVA);  //for 5.02 samples
  trkCh->SetBranchAddress("mtrkDxy1",&mtrkDxy1);
  trkCh->SetBranchAddress("mtrkDxyError1",&mtrkDxyError1);
  trkCh->SetBranchAddress("mtrkDz1",&mtrkDz1);
  trkCh->SetBranchAddress("mtrkDzError1",&mtrkDzError1);
  trkCh->SetBranchAddress("mtrkPfHcal",&mtrkPfHcal); 
  trkCh->SetBranchAddress("mtrkPfEcal",&mtrkPfEcal);
  trkCh->SetBranchAddress("nVtx",&nVtx);
  trkCh->SetBranchAddress("zVtx",&zVtx); 
  
  //centrality and vz
  centCh = new TChain("hiEvtAnalyzer/HiTree");
  for(int i = 0; i<s.nMC; i++)  centCh->Add(s.MCFiles.at(i).c_str());  
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
  jet->SetBranchAddress("rawpt",&rawpt);
  jet->SetBranchAddress("chargedSum",&chargedSum);  
  trkCh->AddFriend(jet);
  
  //evtCh = new TChain("skimanalysis/HltTree");
  //for(int i = 0; i<s.nMC; i++)  evtCh->Add(s.MCFiles.at(i).c_str());
  //evtCh->SetBranchAddress("pcollisionEventSelection",&pcoll);
  //trkCh->AddFriend(evtCh);

  //Histograms to hold stuff...
  TFile * outF = TFile::Open("outputClosures_jan18_differentOrder.root","recreate");
  TH1D *genPre[20], *mrecoPre[20];
  TH2D *genPre2[20], *mrecoPre2[20];
  TH1D * EffNoCorr[10], *FakeNoCorr[10];
  TH2D * EffNoCorr2[10], *FakeNoCorr2[10];
  TH1D * EffCorr[10], *FakeCorr[10];
  TH2D * EffCorr2[10], *FakeCorr2[10];
  TH1D * FinalCorr[10];
  TH2D * FinalCorr2[10];
  
  for(int i = 0; i<8; i++)
  {
    if(i != 1 && i!=7)
    {
      genPre[i] = makeTH1(s,i,"gen");
      mrecoPre[i] = makeTH1(s,i,"mreco");
      EffNoCorr[i] = makeTH1(s,i,"effNoCorr");
      FakeNoCorr[i] = makeTH1(s,i,"fakeNoCorr");
      EffCorr[i] = makeTH1(s,i,"effCorr");
      FakeCorr[i] = makeTH1(s,i,"fakeCorr");
      FinalCorr[i] = makeTH1(s,i,"finalCorr");
    }
    else
    {
      genPre2[i] = makeTH2(s,i,"gen");
      mrecoPre2[i]= makeTH2(s,i,"mreco");
      EffNoCorr2[i] = makeTH2(s,i,"effNoCorr");
      FakeNoCorr2[i] = makeTH2(s,i,"fakeNoCorr");
      EffCorr2[i] = makeTH2(s,i,"effCorr");
      FakeCorr2[i] = makeTH2(s,i,"fakeCorr");
      FinalCorr2[i] = makeTH2(s,i,"finalCorr");
    }
  }
 
  //booking applied corrections
    const int ptBins_corr = 31;
    double ptAxis_corr[ptBins_corr];
    for(int x = 0; x<ptBins_corr;x++) ptAxis_corr[x] = TMath::Power(10,(x*(TMath::Log10(300)-TMath::Log10(0.5))/((float)(ptBins_corr-1))) + TMath::Log10(0.5));
    TH2D * appliedCorrection = new TH2D("appliedCorrection",";p_{T};Correction",ptBins_corr-1,ptAxis_corr,100,0,10);
    TH2D * appliedEffCorrection = new TH2D("appliedEffCorrection",";p_{T};Eff Correction",ptBins_corr-1,ptAxis_corr,100,0,10);
    TH2D * appliedFakeCorrection = new TH2D("appliedFakeCorrection",";p_{T};Fake Correction",ptBins_corr-1,ptAxis_corr,100,0,10);

  //event loop
  std::cout << "starting event loop" << std::endl;
  int numberOfEntries = 2000;
  numberOfEntries = trkCh->GetEntries();
 
  for(int i = 0; i<numberOfEntries; i++)
  { 
    if(i%50000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    if(s.nPb==2)  centCh->GetEntry(i);
    else trkCh->GetEntry(i);
   
    //if(TMath::Abs(vz)>s.vz_window) continue;
    if(s.nPb==2) trkCh->GetEntry(i);
    if(pthat>800) continue;
  
    //getting weight parameters
    int centPU;
    if(s.nPb==2)
    {
      weight = getWeight(s,pthat,zVtx[0],hiBin);    
      centPU = hiBin/2.0;
    }
    if(s.nPb==0) 
    {
      weight = getWeight(s,pthat,zVtx[0],nVtx);
      centPU = nVtx;
    }

	//max jet pt
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
      if( trkPtError[j]/trkPt[j]>0.3 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 ||TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;  
        
      float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
      if(!(trkPt[j]<20 || (Et>0.2*trkPt[j] && Et>trkPt[j]-80))) continue; //Calo Matching 
      
      if(trkPt[j]<0.5 || trkPt[j]>=300) continue;

      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        if(jtpt[k]<50) break;
        if(TMath::Abs(jteta[k])>2 || chargedSum[k]/rawpt[k]<0.01) continue;
        float R = TMath::Power(jteta[k]-trkEta[j],2) + TMath::Power(jtphi[k]-trkPhi[j],2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }


      //fill histograms here
      //for fake
      FakeNoCorr[0]->Fill(trkPt[j],weight);
      FakeNoCorr2[1]->Fill(trkEta[j],trkPhi[j],weight);
      FakeNoCorr[2]->Fill(centPU,weight);
      FakeNoCorr[3]->Fill(maxJetPt,weight);
      FakeNoCorr[4]->Fill(trkEta[j],weight);
      FakeNoCorr[5]->Fill(rmin,weight);
      //FakeNoCorr[6]->Fill(density,weight);
      FakeNoCorr2[7]->Fill(trkEta[j],trkPt[j],weight);   

      float fake = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin,0,2);
      appliedFakeCorrection->Fill(trkPt[j],fake);
      FakeCorr[0]->Fill(trkPt[j],weight*fake);
      FakeCorr2[1]->Fill(trkEta[j],trkPhi[j],weight*fake);
      FakeCorr[2]->Fill(centPU,weight*fake);
      FakeCorr[3]->Fill(maxJetPt,weight*fake);
      FakeCorr[4]->Fill(trkEta[j],weight*fake);
      FakeCorr[5]->Fill(rmin,weight*fake);
      //FakeCorr[6]->Fill(density,weight*fake);
      FakeCorr2[7]->Fill(trkEta[j],trkPt[j],weight*fake);   

      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin,0);
      appliedCorrection->Fill(trkPt[j],correction);
      FinalCorr[0]->Fill(trkPt[j],weight*correction);
      FinalCorr2[1]->Fill(trkEta[j],trkPhi[j],weight*correction);
      FinalCorr[2]->Fill(centPU,weight*correction);
      FinalCorr[3]->Fill(maxJetPt,weight*correction);
      FinalCorr[4]->Fill(trkEta[j],weight*correction);
      FinalCorr[5]->Fill(rmin,weight*correction);
      //FinalCorr[6]->Fill(density,weight*correction);
      FinalCorr2[7]->Fill(trkEta[j],trkPt[j],weight*correction);      
 
      if(trkStatus[j]<-100) continue;
      mrecoPre[0]->Fill(trkPt[j],weight);
      mrecoPre2[1]->Fill(trkEta[j],trkPhi[j],weight); 
      mrecoPre[2]->Fill(centPU,weight);
      mrecoPre[3]->Fill(maxJetPt,weight);
      mrecoPre[4]->Fill(trkEta[j],weight); 
      mrecoPre[5]->Fill(rmin,weight);
      //mrecoPre[6]->Fill(density,weight);
      mrecoPre2[7]->Fill(trkEta[j],trkPt[j],weight);
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
        if(TMath::Abs(jteta[k])>2 || chargedSum[k]/rawpt[k]<0.01) continue;
        float R = TMath::Power(jteta[k]-genEta[j],2) + TMath::Power(jtphi[k]-genPhi[j],2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }
       
      //fill histograms
      //for eff/total
      genPre[0]->Fill(genPt[j],weight);
      genPre2[1]->Fill(genEta[j],genPhi[j],weight);
      genPre[2]->Fill(centPU,weight);
      genPre[3]->Fill(maxJetPt,weight);
      genPre[4]->Fill(genEta[j],weight);
      genPre[5]->Fill(rmin,weight);
      //genPre[6]->Fill(density,weight);
      genPre2[7]->Fill(genEta[j],genPt[j],weight);
	  
      //numerator for efficiency (number of gen tracks matched to highPurity track)
      if(mtrkPtError[j]/mtrkPt[j]>0.3 || TMath::Abs(mtrkDz1[j]/mtrkDzError1[j])>3 || TMath::Abs(mtrkDxy1[j]/mtrkDxyError1[j])>3) mtrkQual[j]=0;  
      float Et = (mtrkPfHcal[j]+mtrkPfEcal[j])/TMath::CosH(genEta[j]);
      if(!(mtrkPt[j]<20 || (Et>0.2*mtrkPt[j] && Et>mtrkPt[j]-80))) mtrkQual[j]=0; //Calo Matching 
      
      if(mtrkQual[j]<1 || mtrkPt[j]<=0) continue;
      EffNoCorr[0]->Fill(genPt[j],weight);
      EffNoCorr2[1]->Fill(genEta[j],genPhi[j],weight);
      EffNoCorr[2]->Fill(centPU,weight);
      EffNoCorr[3]->Fill(maxJetPt,weight);
      EffNoCorr[4]->Fill(genEta[j],weight);
      EffNoCorr[5]->Fill(rmin,weight);
      //EffNoCorr[6]->Fill(density,weight);
      EffNoCorr2[7]->Fill(genEta[j],genPt[j],weight);
	  
      float eff = trkCorr->getTrkCorr(genPt[j],genEta[j],genPhi[j],hiBin,rmin,0,1);
      appliedEffCorrection->Fill(genPt[j],eff);
      EffCorr[0]->Fill(genPt[j],weight*eff);
      EffCorr2[1]->Fill(genEta[j],genPhi[j],weight*eff);
      EffCorr[2]->Fill(centPU,weight*eff);
      EffCorr[3]->Fill(maxJetPt,weight*eff);
      EffCorr[4]->Fill(genEta[j],weight*eff);
      EffCorr[5]->Fill(rmin,weight*eff);
      //EffCorr[6]->Fill(density,weight*eff);
      EffCorr2[7]->Fill(genEta[j],genPt[j],weight*eff);
    }
  }
  std::cout << "About to write" << std::endl;
  outF->Write();
  std::cout << "written" << std::endl;
  outF->Close();
}

void getClosure()
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TrkSettings s("Jan18_PbPb_iterative/TrkCorrInputFile.txt");
  closureTest(s);
  return;
}
