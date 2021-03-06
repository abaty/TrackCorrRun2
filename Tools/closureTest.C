#include "../src/TrkSettings.h"
#include "../src/getWeights.C"
#include "../src/getTrkCorr.h"
#include "chi2Reweighting/Chi2Corrector_PbPb.C"
#include "TMath.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TAxis.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include <cstring>
#include <vector>
#include <iostream>

/*step type key
 * 0 - pt
 * 1 - acceptance
 * 2 - hibin or pu
 * 3 - jt pt
 * 4 - eta
 * 5 - rmin
 * 7 - pt+eta
 * 8 - pt+cent*/

//settings for the histograms used
TH1D * makeTH1(TrkSettings s, int stepType, const char * titlePrefix)
{
  TH1D * hist;
  //set log spacing 
  if(stepType ==0)
  {
    const int ptBins = 46;
    double ptAxis[ptBins];
    ptAxis[0]=0.5;ptAxis[1]=0.6;ptAxis[2]=0.7;ptAxis[3]=0.8;ptAxis[4]=0.9;
    for(int i = 5; i<10; i++){
      s.ptMin = s.ptBinCoarse.at(i);
      s.ptMax = s.ptBinCoarse.at(i+1);
      if(s.ptMin<2) for(int x = 0; x<6;x++) ptAxis[5+x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(6))) + TMath::Log10(s.ptMin));
      else if(s.ptMin<3) for(int x = 0; x<4;x++) ptAxis[11+x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(4))) + TMath::Log10(s.ptMin));
      else for(int x = 0; x<10;x++) ptAxis[15+x+10*(i-7)] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(10))) + TMath::Log10(s.ptMin));
    }
    ptAxis[45]=s.ptMax;
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
  return hist;
}

TH2D * makeTH2(TrkSettings s, int stepType, const char * titlePrefix)
{
  TH2D * hist;
  //if(s.ptMin>=10){s.etaBinFine = s.etaBinFine/2; s.phiBinFine = s.phiBinFine/2;}
  if(stepType ==1)  hist = new TH2D(Form("%s_accept",titlePrefix),";#eta;#phi;",s.etaBinFine,-2.4,2.4,s.phiBinFine,-TMath::Pi(),TMath::Pi());
  if(stepType ==7 || stepType==8)
  {
    const int ptBins = 46;
    double ptAxis[ptBins];
    ptAxis[0]=0.5;ptAxis[1]=0.6;ptAxis[2]=0.7;ptAxis[3]=0.8;ptAxis[4]=0.9;
    for(int i = 5; i<10; i++){
      s.ptMin = s.ptBinCoarse.at(i);
      s.ptMax = s.ptBinCoarse.at(i+1);
      if(s.ptMin<2) for(int x = 0; x<6;x++) ptAxis[5+x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(6))) + TMath::Log10(s.ptMin));
      else if(s.ptMin<3) for(int x = 0; x<4;x++) ptAxis[11+x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(4))) + TMath::Log10(s.ptMin));
      else for(int x = 0; x<10;x++) ptAxis[15+x+10*(i-7)] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(10))) + TMath::Log10(s.ptMin));
    }
    ptAxis[45]=s.ptMax;
    if(stepType==8){
      const int centralityBins = 6;
      double centralityEdges[centralityBins+1] = {0,5,10,30,50,70,90};
      hist = new TH2D(Form("%s_centPt",titlePrefix),";Centrality;#pt;",centralityBins,centralityEdges,ptBins-1,ptAxis);
    }
    else hist = new TH2D(Form("%s_etaPt",titlePrefix),";#eta;#pt;",s.etaBinFine,-2.4,2.4,ptBins-1,ptAxis);
  }
  return hist;
}


void closureTest(const char * in, const char * out,TrkSettings s)
{
  Chi2Corrector_PbPb * chi2corr = new Chi2Corrector_PbPb();
  double trackEtaCut = 1;
  
  TrkCorr* trkCorr = new TrkCorr(Form("%s",in));
//Setup variables for skim
  TFile * inputFile;
  TTree * trkCh;
  TTree * centCh;
  TTree * evtCh;
  TTree * jet;

  //track
  int nTrk;
  int nEv;
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
  float trkChi2[60000];
  unsigned char trkNHit[60000];
  unsigned char trkNlayer[60000];
  unsigned char trkAlgo[60000];
  unsigned char trkOriginalAlgo[60000];
  unsigned char trkNdof[60000];
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
  int   mtrkNHit[60000];
  int   mtrkAlgo[60000];
  int   mtrkOriginalAlgo[60000];
  int   mtrkNlayer[60000];
  float mtrkChi2[60000];
  int   mtrkNdof[60000];
  float pNRec[75000];
  
  //event parameters
  int hiBin;
  float vz;
  float pthat;
  int nref;
  int ngen;
  float jtpt[100];
  float jtphi[100];
  float jteta[100];
  float rawpt[100];
  float chargedSum[100];
  float genpt[100];
  float weight = 1;
  int pcoll;
  //Histograms to hold stuff...
  TFile * outF = TFile::Open(Form("%s%s",in,out),"recreate");
  TH1D *genPre[20], *mrecoPre[20];
  TH2D *genPre2[20], *mrecoPre2[20];
  TH1D * EffNoCorr[10], *FakeNoCorr[10];
  TH2D * EffNoCorr2[10], *FakeNoCorr2[10];
  TH1D * EffCorr[10], *FakeCorr[10];
  TH2D * EffCorr2[10], *FakeCorr2[10];
  TH1D * FinalCorr[10];
  TH2D * FinalCorr2[10];
  TH1D * Eff[10], *Fake[10], *EffClosure[10], *FakeClosure[10], *Closure[10];
  TH2D * Eff2[10], *Fake2[10], *EffClosure2[10], *FakeClosure2[10], *Closure2[10];
  
  for(int i = 0; i<9; i++)
  {
    if(i==6) continue;
    if(i != 1 && i!=7 && i!=8)
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
    const int ptBins_corr = 46;
    double ptAxis_corr[ptBins_corr];
    ptAxis_corr[0]=0.5;ptAxis_corr[1]=0.6;ptAxis_corr[2]=0.7;ptAxis_corr[3]=0.8;ptAxis_corr[4]=0.9;
    for(int i = 5; i<10; i++){
      s.ptMin = s.ptBinCoarse.at(i);
      s.ptMax = s.ptBinCoarse.at(i+1);
      if(s.ptMin<2) for(int x = 0; x<6;x++) ptAxis_corr[5+x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(6))) + TMath::Log10(s.ptMin));
      else if(s.ptMin<3) for(int x = 0; x<4;x++) ptAxis_corr[11+x] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(4))) + TMath::Log10(s.ptMin));
      else for(int x = 0; x<10;x++) ptAxis_corr[15+x+10*(i-7)] = TMath::Power(10,(x*(TMath::Log10(s.ptMax)-TMath::Log10(s.ptMin))/((float)(10))) + TMath::Log10(s.ptMin));
    }
    ptAxis_corr[45]=s.ptMax;
    TH2D * appliedCorrection = new TH2D("appliedCorrection",";p_{T};Correction",ptBins_corr-1,ptAxis_corr,100,0,10);
    TH2D * appliedEffCorrection = new TH2D("appliedEffCorrection",";p_{T};Eff Correction",ptBins_corr-1,ptAxis_corr,100,0,10);
    TH2D * appliedFakeCorrection = new TH2D("appliedFakeCorrection",";p_{T};Fake Correction",ptBins_corr-1,ptAxis_corr,100,0,10);

  //Setup input trees  
  //track tree
  for(int nFiles = 0; nFiles<s.nMC; nFiles++){
  inputFile = TFile::Open(s.MCFiles.at(nFiles).c_str(),"read");
  trkCh = (TTree*) inputFile->Get(Form("%s/trackTree",s.trackTreeName.c_str()));
  trkCh->SetBranchAddress("nTrk",&nTrk); 
  trkCh->SetBranchAddress("nEv",&nEv); 
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
  trkCh->SetBranchAddress("trkNHit",&trkNHit);
  trkCh->SetBranchAddress("trkChi2",&trkChi2); 
  trkCh->SetBranchAddress("trkNlayer",&trkNlayer); 
  trkCh->SetBranchAddress("trkAlgo",&trkAlgo); 
  trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo); 
  trkCh->SetBranchAddress("trkNdof",&trkNdof); 
  trkCh->SetBranchAddress("mtrkNHit",&mtrkNHit); 
  trkCh->SetBranchAddress("mtrkChi2",&mtrkChi2); 
  trkCh->SetBranchAddress("mtrkNlayer",&mtrkNlayer); 
  trkCh->SetBranchAddress("mtrkAlgo",&mtrkAlgo); 
  trkCh->SetBranchAddress("mtrkOriginalAlgo",&mtrkOriginalAlgo); 
  trkCh->SetBranchAddress("mtrkNdof",&mtrkNdof);
  
  //centrality and vz
  centCh = (TTree*) inputFile->Get("hiEvtAnalyzer/HiTree");
  if(s.doCentPU && s.nPb==2) centCh->SetBranchAddress("hiBin",&hiBin);
  trkCh->AddFriend(centCh);  
  
  //pthat and jets
  jet = (TTree*) inputFile->Get(Form("%sJetAnalyzer/t",s.jetDefinition.c_str()));
  jet->SetBranchAddress("pthat", &pthat);
  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("ngen",&ngen);
  jet->SetBranchAddress("jtpt",&jtpt);
  jet->SetBranchAddress("jteta",&jteta);
  jet->SetBranchAddress("jtphi",&jtphi);
  jet->SetBranchAddress("rawpt",&rawpt);
  jet->SetBranchAddress("genpt",&genpt);
  jet->SetBranchAddress("chargedSum",&chargedSum);  
  trkCh->AddFriend(jet);
  
  //evtCh = new TChain("skimanalysis/HltTree");
  //for(int i = 0; i<s.nMC; i++)  evtCh->Add(s.MCFiles.at(i).c_str());
  //evtCh->SetBranchAddress("pcollisionEventSelection",&pcoll);
  //trkCh->AddFriend(evtCh);


  //**************************************************************************************************************************************************************
  //event loop
  std::cout << "starting event loop" << std::endl;
  int numberOfEntries = trkCh->GetEntries();
  numberOfEntries = trkCh->GetEntries();

  for(int i = 0; i<numberOfEntries; i++)
  { 
    if(i%50000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    if(s.nPb==2)  centCh->GetEntry(i);
    trkCh->GetEntry(i);
    if(s.doSplit && nEv%2==0) continue; 
    if(TMath::Abs(zVtx[0])>s.vz_window) continue;
    
    if(pthat>800) continue;
  
    //getting weight parameters
    float centPU;
    if(s.nPb==2)
    {
      weight = getWeight(s,pthat,zVtx[0],hiBin,in);    
      centPU = hiBin/2.0;
    }
    if(s.nPb==0) 
    {
      weight = getWeight(s,pthat,zVtx[0],nVtx,in);
      centPU = nVtx;
    }

	//max jet pt
    float maxJetPt = -999;
    for(int k = 0; k<nref; k++)
    {
      if(TMath::Abs(jteta[k])>2) continue;
      if(jtpt[k]>maxJetPt) maxJetPt=jtpt[k];
    }
    float maxGenJetPt = -999;
    for(int k = 0; k<ngen; k++)
    {
      if(genpt[k]>maxGenJetPt) maxGenJetPt = genpt[k];
    }
 
    //track loop 
    for(int j = 0; j<nTrk; j++)
    {
      if(trkPt[j]>pthat/1.5) continue;
      if(TMath::Abs(trkEta[j])>=trackEtaCut) continue;
      if(highPurity[j]!=1) continue;
      if( trkPtError[j]/trkPt[j]>0.3 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 ||TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue; 
      if(s.nPb==2 || s.doTrackTriggerCuts){
        if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>((s.nPb==2)?(1.0/chi2corr->getChi2Scale(hiBin,trkPt[j])):1)*0.15) continue; 
        if(trkNHit[j]<11 && trkPt[j]>0.7) continue;
        if( trkPtError[j]/trkPt[j]>0.1) continue; 
      }
      //if((maxJetPt>50 && trkPt[j]>maxJetPt) || (maxJetPt<50 && trkPt[j]>50)) continue;
      if(s.doTrackTriggerCuts && (trkNHit[j]<11 || (int)trkOriginalAlgo[j]<4 || (int)trkOriginalAlgo[j]>7)) continue; //track trigger cuts
        
      float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
      if(s.doCaloMatch && !(trkPt[j]<20 || (Et>0.5*trkPt[j]))) continue; //Calo Matching       
      if(trkPt[j]<0.5 || trkPt[j]>=400) continue;

      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        if(jtpt[k]<50) break;
        if(TMath::Abs(jteta[k])>2 || chargedSum[k]/rawpt[k]<0.01) continue;
        float R = TMath::Power(jteta[k]-trkEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[k]-trkPhi[j])),2);
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
      FakeNoCorr2[7]->Fill(trkEta[j],trkPt[j],weight);   
      FakeNoCorr2[8]->Fill(centPU,trkPt[j],weight);   

      float fake = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin,0,2);
      appliedFakeCorrection->Fill(trkPt[j],fake);
      FakeCorr[0]->Fill(trkPt[j],weight*fake);
      FakeCorr2[1]->Fill(trkEta[j],trkPhi[j],weight*fake);
      FakeCorr[2]->Fill(centPU,weight*fake);
      FakeCorr[3]->Fill(maxJetPt,weight*fake);
      FakeCorr[4]->Fill(trkEta[j],weight*fake);
      FakeCorr[5]->Fill(rmin,weight*fake);
      FakeCorr2[7]->Fill(trkEta[j],trkPt[j],weight*fake);   
      FakeCorr2[8]->Fill(centPU,trkPt[j],weight*fake);   

      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin,0);
      appliedCorrection->Fill(trkPt[j],correction);
      FinalCorr[0]->Fill(trkPt[j],weight*correction);
      FinalCorr2[1]->Fill(trkEta[j],trkPhi[j],weight*correction);
      FinalCorr[2]->Fill(centPU,weight*correction);
      FinalCorr[3]->Fill(maxJetPt,weight*correction);
      FinalCorr[4]->Fill(trkEta[j],weight*correction);
      FinalCorr[5]->Fill(rmin,weight*correction);
      FinalCorr2[7]->Fill(trkEta[j],trkPt[j],weight*correction);      
      FinalCorr2[8]->Fill(centPU,trkPt[j],weight*correction);      
 
      if(trkStatus[j]<-100) continue;
      mrecoPre[0]->Fill(trkPt[j],weight);
      mrecoPre2[1]->Fill(trkEta[j],trkPhi[j],weight); 
      mrecoPre[2]->Fill(centPU,weight);
      mrecoPre[3]->Fill(maxJetPt,weight);
      mrecoPre[4]->Fill(trkEta[j],weight); 
      mrecoPre[5]->Fill(rmin,weight);
      mrecoPre2[7]->Fill(trkEta[j],trkPt[j],weight);
      mrecoPre2[8]->Fill(centPU,trkPt[j],weight);
    }
 
    //gen 
    for(int j = 0; j<nParticle; j++)
    {
      if(genPt[j]>pthat/1.5) continue;
      if(TMath::Abs(genEta[j])>=trackEtaCut) continue;
      if(genPt[j]<0.5 || genPt[j]>400) continue;
    
      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        if(jtpt[k]<50) break;
        if(TMath::Abs(jteta[k])>2 || chargedSum[k]/rawpt[k]<0.01) continue;
        float R = TMath::Power(jteta[k]-genEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[k]-genPhi[j])),2);
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
      genPre2[7]->Fill(genEta[j],genPt[j],weight);
      genPre2[8]->Fill(centPU,genPt[j],weight);
	  
      //numerator for efficiency (number of gen tracks matched to highPurity track)
      if(mtrkPtError[j]/mtrkPt[j]>0.3 || TMath::Abs(mtrkDz1[j]/mtrkDzError1[j])>3 || TMath::Abs(mtrkDxy1[j]/mtrkDxyError1[j])>3) mtrkQual[j]=0;  
      
      if(s.nPb==2 || s.doTrackTriggerCuts){
        if(mtrkPtError[j]/mtrkPt[j]>0.1) mtrkQual[j]=0;  
        if(mtrkChi2[j]/(float)mtrkNdof[j]/(float)mtrkNlayer[j]>((s.nPb==2)?(1.0/chi2corr->getChi2Scale(hiBin,mtrkPt[j])):1)*0.15) mtrkQual[j]=0; 
        if(mtrkNHit[j]<11 && mtrkPt[j]>0.7) mtrkQual[j]=0;
      }
      //if((maxJetPt>50 && mtrkPt[j]>maxJetPt) || (maxJetPt<=50 && mtrkPt[j]>50)) mtrkQual[j]=0;
      if(s.doTrackTriggerCuts && (mtrkNHit[j]<11 || (int)mtrkOriginalAlgo[j]<4 || (int)mtrkOriginalAlgo[j]>7)) mtrkQual[j]=0;   //track trigger cuts
      float Et = (mtrkPfHcal[j]+mtrkPfEcal[j])/TMath::CosH(genEta[j]);
      if(s.doCaloMatch && !(mtrkPt[j]<20 || (Et>0.5*mtrkPt[j]))) mtrkQual[j]=0; //Calo Matching 
      
      if(mtrkQual[j]<1 || mtrkPt[j]<=0) continue;
      EffNoCorr[0]->Fill(genPt[j],weight);
      EffNoCorr2[1]->Fill(genEta[j],genPhi[j],weight);
      EffNoCorr[2]->Fill(centPU,weight);
      EffNoCorr[3]->Fill(maxJetPt,weight);
      EffNoCorr[4]->Fill(genEta[j],weight);
      EffNoCorr[5]->Fill(rmin,weight);
      EffNoCorr2[7]->Fill(genEta[j],genPt[j],weight);
      EffNoCorr2[8]->Fill(centPU,genPt[j],weight);
	  
      float eff = trkCorr->getTrkCorr(genPt[j],genEta[j],genPhi[j],hiBin,rmin,0,1);
      appliedEffCorrection->Fill(genPt[j],eff);
      EffCorr[0]->Fill(genPt[j],weight*eff);
      EffCorr2[1]->Fill(genEta[j],genPhi[j],weight*eff);
      EffCorr[2]->Fill(centPU,weight*eff);
      EffCorr[3]->Fill(maxJetPt,weight*eff);
      EffCorr[4]->Fill(genEta[j],weight*eff);
      EffCorr[5]->Fill(rmin,weight*eff);
      EffCorr2[7]->Fill(genEta[j],genPt[j],weight*eff);
      EffCorr2[8]->Fill(centPU,genPt[j],weight*eff);
    }
  }
  inputFile->Close();
  }

  //rebin high pt eta regions to reflect table bins
  const int finalEtaBins = 9;
  int etaBinsToMerge[finalEtaBins] = {4,3,2,2,2,2,2,3,4};
  int Lbin = 1;
  for(int i = 0; i<finalEtaBins;i++){
    for(int j = 1; j<genPre2[7]->GetYaxis()->GetNbins()+1;j++){
      if(genPre2[7]->GetYaxis()->GetBinLowEdge(j)<10) continue;
      float sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += genPre2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) genPre2[7]->SetBinContent(Lbin+s,j,sum);
      sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += EffNoCorr2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) EffNoCorr2[7]->SetBinContent(Lbin+s,j,sum);
      sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += EffCorr2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) EffCorr2[7]->SetBinContent(Lbin+s,j,sum);
      sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += mrecoPre2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) mrecoPre2[7]->SetBinContent(Lbin+s,j,sum);
      sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += FakeCorr2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) FakeCorr2[7]->SetBinContent(Lbin+s,j,sum);
      sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += FakeNoCorr2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) FakeNoCorr2[7]->SetBinContent(Lbin+s,j,sum);
      sum = 0;
      for(int s = 0; s<etaBinsToMerge[i]; s++) sum += FinalCorr2[7]->GetBinContent(Lbin+s,j);
      for(int s = 0; s<etaBinsToMerge[i]; s++) FinalCorr2[7]->SetBinContent(Lbin+s,j,sum);
    }
    Lbin += etaBinsToMerge[i];
  } 

  //do divisions
  for(int i = 0; i<9; i++)
  {
    if(i==6) continue;
    if(i != 1 && i!=7 && i!=8)
    {
      Eff[i] = (TH1D*)EffNoCorr[i]->Clone(Form("Efficiency_%d",i));
      Eff[i]->Divide(genPre[i]);
      Eff[i]->SetDirectory(outF);
      Fake[i] = (TH1D*)FakeNoCorr[i]->Clone(Form("Fake_%d",i));
      Fake[i]->Divide(mrecoPre[i]);
      Fake[i]->SetDirectory(outF);
      EffClosure[i] = (TH1D*)EffCorr[i]->Clone(Form("EfficiencyClosure_%d",i));
      EffClosure[i]->Divide(genPre[i]);
      EffClosure[i]->SetDirectory(outF);
      FakeClosure[i] = (TH1D*)FakeCorr[i]->Clone(Form("FakeClosure_%d",i));
      FakeClosure[i]->Divide(mrecoPre[i]);
      FakeClosure[i]->SetDirectory(outF);
      Closure[i] =(TH1D*)FinalCorr[i]->Clone(Form("netClosure_%d",i));
      Closure[i]->Divide(genPre[i]);
      Closure[i]->SetDirectory(outF);
    }
    else
    {
      Eff2[i] = (TH2D*)EffNoCorr2[i]->Clone(Form("Efficiency_%d",i));
      Eff2[i]->Divide(genPre2[i]);
      Eff2[i]->SetDirectory(outF);
      Fake2[i] = (TH2D*)FakeNoCorr2[i]->Clone(Form("Fake_%d",i));
      Fake2[i]->Divide(mrecoPre2[i]);
      Fake2[i]->SetDirectory(outF);
      EffClosure2[i] = (TH2D*)EffCorr2[i]->Clone(Form("EfficiencyClosure_%d",i));
      EffClosure2[i]->Divide(genPre2[i]);
      EffClosure2[i]->SetDirectory(outF);
      FakeClosure2[i] = (TH2D*)FakeCorr2[i]->Clone(Form("FakeClosure_%d",i));
      FakeClosure2[i]->Divide(mrecoPre2[i]);
      FakeClosure2[i]->SetDirectory(outF);
      Closure2[i] =(TH2D*)FinalCorr2[i]->Clone(Form("netClosure_%d",i));
      Closure2[i]->Divide(genPre2[i]);
      Closure2[i]->SetDirectory(outF);
    }
  }


  std::cout << "About to write" << std::endl;
  outF->Write();
  std::cout << "written" << std::endl;
  outF->Close();
}

void getClosure(const char * inputDirectory, const char * outputFile = "Closure.root")
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TrkSettings s(Form("%sTrkCorrInputFile.txt",inputDirectory));
  closureTest(inputDirectory,outputFile,s);
  return;
}
