#include "TrkSettings.h"
#include "getWeights.C"
#include "TMath.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "chi2Reweighting/Chi2Corrector_PbPb.C"
#include <cstring>
#include <vector>
#include <iostream>

void makeSkim(TrkSettings s, bool doCondor)
{
  Chi2Corrector_PbPb * chi2corr = new Chi2Corrector_PbPb();
  std::cout << "\nJob number: " << s.job << "\nCorresponds to the following parameters\nnSkip: " << s.nSkip
  << "\nptMin: " << s.ptMin << "\nptMax: " << s.ptMax << "\ncentMin: " << s.centPUMin << "\ncentMax " << s.centPUMax << std::endl;

  std::cout << "\nCreating initial skim of important information" << std::endl;

//Setup variables for skim
  TTree * trkCh;
  TTree * centCh;
  TTree * evtCh;
  TTree * jet;

  //track
  int nTrk;
  float trkPt[60000];
  float trkEta[60000];
  float trkPhi[60000];
  float trkStatus[60000]; //for trkStatus, -999 = fake, -99 = secondary, 1 & 2 are matched tracks
  bool highPurity[60000];
  float trkMVA[60000];
  float trkPtError[60000];
  float trkDxy1[60000];
  float trkDxyError1[60000];
  float trkDz1[60000];
  float trkDzError1[60000];
  float pfHcal[60000];
  float pfEcal[60000];
  float trkChi2[60000];
  unsigned char trkNHit[60000];
  unsigned char trkNlayer[60000];
  unsigned char trkAlgo[60000];
  unsigned char trkOriginalAlgo[60000];
  unsigned char trkNdof[60000];
  int nEv;
  int nVtx;

  //gen parameters
  int nParticle;
  float genPt[60000];
  float mtrkPt[60000];
  float genEta[60000];
  float genPhi[60000];
  bool   mtrkQual[60000];//for 5.02 samples
  float   mtrkMVA[60000];
  float   mtrkPtError[60000];
  int   mtrkNHit[60000];
  float mtrkDxy1[60000];
  float mtrkDxyError1[60000];
  float mtrkDz1[60000];
  float mtrkDzError1[60000];
  float mtrkPfHcal[60000];
  float mtrkPfEcal[60000];
  int   mtrkAlgo[60000];
  int   mtrkOriginalAlgo[60000];
  int   mtrkNlayer[60000];
  float mtrkChi2[60000];
  int   mtrkNdof[60000];
  float pNRec[60000];

  
  //event parameters
  int hiBin;
  float zVtx[100];
  float pthat;
  int nref;
  int ngen;
  float jtpt[100];
  float jtphi[100];
  float jteta[100];
  float rawpt[100];
  float genpt[100];
  float chargedSum[100];
  float weight = 1;
  int pHBHENoiseFilterResultProducer , pPAprimaryVertexFilter , pBeamScrapingFilter;
  int pClusterCompaitiblityFilter, pprimaryVertexFilter, phfCoincFilter3;

  //Setup input trees  
  //track tree    
  TFile * inputFile;
  TFile * skimOut;
  float etaCut = 2.4;
  //Setup output Ntuples
  std::string trackVars;
  std::string particleVars;
  particleVars="genPt:genEta:genPhi:weight:centPU:rmin:jtpt:pNRec:mtrkPt:mtrkQual:nEv";
  trackVars=   "trkPt:trkEta:trkPhi:weight:centPU:rmin:jtpt:trkStatus:nEv";
  std::string ifPP = "";
  if(s.nPb==0) ifPP = "pp_";
  if(doCondor) skimOut = TFile::Open(Form("%strackSkim_job%d.root",ifPP.c_str(),s.job),"recreate");
  else         skimOut = TFile::Open(Form("/export/d00/scratch/abaty/trackingEff/ntuples/%strackSkim_job%d.root",ifPP.c_str(),s.job),"recreate");
  TNtuple * gen  = new TNtuple("Gen","",particleVars.data()); 
  TNtuple * reco = new TNtuple("Reco","",trackVars.data());

  std::cout << "starting skim loop" << std::endl;
  for(int nFile = 0; nFile<s.nMC; nFile++){ 
  inputFile = TFile::Open(s.MCFiles.at(nFile).c_str(),"read");
  trkCh = (TTree*)inputFile->Get(Form("%s/trackTree",s.trackTreeName.c_str()));
  trkCh->SetBranchAddress("nTrk",&nTrk); 
  trkCh->SetBranchAddress("nEv",&nEv); 
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkStatus",&trkStatus);
  //trkCh->SetBranchAddress("trkMVA",&trkMVA);
  trkCh->SetBranchAddress("trkPtError",&trkPtError);
  trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
  trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
  trkCh->SetBranchAddress("trkDz1",&trkDz1);
  trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
  trkCh->SetBranchAddress("nVtx",&nVtx);
  trkCh->SetBranchAddress("zVtx",&zVtx); 
  
  trkCh->SetBranchAddress("nParticle",&nParticle);
  trkCh->SetBranchAddress("pPt",&genPt);
  trkCh->SetBranchAddress("pEta",&genEta);
  trkCh->SetBranchAddress("pPhi",&genPhi);
  trkCh->SetBranchAddress("pNRec",&pNRec);
  trkCh->SetBranchAddress("mtrkPt",&mtrkPt);
  trkCh->SetBranchAddress("mhighPurity",&mtrkQual); 
  //trkCh->SetBranchAddress("mtrkMVA",&mtrkMVA); 
  trkCh->SetBranchAddress("mtrkPtError",&mtrkPtError); 
  trkCh->SetBranchAddress("mtrkDxy1",&mtrkDxy1);
  trkCh->SetBranchAddress("mtrkDxyError1",&mtrkDxyError1);
  trkCh->SetBranchAddress("mtrkDz1",&mtrkDz1);
  trkCh->SetBranchAddress("mtrkDzError1",&mtrkDzError1);
  trkCh->SetBranchAddress("trkChi2",&trkChi2); 
  trkCh->SetBranchAddress("trkNlayer",&trkNlayer); 
  trkCh->SetBranchAddress("trkNdof",&trkNdof); 
  trkCh->SetBranchAddress("mtrkChi2",&mtrkChi2); 
  trkCh->SetBranchAddress("mtrkNdof",&mtrkNdof);
  trkCh->SetBranchAddress("mtrkNlayer",&mtrkNlayer); 
  trkCh->SetBranchAddress("trkNHit",&trkNHit);
  trkCh->SetBranchAddress("mtrkNHit",&mtrkNHit); 
  if(s.doTrackTriggerCuts)
  {
    trkCh->SetBranchAddress("trkAlgo",&trkAlgo); 
    trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo); 
    trkCh->SetBranchAddress("mtrkAlgo",&mtrkAlgo); 
    trkCh->SetBranchAddress("mtrkOriginalAlgo",&mtrkOriginalAlgo); 
  } 
  if(s.doCaloMatch)
  {
    trkCh->SetBranchAddress("pfHcal",&pfHcal); 
    trkCh->SetBranchAddress("pfEcal",&pfEcal); 
    trkCh->SetBranchAddress("mtrkPfHcal",&mtrkPfHcal); 
    trkCh->SetBranchAddress("mtrkPfEcal",&mtrkPfEcal); 
  }

  //centrality
  if(s.doCentPU && s.nPb==2)
  {
    centCh = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    centCh->SetBranchAddress("hiBin",&hiBin);
    trkCh->AddFriend(centCh);  
  }
  
  //pthat and jets
  jet = (TTree*)inputFile->Get(Form("%sJetAnalyzer/t",s.jetDefinition.c_str()));
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
  
  evtCh = (TTree*)inputFile->Get("skimanalysis/HltTree");
  if(s.nPb==0)
  {
    evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
    evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);  
  }
  else if(s.nPb==2)
  {
    evtCh->SetBranchAddress("pClusterCompaitiblityFilter",&pClusterCompaitiblityFilter);  
    evtCh->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);  
    evtCh->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);  
  }
  trkCh->AddFriend(evtCh);


  //Actual skimming
  int processed = 0;
 
  int numberOfEntries = 0;
  if(doCondor) numberOfEntries = trkCh->GetEntries();
  else         numberOfEntries = trkCh->GetEntries()/8; 
  for(int i = 0; i<numberOfEntries; i++)
  {
  //  if(i%100!=0) continue;
    if(i%500==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    if(s.nPb==2)  centCh->GetEntry(i);
    trkCh->GetEntry(i);
   
    //some event selections on centrality, vz, or just throwing away some events because stats not needed 
    if((s.nPb==2) && ((hiBin/2.0 < s.centPUMin) || (hiBin/2.0 >= s.centPUMax))) continue;
    if((s.nPb==0) && ((nVtx < s.centPUMin) || (nVtx >= s.centPUMax))) continue;
    if(TMath::Abs(zVtx[0])>s.vz_window) continue;
    if(processed%(s.nSkip) !=0)
    {
      processed++;
      continue;
    }
    if(pthat>800) continue;
    //if(s.nPb==0 && (pPAprimaryVertexFilter==0 || pBeamScrapingFilter==0)) continue;
    //if(s.nPb==2 && (pClusterCompaitiblityFilter==0 || pprimaryVertexFilter==0 || phfCoincFilter3==0)) continue;
    if(s.nPb==2 && (pprimaryVertexFilter==0 || phfCoincFilter3==0)) continue;
  
    //getting weight parameters
    float centPU;
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
      if(TMath::Abs(trkEta[j])>etaCut) continue;
      if(trkPt[j]<s.ptMin || trkPt[j]>=s.ptMax) continue;
      if(highPurity[j]!=1) continue;
      if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
      if(trkPtError[j]/trkPt[j]>0.3) continue;
      if(s.nPb==2 || s.doTrackTriggerCuts){
        if(trkPtError[j]/trkPt[j]>0.1) continue;
        if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>((s.nPb==2)?(1.0/chi2corr->getChi2Scale(hiBin,trkPt[j])):1)*0.15) continue;
        if(trkNHit[j]<11 && trkPt[j]>0.7) continue;
      }
      //if((maxJetPt>50 && trkPt[j]>maxJetPt) || (maxJetPt<=50 && trkPt[j]>50)) continue;
      if(s.doCaloMatch)
      {
        float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
        if(!(trkPt[j]<20 || (Et>0.5*trkPt[j]))) continue; //Calo Matching       
      }
      
      if(s.doTrackTriggerCuts && ((int)trkOriginalAlgo[j]<4 || (int)trkOriginalAlgo[j]>7)) continue; //track trigger cuts


      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        //40 for fragmetnation functions
        //if(jtpt[k]<50) break;
        if(jtpt[k]<40) break;
        if(chargedSum[k]/rawpt[k]<0.01 || TMath::Abs(jteta[k])>2) continue;
        float R = TMath::Power(jteta[k]-trkEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[k]-trkPhi[j])),2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }
      //for fragmentation functions
      if(rmin>0.3) continue;

      
      float trkEntry[] = {trkPt[j],trkEta[j],trkPhi[j],weight,(float)centPU,rmin,maxJetPt,(float)trkStatus[j],(float)(nEv%2)};
      reco->Fill(trkEntry);
    }
 
    //gen 
    for(int j = 0; j<nParticle; j++)
    {
      if(genPt[j]>pthat/1.5) continue;
      if(TMath::Abs(genEta[j])>etaCut) continue;
      if(genPt[j]<s.ptMin || genPt[j]>=s.ptMax) continue;
      if(mtrkPtError[j]/mtrkPt[j]>0.3 || TMath::Abs(mtrkDz1[j]/mtrkDzError1[j])>3 || TMath::Abs(mtrkDxy1[j]/mtrkDxyError1[j])>3) mtrkQual[j]=0;  

      if(s.nPb==2 || s.doTrackTriggerCuts){
        if(mtrkPtError[j]/mtrkPt[j]>0.1) mtrkQual[j]=0;  
        if(mtrkChi2[j]/(float)mtrkNdof[j]/(float)mtrkNlayer[j]>((s.nPb==2)?(1.0/chi2corr->getChi2Scale(hiBin,mtrkPt[j])):1)*0.15) mtrkQual[j]=0;
        if(mtrkNHit[j]<11 && mtrkPt[j]>0.7) mtrkQual[j]=0;
      }
      //if((maxJetPt>50 && mtrkPt[j]>maxJetPt) || (maxJetPt<=50 && mtrkPt[j]>50)) mtrkQual[j]=0;
      if(s.doCaloMatch)
      {
        float Et = (mtrkPfHcal[j]+mtrkPfEcal[j])/TMath::CosH(genEta[j]);
        if(!(mtrkPt[j]<20 || (Et>0.5*mtrkPt[j]))) mtrkQual[j]=0; //Calo Matching 
      }
      if(s.doTrackTriggerCuts && ((int)mtrkOriginalAlgo[j]<4 || (int)mtrkOriginalAlgo[j]>7)) mtrkQual[j]=0;   //track trigger cuts

      //find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        //40 for fragmetnation functions
        //if(jtpt[k]<50) break;
        if(jtpt[k]<40) break;
        if(chargedSum[k]/rawpt[k]<0.01 || TMath::Abs(jteta[k])>2) continue;
        float R = TMath::Power(jteta[k]-genEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[k]-genPhi[j])),2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }
      //for fragmentation functions
      if(rmin>0.3) continue;

      float genEntry[] = {genPt[j],genEta[j],genPhi[j],weight,(float)centPU,rmin,maxJetPt,pNRec[j],mtrkPt[j],(float)mtrkQual[j],(float)(nEv%2)};
      gen->Fill(genEntry); 
    }
  processed++;
  }
  inputFile->Close(); 
  }//end file loop

  std::cout << "Writing skim..." << std::endl;
  skimOut->Write();
  skimOut->Close();   
  std::cout << "Done skimming." << std::endl;
}
