//important things needed have comments with !!!


#include "getTrkCorr.h"  //include this!!!
#include "TFile.h"
#include "TTree.h"
#include <iostream>



int example()
{
 TrkCorr t;              //make a TrkCorr Object!!!
 
 float pt[100000];
 float eta[100000];
 float phi[100000];
 int nTrk = -1;
 //int hiBin = -1;

 //need to change this to a pp file...
 TFile * f = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root","read");
 TTree * tree = (TTree*)f->Get("ppTrack/trackTree");
 tree->SetBranchAddress("trkPt",&pt);
 tree->SetBranchAddress("trkEta",&eta);
 tree->SetBranchAddress("trkPhi",&phi);
 tree->SetBranchAddress("nTrk",&nTrk);
 //TTree * hi = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
 //hi->SetBranchAddress("hiBin",&hiBin);
 //tree->AddFriend(hi);

 for(int i = 0; i<5; i++)
 {
   tree->GetEntry(i);
   t.UpdateEventInfo(pt,eta,phi,nTrk);           //call this in the event loop, feed it an array of pt,eta,phi, and nTrk for all the tracks in the event!!!
   for(int j=0; j<nTrk; j++)
   {
     if(j>10) break;
     std::cout << i << " " << j << " " << t.getTrkCorr(pt[j],eta[j],phi[j]) <<"\n"<< std::endl;  //call getTrkCorr to get the correction.!!!
   }
 }

}
