#ifndef GETWEIGHTS
#define GETWEIGHTS

#include <iostream>
#include "Settings.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"

bool isWeightFileInitialized = false;
TH1D * PthatWeight;
TH1D * VertexWeight;
TH1D * CentPUWeight; 

void produceWeights(Settings s)
{
  TH1::SetDefaultSumw2();
  std::cout << "Starting weight calculation" << std::endl;
  std::cout << "Calculating data distributions..." << std::endl;
  
  float pthat;
  float zVtx;
  int nVtx, hiBin;
  int pcoll, noiseFilter;
  int pHBHENoiseFilterResultProducer , pPAprimaryVertexFilter , pBeamScrapingFilter;
  int pClusterCompaitiblityFilter, pprimaryVertexFilter, phfCoincFilter3;

  TH1D * dVz;
  TH1D * dCentPU;
  TCanvas * c1 = new TCanvas("c1","",800,600);

  TFile * f;
  TTree * evtSel;
  TTree * centTree;
  TTree * trk;
  if(s.doVtx || s.doCentPU)
  {
    f = TFile::Open(s.DataFile.c_str(),"read");

    //stuff needed for event selections
    evtSel = (TTree*) f->Get("skimanalysis/HltTree");
    //evtSel->SetBranchAddress("pHBHENoiseFilterResultProducer",&pHBHENoiseFilterResultProducer);
    if(s.nPb==0)
    {
      evtSel->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
      evtSel->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);  
    }
    else if(s.nPb==2)
    {
      evtSel->SetBranchAddress("pClusterCompaitiblityFilter",&pClusterCompaitiblityFilter);  
      evtSel->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);  
      evtSel->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);  
    }

    //stuff needed for pileup RW
    trk = (TTree*) f->Get(Form("%s/trackTree",s.trackTreeName.c_str()));
    trk->SetBranchAddress("nVtx",&nVtx);
    trk->SetBranchAddress("zVtx",&zVtx);
    evtSel->AddFriend(trk);
    
    //stuff needed for centrality RW and vtx
    if(s.doCentPU && s.nPb==2)  
    {
      centTree = (TTree*) f->Get("hiEvtAnalyzer/HiTree");
      centTree->SetBranchAddress("hiBin",&hiBin);
      evtSel->AddFriend(centTree);
    }

    //Calculating data distributions
    std::string eventSelectionString;
    if(s.nPb==0) eventSelectionString = "pPAprimaryVertexFilter && pBeamScrapingFilter";
    else if(s.nPb==2) eventSelectionString = "pClusterCompaitiblityFilter && pprimaryVertexFilter && phfCoincFilter3"; 
    if(s.doVtx)
    {
      dVz = new TH1D("dVz",";vz;Events",30,-15,15);
      evtSel->Draw("zVtx>>dVz",Form("%s && TMath::Abs(zVtx)<%d",eventSelectionString.c_str(),s.vz_window));
      dVz->SetDirectory(0);
      dVz->Scale(1.0/(double)dVz->GetEntries());
      dVz->Draw();
      c1->SaveAs("../../evalPlots/dataVz.png");
    }
    if(s.doCentPU && s.nPb==2)
    {
      dCentPU = new TH1D("dCentPU",";hiBin;Events",100,0,200);
      evtSel->Draw("hiBin>>dCentPU",Form("%s && TMath::Abs(zVtx)<%d",eventSelectionString.c_str(),s.vz_window));
      dCentPU->SetDirectory(0);
      dCentPU->Scale(1.0/(double)dCentPU->GetEntries());
      dCentPU->Draw();
      c1->SaveAs("../../evalPlots/dataCentPU.png");
    }
    else if(s.doCentPU && s.nPb==0)
    {
      dCentPU = new TH1D("dCentPU",";nVtx;Events",30,0,30);
      evtSel->Draw("nVtx>>dCentPU",Form("%s && TMath::Abs(zVtx)<%d",eventSelectionString.c_str(),s.vz_window));
      dCentPU->SetDirectory(0);
      dCentPU->Scale(1.0/(double)dCentPU->GetEntries());
      dCentPU->Draw();
      c1->SaveAs("../../evalPlots/dataCentPU.png");
    }
    f->Close();
  }

  std::cout << "Calculating MC distributions..." << std::endl;
  TChain * jet;
  TChain * trkCh;
  TChain * centCh;
  TChain * evtCh;
  TH1D * MCVz;
  TH1D * MCCentPU;
  TH1D * MCPthat;
    
  if(s.doPthat || s.doVtx || s.doCentPU)
  {
    if(s.doPthat) MCPthat = new TH1D("MCPthat",";pthat;Events",320,0,800);
    if(s.doVtx) MCVz = new TH1D("MCVz",";vz;Events",30,-15,15);
    if(s.doCentPU && s.nPb==2) MCCentPU = new TH1D("MCCentPU",";hiBin;Events",100,0,200);
    if(s.doCentPU && s.nPb==0) MCCentPU = new TH1D("MCCentPU",";nVtx;Events",30,0,30);

    //stuff for pthat
    jet = new TChain(Form("%sJetAnalyzer/t",s.jetDefinition.c_str()));
    for(int i = 0; i<s.nMC; i++)  jet->Add(s.MCFiles.at(i).c_str());  
    jet->SetBranchAddress("pthat", &pthat);
     
    //stuff for pcoll
    evtCh = new TChain("skimanalysis/HltTree");
    for(int i = 0; i<s.nMC; i++)  evtCh->Add(s.MCFiles.at(i).c_str());
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
    jet->AddFriend(evtCh);
 
    //stuff needed for pileup RW
    trkCh = new TChain(Form("%s/trackTree",s.trackTreeName.c_str()));
    for(int i = 0; i<s.nMC; i++)  trkCh->Add(s.MCFiles.at(i).c_str());  
    trkCh->SetBranchAddress("nVtx",&nVtx);
    trkCh->SetBranchAddress("zVtx",&zVtx);
    jet->AddFriend(trkCh);
    
    //stuff needed for centrality RW and vtx
    if(s.doCentPU && s.nPb==2)
    {
      centCh = new TChain("hiEvtAnalyzer/HiTree");
      for(int i = 0; i<s.nMC; i++)  centCh->Add(s.MCFiles.at(i).c_str());  
      centCh->SetBranchAddress("hiBin",&hiBin);
      jet->AddFriend(centCh);
    }
    
    std::string eventSelectionString_MC;
    //to be used with updated samples
    //if(s.nPb==0) eventSelectionString_MC = "pPAprimaryVertexFilter && pBeamScrapingFilter";
    //else if(s.nPb==2) eventSelectionString_MC = "pprimaryVertexFilter && phfCoincFilter3";       
    if(s.nPb==0) eventSelectionString_MC = "1";
    else if(s.nPb==2) eventSelectionString_MC = "1";       
    //Calculating MC distributions
    if(s.doPthat)
    {
      MCPthat->GetDirectory()->cd();
      jet->Draw("pthat>>MCPthat", Form("TMath::Abs(zVtx)<%d && %s",s.vz_window,eventSelectionString_MC.c_str()));
      MCPthat->SetDirectory(0);
      MCPthat->Draw();
      c1->SaveAs("../../evalPlots/MCPthat.png");
    }
    if(s.doVtx)
    {
      MCVz->GetDirectory()->cd();
      jet->Draw("zVtx>>MCVz",Form("TMath::Abs(zVtx)<%d && %s",s.vz_window,eventSelectionString_MC.c_str()));
      MCVz->SetDirectory(0);
      MCVz->Scale(1.0/(double)MCVz->GetEntries());
      MCVz->Draw();
      c1->SaveAs("../../evalPlots/MCVz.png");
    }
    if(s.doCentPU && s.nPb==2)
    {
      MCCentPU->GetDirectory()->cd();
      jet->Draw("hiBin>>MCCentPU",Form("TMath::Abs(zVtx)<%d && %s",s.vz_window,eventSelectionString_MC.c_str()));
      MCCentPU->SetDirectory(0);
      MCCentPU->Scale(1.0/(double)MCCentPU->GetEntries());
      MCCentPU->Draw();
      c1->SaveAs("../../evalPlots/MCCentPU.png");
    }
    else if(s.doCentPU && s.nPb==0)
    {
      MCCentPU->GetDirectory()->cd();
      jet->Draw("nVtx>>MCCentPU",Form("TMath::Abs(zVtx)<%d && %s",s.vz_window,eventSelectionString_MC.c_str()));
      MCCentPU->SetDirectory(0);
      MCCentPU->Scale(1.0/(double)MCCentPU->GetEntries());
      MCCentPU->Draw();
      c1->SaveAs("../../evalPlots/MCCentPU.png");
    } 
  }
  
  //calculating pthat distributions
  TFile * out = TFile::Open(Form("%s_Weights.root",s.jobName.c_str()),"recreate");
  if(s.doPthat)
  {
    std::vector<int> numberOfEvents;
    for(int i = 1; i<MCPthat->GetSize()-1; i++)
    {
      for(int j = 0; j<s.nMC; j++)
      {
        if(i==1) numberOfEvents.push_back(0);
        if(MCPthat->GetBinCenter(i)>s.pthatBins.at(j) && MCPthat->GetBinCenter(i)<s.pthatBins.at(j+1)) numberOfEvents.at(j) += MCPthat->GetBinContent(i);
      }  
    }


    //outputting weight file 
    TH1D * pthatWeight = (TH1D*) MCPthat->Clone("pthatWeight");
    for(int i = 1; i<pthatWeight->GetSize()-1; i++)
    {
      for(int j = 0; j<s.nMC; j++)
      {
        if(pthatWeight->GetBinCenter(i)>s.pthatBins.at(j) && pthatWeight->GetBinCenter(i)<s.pthatBins.at(j+1))
        {
          pthatWeight->SetBinContent(i,s.pthatCrossSection.at(j)/numberOfEvents.at(j));
        }
        pthatWeight->SetBinError(i,0);
      }
    }
    pthatWeight->Draw();
    c1->SaveAs("../../evalPlots/pthatWeights.png");
    pthatWeight->Write();

    //outputting new distribution
    for(int i = 1; i<MCPthat->GetSize()-1; i++)
    {
      for(int j = 0; j<s.nMC; j++)
      {
        if(MCPthat->GetBinCenter(i)>s.pthatBins.at(j) && MCPthat->GetBinCenter(i)<s.pthatBins.at(j+1))
        {
          MCPthat->SetBinContent(i,MCPthat->GetBinContent(i)/numberOfEvents.at(j)*s.pthatCrossSection.at(j));
          MCPthat->SetBinError(i,MCPthat->GetBinError(i)/numberOfEvents.at(j)*s.pthatCrossSection.at(j));
        }
      }
    }
    c1->SetLogy();
    MCPthat->Draw();
    MCPthat->GetXaxis()->SetRangeUser(30,800);
    MCPthat->Draw();
    c1->SaveAs("../../evalPlots/PthatDist.png");
    c1->SetLogy(0);
  }

  if(s.doVtx)
  {
    dVz->Divide(MCVz);
    dVz->Draw();
    c1->SaveAs("../../evalPlots/VzWeight.png");
    dVz->SetName("VertexWeight");
    dVz->Write();
  }
  if(s.doCentPU)
  {
    dCentPU->Divide(MCCentPU);
    dCentPU->Draw();
    c1->SaveAs("../../evalPlots/CentPUWeight.png");
    dCentPU->SetName("CentPUWeight");
    dCentPU->Write();
  }

  out->Close();

  delete c1;
}

double getWeight(Settings s, double pthat, double vz, double centPU)
{
  if(!isWeightFileInitialized)
  {
    TFile *  weightFile = TFile::Open(Form("%s_Weights.root",s.jobName.c_str()),"read");
    if(s.doPthat)
    {
      PthatWeight = (TH1D*)weightFile->Get("pthatWeight");
      PthatWeight->SetDirectory(0);
    }
    if(s.doVtx)
    {
      VertexWeight = (TH1D*)weightFile->Get("VertexWeight");
      VertexWeight->SetDirectory(0);
    }
    if(s.doCentPU)
    {
      CentPUWeight = (TH1D*)weightFile->Get("CentPUWeight");   
      CentPUWeight->SetDirectory(0);
    }
    weightFile->Close();
    isWeightFileInitialized=true;
  }
  
  double weight = 1;
  if(s.doPthat)  weight = weight*PthatWeight->GetBinContent(PthatWeight->FindBin(pthat));
  if(s.doVtx)    weight = weight*VertexWeight->GetBinContent(VertexWeight->FindBin(vz));
  if(s.doCentPU) weight = weight*CentPUWeight->GetBinContent(CentPUWeight->FindBin(centPU));
  return weight;
}

#endif
