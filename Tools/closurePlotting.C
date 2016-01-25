#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TAttAxis.h"
#include "TCanvas.h"

void closurePlotting()
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  
  const int nFiles = 3;

  TFile * inf[nFiles];
  TH1D *effNum[nFiles],*effDen[nFiles],*fakeNum[nFiles],*fakeDen[nFiles],*effCorrNum[nFiles],*fakeCorrNum[nFiles];
  
  for(int i = 0; i<nFiles; i++)
  {
    if(i==0) inf[i]  = TFile::Open("Jan18_PbPb_iterative/outputClosures_jan18.root");
    if(i==1) inf[i]  = TFile::Open("Jan18_PbPb_Iterative_CrosscheckDifferentOrder/outputClosures_jan18_differentOrder.root");
    if(i==2) inf[i]  = TFile::Open("Jan18_PbPb_Simple/outputClosures_jan18_simple.root");

    effNum[i] = (TH1D*)inf[i]->Get("effNoCorr_pt");
    effDen[i] = (TH1D*)inf[i]->Get("gen_pt");
    effCorrNum[i] = (TH1D*)inf[i]->Get("effCorr_pt");
    fakeNum[i] = (TH1D*)inf[i]->Get("fakeNoCorr_pt");
    fakeDen[i] = (TH1D*)inf[i]->Get("mreco_pt");
    fakeCorrNum[i] = (TH1D*)inf[i]->Get("fakeCorr_pt");

    effNum[i]->Divide(effDen[i]);
    effCorrNum[i]->Divide(effDen[i]);
    fakeNum[i]->Divide(fakeDen[i]);
    fakeCorrNum[i]->Divide(fakeDen[i]);
    effNum[i]->SetDirectory(0);
    effCorrNum[i]->SetDirectory(0);
    fakeNum[i]->SetDirectory(0);
    fakeCorrNum[i]->SetDirectory(0);
    inf[i]->Close();
  }
  for(int i = 0; i<nFiles; i++)
  {
    effNum[i]->SetMarkerColor(i+1);
    effCorrNum[i]->SetMarkerColor(i+1);
    fakeNum[i]->SetMarkerColor(i+1);
    fakeCorrNum[i]->SetMarkerColor(i+1);
    effNum[i]->SetLineColor(i+1);
    effCorrNum[i]->SetLineColor(i+1);
    fakeNum[i]->SetLineColor(i+1);
    fakeCorrNum[i]->SetLineColor(i+1);
    if(i==0){
      effNum[i]->SetMarkerStyle(21);
      effCorrNum[i]->SetMarkerStyle(21);
      fakeNum[i]->SetMarkerStyle(21);
      fakeCorrNum[i]->SetMarkerStyle(21);
      effNum[i]->SetLineWidth(2);
      effCorrNum[i]->SetLineWidth(2);
      fakeNum[i]->SetLineWidth(2);
      fakeCorrNum[i]->SetLineWidth(2);
    }else{
      effNum[i]->SetLineWidth(1);
      effCorrNum[i]->SetLineWidth(1);
      fakeNum[i]->SetLineWidth(1);
      fakeCorrNum[i]->SetLineWidth(1);
    }
  }

  TLegend * leg = new TLegend(0.2,0.65,0.6,0.9);
  leg->AddEntry(effNum[0],"Iterative Table","p");
  leg->AddEntry(effNum[1],"Iterative (Different Order)","p");
  leg->AddEntry(effNum[2],"Simple Table","p");


  c1->SetLogx();
  effNum[0]->Draw();
  effNum[0]->GetYaxis()->SetTitle("Efficiency");
  //for(int i = 1; i<nFiles; i++) effNum[i]->Draw("same");
  c1->SaveAs("../evalPlots/efficiency.png");
  c1->SaveAs("../evalPlots/efficiency.pdf");
  c1->Clear();
  effCorrNum[0]->Draw();
  effCorrNum[0]->GetYaxis()->SetTitle("Efficiency Closure");
  for(int i = 1; i<nFiles; i++) effCorrNum[i]->Draw("same");
  leg->Draw("same");
  c1->SaveAs("../evalPlots/efficiencyClosure.png");
  c1->SaveAs("../evalPlots/efficiencyClosure.pdf");
  c1->Clear();
  fakeNum[0]->Draw();
  fakeNum[0]->GetYaxis()->SetTitle("1/(1-fake)");
  fakeNum[0]->GetYaxis()->SetRangeUser(0.8,1.2);
  //for(int i = 1; i<nFiles; i++) fakeNum[i]->Draw("same");
  c1->SaveAs("../evalPlots/fake.png");
  c1->SaveAs("../evalPlots/fake.pdf");
  c1->Clear();
  fakeCorrNum[0]->GetYaxis()->SetRangeUser(0.8,1.2);
  fakeCorrNum[0]->GetYaxis()->SetTitle("fake Closure");
  fakeCorrNum[0]->Draw();
  leg->Draw("same");
  for(int i = 1; i<nFiles; i++) fakeCorrNum[i]->Draw("same");
  c1->SaveAs("../evalPlots/fakeClosure.png");
  c1->SaveAs("../evalPlots/fakeClosure.pdf");
}
