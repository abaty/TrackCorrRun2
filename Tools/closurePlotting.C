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
#include <cstring>

void closurePlotting()
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  float etaLimit = 1.0;
  float minPtLimit = 0.72;
  float maxPtLimit = 390;

  TCanvas * c1 = new TCanvas("c1","c1",800,600);
 
  /* 
  TLegend * leg = new TLegend(0.2,0.65,0.6,0.9);
  leg->AddEntry(effNum[0],"Iterative Table","p");
  leg->AddEntry(effNum[1],"Iterative (Different Order)","p");
  leg->AddEntry(effNum[2],"Simple Table","p");
  */
  const int nFiles = 2;
  std::string inDir[nFiles];
  inDir[0] = "TrkCorr_Mar15_Iterative_PbPb";
  inDir[1] = "TrkCorr_May12_Iterative_RpPbCuts";
  /*inDir[0] = "Jan_24_Iterative";
  inDir[1] = "Jan_24_Iterative_Alternate";
  inDir[2] = "Jan_24_Simple";
  inDir[3] = "Jan_24_pp_Iterative";
  inDir[4] = "Jan_24_pp_Iterative_Alternate";
  inDir[5] = "Jan_24_pp_Simple";
  inDir[6] = "Jan_24_Iterative_Split";
  inDir[7] = "Jan_24_Iterative_Alternate_Split";
  inDir[8] = "Jan_24_Simple_Split";
  inDir[9] = "Jan_24_pp_Iterative_Split";
  inDir[10] = "Jan_24_pp_Iterative_Alternate_Split";
  inDir[11] = "Jan_24_pp_Simple_Split";
  */

  std::string outDir[2];
  outDir[0] = "PbPb";
  outDir[1] = "pp";

  const int n1DPlots = 2;
  TH1D *Eff[nFiles][n1DPlots], *Fake[nFiles][n1DPlots], *EffClos[nFiles][n1DPlots], *FakeClos[nFiles][n1DPlots];
  const int n2DPlots = 2;
  TH2D *Eff2[nFiles][n2DPlots], *Fake2[nFiles][n2DPlots], *EffClos2[nFiles][n2DPlots], *FakeClos2[nFiles][n2DPlots];

  for(int i = 0; i<nFiles; i++)
  {
    TFile * f = TFile::Open(Form("../Corrections/%s/Closure.root",inDir[i].c_str()));
    for(int j = 0; j<n1DPlots; j++)
    {
      Eff[i][j] = (TH1D*)f->Get(Form("Efficiency_%d",j==0?0:2));
      Eff[i][j]->SetDirectory(0);
      Eff[i][j]->GetYaxis()->SetTitle("Efficiency");
      EffClos[i][j] = (TH1D*)f->Get(Form("EfficiencyClosure_%d",j==0?0:2));
      EffClos[i][j]->SetDirectory(0);
      EffClos[i][j]->GetYaxis()->SetTitle("Eff. Closure");
      Fake[i][j] = (TH1D*)f->Get(Form("Fake_%d",j==0?0:2));
      Fake[i][j]->SetDirectory(0);
      Fake[i][j]->GetYaxis()->SetTitle("1/(1-Fake)");
      FakeClos[i][j] = (TH1D*)f->Get(Form("FakeClosure_%d",j==0?0:2));
      FakeClos[i][j]->SetDirectory(0);
      FakeClos[i][j]->GetYaxis()->SetTitle("Fake Closure");
    }
    for(int j = 0; j<n2DPlots; j++)
    {
      Eff2[i][j] = (TH2D*)f->Get(Form("Efficiency_%d",j==0?7:8));
      Eff2[i][j]->SetDirectory(0);
      Eff2[i][j]->GetYaxis()->SetTitle("p_{T}");
      EffClos2[i][j] = (TH2D*)f->Get(Form("EfficiencyClosure_%d",j==0?7:8));
      EffClos2[i][j]->SetDirectory(0);
      EffClos2[i][j]->GetYaxis()->SetTitle("p_{T}");
      Fake2[i][j] = (TH2D*)f->Get(Form("Fake_%d",j==0?7:8));
      Fake2[i][j]->SetDirectory(0);
      Fake2[i][j]->GetYaxis()->SetTitle("p_{T}");
      FakeClos2[i][j] = (TH2D*)f->Get(Form("FakeClosure_%d",j==0?7:8));
      FakeClos2[i][j]->SetDirectory(0);
      FakeClos2[i][j]->GetYaxis()->SetTitle("p_{T}");
    }
    f->Close();
  }

  for(int i = 0; i<nFiles; i++)
  {
    c1->SetLogy(0);
    for(int j = 0; j<n1DPlots; j++)
    {
      if(j==0) c1->SetLogx();
      else     c1->SetLogx(0);
      Eff[i][j]->GetYaxis()->SetRangeUser(0,1);
      if(j==0)Eff[i][j]->GetXaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      Eff[i][j]->Draw();
      c1->SaveAs(Form("plots/May12_pp/%s%d__Efficiency.png",inDir[i].c_str(),j));
      c1->SaveAs(Form("plots/May12_pp/%s%d__Efficiency.pdf",outDir[i].c_str(),j));
      EffClos[i][j]->GetYaxis()->SetRangeUser(0.8,1.2);
      if(j==0)EffClos[i][j]->GetXaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      EffClos[i][j]->Draw();
      c1->SaveAs(Form("plots/May12_pp/%s%d__EfficiencyClosure.png",inDir[i].c_str(),j));
      c1->SaveAs(Form("plots/May12_pp/%s%d__EfficiencyClosure.pdf",outDir[i].c_str(),j));
      Fake[i][j]->GetYaxis()->SetRangeUser(0.8,1.2);
      if(j==0)Fake[i][j]->GetXaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      Fake[i][j]->Draw();
      c1->SaveAs(Form("plots/May12_pp/%s%d__Fake.png",inDir[i].c_str(),j));
      c1->SaveAs(Form("plots/May12_pp/%s%d__Fake.pdf",outDir[i].c_str(),j));
      FakeClos[i][j]->GetYaxis()->SetRangeUser(0.8,1.2);
      if(j==0)FakeClos[i][j]->GetXaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      FakeClos[i][j]->Draw();
      c1->SaveAs(Form("plots/May12_pp/%s%d__FakeClosure.png",inDir[i].c_str(),j));
      c1->SaveAs(Form("plots/May12_pp/%s%d__FakeClosure.pdf",outDir[i].c_str(),j));
    }
    c1->Clear();
    c1->SetLogx(0);
    c1->SetLogy();
    for(int j = 0; j<n2DPlots; j++)
    {
      Eff2[i][j]->GetZaxis()->SetRangeUser(0,1);
      Eff2[i][j]->GetYaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      if(j==0)Eff2[i][j]->GetXaxis()->SetRangeUser(-etaLimit,etaLimit);
      Eff2[i][j]->Draw("colz");
      c1->SaveAs(Form("plots/May12_pp/%s%d__Efficiency.png",inDir[i].c_str(),j+7));
      c1->SaveAs(Form("plots/May12_pp/%s%d__Efficiency.pdf",outDir[i].c_str(),j+7));
      EffClos2[i][j]->GetZaxis()->SetRangeUser(0.8,1.2);
      EffClos2[i][j]->GetYaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      if(j==0)EffClos2[i][j]->GetXaxis()->SetRangeUser(-etaLimit,etaLimit);
      EffClos2[i][j]->Draw("colz");
      c1->SaveAs(Form("plots/May12_pp/%s%d__EfficiencyClosure.png",inDir[i].c_str(),j+7));
      c1->SaveAs(Form("plots/May12_pp/%s%d__EfficiencyClosure.pdf",outDir[i].c_str(),j+7));
      Fake2[i][j]->GetZaxis()->SetRangeUser(1,1.2);
      Fake2[i][j]->GetYaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      if(j==0)Fake2[i][j]->GetXaxis()->SetRangeUser(-etaLimit,etaLimit);
      Fake2[i][j]->Draw("colz");
      c1->SaveAs(Form("plots/May12_pp/%s%d__Fake.png",inDir[i].c_str(),j+7));
      c1->SaveAs(Form("plots/May12_pp/%s%d__Fake.pdf",outDir[i].c_str(),j+7));
      FakeClos2[i][j]->GetZaxis()->SetRangeUser(0.8,1.2);
      FakeClos2[i][j]->GetYaxis()->SetRangeUser(minPtLimit,maxPtLimit);
      if(j==0)FakeClos2[i][j]->GetXaxis()->SetRangeUser(-etaLimit,etaLimit);
      FakeClos2[i][j]->Draw("colz");
      c1->SaveAs(Form("plots/May12_pp/%s%d__FakeClosure.png",inDir[i].c_str(),j+7));
      c1->SaveAs(Form("plots/May12_pp/%s%d__FakeClosure.pdf",outDir[i].c_str(),j+7));
    }
  }
}
