#include <string>
#include <iostream>

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

class Chi2Corrector_PbPb
{
  public:
    Chi2Corrector_PbPb();
    double getChi2Scale(int cent, double pt);
    virtual ~Chi2Corrector_PbPb();

  private:
  
    TFile * file;
    TH1D *h_0_5;
    TH1D *h_5_10;
    TH1D *h_10_30;
    TH1D *h_30_50;
    TH1D *h_50_70;
    TH1D *h_70_90;
    TF1 *f_0_5;
    TF1 *f_5_10;
    TF1 *f_10_30;
    TF1 *f_30_50;
    TF1 *f_50_70;
    TF1 *f_70_90;
};


Chi2Corrector_PbPb::Chi2Corrector_PbPb()
{
  file = new TFile("PlotPlotChi2Scaling_PbPb.root");

  h_0_5 = (TH1D*)file->Get("hScaling_0_5");
  h_5_10 = (TH1D*)file->Get("hScaling_5_10");
  h_10_30 = (TH1D*)file->Get("hScaling_10_30");
  h_30_50 = (TH1D*)file->Get("hScaling_30_50");
  h_50_70 = (TH1D*)file->Get("hScaling_50_70");
  h_70_90 = (TH1D*)file->Get("hScaling_70_90");

  f_0_5 = (TF1*)file->Get("fit_0_5_lemma");
  f_5_10 = (TF1*)file->Get("fit_5_10_lemma");
  f_10_30 = (TF1*)file->Get("fit_10_30_lemma");
  f_30_50 = (TF1*)file->Get("fit_30_50_lemma");
  f_50_70 = (TF1*)file->Get("fit_50_70_lemma");
  f_70_90 = (TF1*)file->Get("fit_70_90_lemma");
}

Chi2Corrector_PbPb::~Chi2Corrector_PbPb()
{
/*
  delete h_0_5;
  delete h_5_10;
  delete h_10_30;
  delete h_30_50;
  delete h_50_70;

  delete f_0_5;
  delete f_5_10;
  delete f_10_30;
  delete f_30_50;
  delete f_50_70;
*/
  delete file;
}

double Chi2Corrector_PbPb::getChi2Scale(int cent, double pt)
{
  if(cent<10) {
    if(pt<10)
      return h_0_5->GetBinContent(h_0_5->FindBin(pt));
    else
      return f_0_5->Eval(60.); //constant value
  }
  if(cent>=10 && cent<20) {
    if(pt<20)
      return h_5_10->GetBinContent(h_5_10->FindBin(pt));
    else
      return f_5_10->Eval(60.); //constant value
  }
  if(cent>=20 && cent<60) {
    if(pt<20)
      return h_10_30->GetBinContent(h_10_30->FindBin(pt));
    else
      return f_10_30->Eval(60.); //constant value
  }
  if(cent>=60 && cent<100) {
    if(pt<20)
      return h_30_50->GetBinContent(h_30_50->FindBin(pt));
    else
      return f_30_50->Eval(60.); //constant value
  }
  if(cent>=100 && cent<140) {
    if(pt<20)
      return h_50_70->GetBinContent(h_50_70->FindBin(pt));
    else
      return f_50_70->Eval(60.); //constant value
  }
  if(cent>=140 && cent<180) {
    if(pt<20)
      return h_70_90->GetBinContent(h_70_90->FindBin(pt));
    else
      return f_70_90->Eval(60.); //constant value
  }
  if(cent>=180) { //use 70-90%
    if(pt<20)
      return h_70_90->GetBinContent(h_70_90->FindBin(pt));
    else
       return f_70_90->Eval(60.); //constant value
  }

  return 1.;
}

