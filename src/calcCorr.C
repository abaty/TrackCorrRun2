#include "TH1D.h"
#include "Settings.h"
#include "iterate.C"
#include "makeSkim.C"
#include <vector>
#include <iostream>

void calcCorr(int job, int effOrFake)
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  //***********************************************************************************************
  //pick up Global settings and then set local process settings
  Settings s("TrkCorrInputFile.txt");
  s.job = job;

  //figuring out which subset of data to skim based on the job number
  int jobTemp = s.job;
  for(int i = 0; i<s.nPtBinCoarse; i++)
  {
    if(jobTemp>=s.nCentPUBinCoarse.at(i)) jobTemp = jobTemp-s.nCentPUBinCoarse.at(i);
    else
    { 
      s.nSkip = s.eventSkip.at(i).at(jobTemp);
      s.ptMin = s.ptBinCoarse.at(i);
      s.ptMax = s.ptBinCoarse.at(i+1);
      s.centPUMin = s.centPUBinCoarse.at(i).at(jobTemp);
      s.centPUMax = s.centPUBinCoarse.at(i).at(jobTemp+1);
      break;
    }
  }

  //************************************************************************************************
  //make skims needed from main data files to parse down info
  //will need to do something with 'job number' for fake events eventually depending on the condor implementation... 
  if(effOrFake==0) makeSkim(s,"Eff");
  if(effOrFake==1) makeSkim(s,"Fake");

  //***********************************************************************************************
  // iteration procedure
  std::cout << "Beginning iteration procedure..." << std::endl;
  for(int i = 0; i<3 ; i++)
  {
    std::cout << "Starting iteration " << i << std::endl;
    //stepTypes are pt, acceptance, centrality, jet pt, eta, rmin, occupancy
    int stepType = s.stepOrder.at(i%s.nStep); 
    if(effOrFake==0) iterate(s,i,stepType,"Eff");
    if(effOrFake==1) iterate(s,i,stepType,"Fake");

    //check to see if we are done yet 
    if(i>=s.nStep*s.fullIterations && s.terminateStep==stepType) break;
  } 
}

//***********************************************************************************************
int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: calcCorr <job#> <eff vs fake (0 or 1)>" << std::endl;
    return 1;
  }

  int job = std::atoi(argv[1]);
  int effVsFake = std::atoi(argv[2]);
  calcCorr(job,effVsFake);
 
  return 0;
}
