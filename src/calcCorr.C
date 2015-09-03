#include "TH1D.h"
#include "Settings.h"
#include "makeSkim.C"

void calcCorr(int job, int effVsFake)
{
  Settings s("TrkCorrInputFile.txt");
  if(effVsFake==0) makeSkim(s,job,"Eff");

//will need to do something with 'job number' for fake events eventually depending on the condor implementation... 
  if(effVsFake==1) makeSkim(s,job,"Fake");
}

int main(int argc, const char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: calculateEff <job#> <eff vs fake (0 or 1)>" << std::endl;
    return 1;
  }

  int job = std::atoi(argv[1]);
  int effVsFake = std::atoi(argv[2]);
  calcCorr(job,effVsFake);
 
  return 0;
}
