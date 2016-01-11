#include "TrkSettings.h"
#include "getWeights.C"

void prepStep()
{
  TrkSettings s("TrkCorrInputFile.txt");
  produceWeights(s);
}

int main()
{
  prepStep();
  return 1;
}
