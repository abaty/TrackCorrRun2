#include "Settings.h"
#include "getWeights.C"

void prepStep()
{
  Settings s("TrkCorrInputFile.txt");
  produceWeights(s);
}

int main()
{
  prepStep();
  return 1;
}
