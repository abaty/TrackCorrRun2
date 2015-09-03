#include "Settings.h"
#include "getWeights.C"

void prepStep()
{
  Settings test("../TrkCorrInputFile.txt");
  produceWeights(test);
}
