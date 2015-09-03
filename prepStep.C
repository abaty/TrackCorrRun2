#include "src/Settings.h"
#include "src/getWeights.C"

void prepStep()
{
  Settings test("TrkCorrInputFile.txt");
  produceWeights(test);
}
