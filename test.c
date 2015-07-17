#include "src/Settings.h"
#include "src/getWeights.C"

void test()
{
  Settings test("TrkCorrInputFile.txt");
  produceWeights(test);
}
