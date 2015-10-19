#include "getTrkCorr.h"
#include <iostream>

int example()
{
 TrkCorr t;
 
 float pt[1000];
 float eta[1000];
 float phi[1000];
 bool highPurity[1000];
 int nTrk = 100;

 for(int i = 0; i<1000; i++)
 {
   pt[i] = 0.5+i/10.0;
   eta[i] = -2+i/100.0;
   phi[i] = -2+i/100.0;
   highPurity[i] = true;
 }

 t.UpdateEventInfo(pt,eta,phi,highPurity,101);
 std::cout <<  t.getTrkCorr(0.5,1,-2,30) << std::endl;
}
