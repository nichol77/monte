///////////////////////////////////////////////////////////////////////////////
/////   UtilityFuncs.cxx
/////  Namespace for holding utility functions
//////////////////////////////////////////////////////////////////////////////

#include "UtilityFuncs.h"
#include "TMath.h"

Double_t UtilityFuncs::linearInterpolation(Double_t x, Double_t *iX, Double_t *iY, Int_t numPoints)
{
    Int_t i1,i2;
    
    if(x>iX[numPoints-1]) {
	i1=numPoints-2;
	i2=numPoints-1;
    }
    else if(x<iX[0]) {
	i1=0;
	i2=1;    
    } 
    else {
	i1=TMath::BinarySearch(numPoints,iX,x);
	i2=i1+1;
    }
    Double_t m=(iY[i1]-iY[i2])/(iX[i1]-iX[i2]);
    Double_t retVal=iY[i1]-m*(iX[i1]-x);
    
    return retVal;
    
            	    
}
