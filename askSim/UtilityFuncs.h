///////////////////////////////////////////////////////////////////////////////
/////   UtilityFuncs.h
/////  Namespace for holding utility functions
//////////////////////////////////////////////////////////////////////////////

#ifndef UTILITYFUNCS_H
#define UTILITYFUNCS_H

#include "TObject.h"

#include "AskConventions.h"


namespace UtilityFuncs  {
    //UtilityFunction
    Double_t linearInterpolation(Double_t x, Double_t *iX, Double_t *iY, Int_t numPoints);
};


#endif //UTILITYFUNCS_H
