///////////////////////////////////////////////////////////////////////////////
/////   RFPulse.cxx
/////   Describes the location and characteristics of an RFPulse within
/////   an RFStation
//////////////////////////////////////////////////////////////////////////////
#include "RFPulse.h"
#include "RFPulse.h"

#include <iostream>

ClassImp(RFPulse)

RFPulse::RFPulse()
{
//Default (zero) constructor 
  
}

RFPulse::~RFPulse()
{
  //Default destructor
}

RFPulse::RFPulse(TVector3 location, TVector3 direction)
  :fLocation(location),fDirection(direction)
{
  //Location asignment constructor
}

