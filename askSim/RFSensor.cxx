///////////////////////////////////////////////////////////////////////////////
/////   RFSensor.cxx
/////   Describes the location and characteristics of an RFSensor within
/////   an RFStation
//////////////////////////////////////////////////////////////////////////////
#include "RFSensor.h"
#include "RFPulse.h"

#include <iostream>

ClassImp(RFSensor)

RFSensor::RFSensor()
{
//Default (zero) constructor 
  
}

RFSensor::~RFSensor()
{
  //Default destructor
}

RFSensor::RFSensor(TVector3 relPos)
  :fRelPos(relPos)
{
  //Location asignment constructor
}


Bool_t RFSensor::testTrigger(RFPulse  &pulsey)
{
  std::cerr << "Using Dummy Test Trigger Routine\n";
}
