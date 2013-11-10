//////////////////////////////////////////////////////////////////////////////
/////  AskConventions.h        ASK Convetnions                       /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A handy file full of enumerations and the like identifying     /////
/////     the conventions we try and use in the ASK offline code.      /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////


#ifndef ASKCONVENTIONS_H
#define ASKCONVENTIONS_H

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

#include "PhysicalConstants.h"

namespace AskCons {

   enum WhichFlux_t {
      kStandardGzk = 1,
      kMonoEnergetic = 2
   };


    enum WorldModelMaterial_t {
	kWMFirn = 0x01,
	kWMIce,
	kWMWater,
	kWMSoftSed,
	kWMHardSed,
	kWMUpperCrust,
	kWMMiddleCrust,
	kWMLowerCrust,
	kWMInside, //Might improve at future date
	kWMStdRock
    };


    enum WorldModelType_t { //Absolute values are relative to centre of Earth
	kWMGeoid = 0x01,
	kWMSurface, //actual surface (relative to geoid)
	kWMSurfaceAbsolute,
	kWMFirnBottom,
	kWMFirnBottomAbsolute,
	kWMIceBottom,
	kWMIceBottomAbsolute,
	kWMWaterBottom,
	kWMWaterBottomAbsolute,
	kWMSoftSedBottom,
	kWMSoftSedBottomAbsolute,
	kWMHardSedBottom,
	kWMHardSedBottomAbsolute,
	kWMUpperCrustBottom,
	kWMUpperCrustBottomAbsolute,
	kWMMiddleCrustBottom,
	kWMMiddleCrustBottomAbsolute,
	kWMLowerCrustBottom,
	kWMLowerCrustBottomAbsolute,
	kWMFirnDepth,
	kWMIceThickness,
	kWMWaterDepth
    };
	



    enum BedmapType_t {
	kBathymetry = 0x01,
	kBedElevation,
	kGroundBed,
	kIceThickness,
	kSurfaceHeight,
	kWaterDepth
    };

    enum Crust2Type_t {
	kTopOfIce=0x01,
	kTopOfIceAbsolute,
	kBottomOfIce,
	kBottomOfWater,
	kBottomOfSoftSed,
	kBottomOfHardSed,
	kBottomOfUpperCrust,
	kBottomOfMiddleCrust,
	kBottomOfLowerCrust,
	kGeoid,
	kElevation,
	kWaterDensity,
	kIceDensity,
	kSoftSedDensity,
	kHardSedDensity,
	kUpperCrustDensity,
	kMiddleCrustDensity,
	kLowerCrustDensity,
	kThicknessOfIce,
	kThicknessOfWater,
	kThicknessOfSoftSed,
	kThicknessOfHardSed,
	kThicknessOfUpperCrust,
	kThicknessOfMiddleCrust,
	kThicknessOfLowerCrust,
	kTotalCrustThickness
    };

    enum InteractionType_t {
	kNC       = 0x01,
	kCC,
	kBrem,
	kPair,
	kPhoto,
	kKnockOn,
	kDecay
    };
    
    
    enum ParticleType_t {
	kNothing  = 0,
	kElectron = 11,
	kNue      = 12,
	kMuon     = 13,
	kNumu     = 14,
	kTau      = 15,
	kNutau    = 16,
	kPositron = -11,
	kEbar     = -11,
	kANue     = -12,
	kNuebar   = -12,
	kMuonbar  = -13,
	kAMuon    = -13,
	kNumubar  = -14,
	kANumu    = -14,
	kTaubar   = -15,
	kATau     = -15,
	kNutaubar = -16,
	kANutau   = -16	
    };

    inline bool isAChargedLepton(ParticleType_t type)
    {
	switch (type) {
	    case kMuon:
	    case kTau:
	    case kMuonbar:
	    case kTaubar:
		return true;
	    default:
		return false;
	}
    }

    inline double isItATau(ParticleType_t type)
    {
	switch (type) {
	    case kTau:
	    case kTaubar:
		return 1;
	    default:
		return 0;
	}
    }

    inline ParticleType_t getCCPartner(ParticleType_t type)
    {    
	switch (type) {
	    case kElectron: return kNue;
	    case kMuon: return kNumu;
	    case kTau: return kNutau;
	    case kNue: return kElectron;
	    case kNumu: return kMuon;
	    case kNutau: return kTau;
	    case kEbar: return kNuebar;
	    case kMuonbar: return kNumubar;
	    case kTaubar: return kNutaubar;
	    case kNuebar: return kEbar;
	    case kNumubar: return kMuonbar;
	    case kNutaubar: return kTaubar;
	    default : return kNothing;
	}	    	    
    }

    inline double getMass(ParticleType_t type) 
    {  
	switch (type) {
	    case kElectron: return electronMass;
	    case kMuon: return muonMass;
	    case kTau: return tauMass;
	    case kNue: return 0;
	    case kNumu: return 0;
	    case kNutau: return 0;
	    case kEbar: return electronMass;
	    case kMuonbar: return muonMass;
	    case kTaubar: return tauMass;
	    case kNuebar: return 0;
	    case kNumubar: return 0;
	    case kNutaubar: return 0;
	    default : return 0;
	}	    	    
    }
	


    enum MaterialType_t {
	kWater    =  0,
	kRock     =  1,
	kIce      =  2,
	kFirn     =  3
    };

    inline double getZForMaterial(MaterialType_t type)
    {
	switch(type) {
	    case kWater:
		return 6.6;
	    case kRock:
		return 11;
	    case kIce:
		return 6.6;
	    case kFirn:
		return 6.6;
	    default:
		return 6.6; //Default to water
	}
    }

    inline double getAForMaterial(MaterialType_t type)
    {
	switch(type) {
	    case kWater:
		return 11.89;
	    case kRock:
		return 22;
	    case kIce:
		return 11.89;
	    case kFirn:
		return 11.89;
	    default:
		return 11.89; //Default to water
	}
    }

    inline double getDensityForMaterial(MaterialType_t type)
    {
	switch(type) {
	    case kWater:
		return 1;
	    case kRock:
		return 2.2;
	    case kIce:
		return 0.92;
	    case kFirn:
		return 0.92; //Will change to better value at some point
	    default:
		return 1; //Default to water
	}
    }

 


}




#endif //ASKCONVENTIONS_H
