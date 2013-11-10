///////////////////////////////////////////////////////////////////////////////
/////   WorldModel.h
/////   Is a model of the world, strangely enough
/////   The Crust2 stuff is taken from Amy's UCLA sim,
/////   the bedmap stuff I added.
/////   rjn@mps.ohio-state.edu
/////   4th May 2005
//////////////////////////////////////////////////////////////////////////////

#ifndef WORLDMODEL_H
#define WORLDMODEL_H

#include "TObject.h"
#include "TVector3.h"

#include "AskConventions.h"
#include "AskGeom.h"
#include "BedmapTable.h"

class WorldModel : public TObject {
public:
    WorldModel(Int_t useBedmap=1);
    ~WorldModel();
 
  //Instance generator
  static WorldModel*  Instance(Int_t useBedmap=1);

    Double_t getCrust2Data(Double_t x, Double_t y, AskCons::Crust2Type_t type);
    Double_t getCrust2Data(Double_t x, Double_t y, Double_t z, AskCons::Crust2Type_t type);
   
//    Double_t getData(Double_t x, Double_t y, AskCons::WorldModelType_t type);
//    AskCons::WorldModelMaterial_t getMaterial(Double_t point[3]);
    Double_t getSurface(Double_t x, Double_t y, Double_t z=6e6);
    Double_t getSurface(TVector3 &thePos);
    Double_t getIceThickness(Double_t x, Double_t y, Double_t z=6e6);
    Double_t getIceThickness(TVector3 &thepos);

    void getSurfaceNormal(Double_t x, Double_t y, TVector3 *normal);

 protected:
    static WorldModel *fgInstance;  
    

private:

    void loadCrustData();
    void defineGeoid();
    void initBedMap();

    Int_t fUseBedmap;

    BedmapTable *fBedMapBedElevation;
    BedmapTable *fBedMapWaterDepth;
    BedmapTable *fBedMapIceThickness;



// parameters of the crust 2 earth model
    double fSurface[crust2::N_LON][crust2::N_LAT]; // elevation at the surface (top of ice) compared to geoid (in meters)
    double fIceBottom[crust2::N_LON][crust2::N_LAT];  // elevation at the *bottom* of ice layer (in meters)
    double fWaterBottom[crust2::N_LON][crust2::N_LAT]; // elevation at the bottom of water layer (in meters)
    double fSoftSedBottom[crust2::N_LON][crust2::N_LAT]; // elevation at the bottom of soft set layer (in meters)
    double fHardSedBottom[crust2::N_LON][crust2::N_LAT]; // elev at bottom of hard sed layer (in meters)
    double fUpperCrustBottom[crust2::N_LON][crust2::N_LAT]; // elev at bottom of upper crust layer (in meters)
    double fMiddleCrustBottom[crust2::N_LON][crust2::N_LAT]; // elev at bottom of middle crust layer (in meters)
    double fLowerCrustBottom[crust2::N_LON][crust2::N_LAT]; // elev at bottom of lower crust layer (in meters)
    double fGeoid[crust2::N_LAT];  // realistic shape of earth-radius at each latitude (in meters)
    double fElevation[crust2::N_LON][crust2::N_LAT]; // If no water, measures the elevation (relative to geoid, in meters) of the top of the ice or rock (i.e., air interface).  If there is water, measures elevation to bottom of water. (There may or may not be ice on top of the water.) 
    
    double fWaterDensity[crust2::N_LON][crust2::N_LAT]; // density of water layer bin by bin
    double fIceDensity[crust2::N_LON][crust2::N_LAT]; // density of ice layer bin by bin
    double fSoftSedDensity[crust2::N_LON][crust2::N_LAT]; // density of soft sed layer
    double fHardSedDensity[crust2::N_LON][crust2::N_LAT]; // density of hard sed layer
    double fUpperCrustDensity[crust2::N_LON][crust2::N_LAT]; // density of upper crust layer
    double fMiddleCrustDensity[crust2::N_LON][crust2::N_LAT]; // density of middle crust layer
    double fLowerCrustDensity[crust2::N_LON][crust2::N_LAT]; // density of lower crust layer
    double fWaterThickness[crust2::N_LON][crust2::N_LAT];  // thickness of water layer (in km)
    double fIceThickness[crust2::N_LON][crust2::N_LAT]; // thickness of ice layer (in km)
    double fSoftSedThickness[crust2::N_LON][crust2::N_LAT]; // thickness of soft sed layer (in km)
    double fHardSedThickness[crust2::N_LON][crust2::N_LAT]; // thickness of hard sed layer (in km)
    double fUpperCrustThickness[crust2::N_LON][crust2::N_LAT]; // thickness of upper crust layer (in km)
    double fMiddleCrustThickness[crust2::N_LON][crust2::N_LAT]; // thickness of middle crust layer (in km)
    double fLowerCrustThickness[crust2::N_LON][crust2::N_LAT]; // thickness of lower crust layer (in km)
    double fTotalCrustThickness[crust2::N_LON][crust2::N_LAT]; // total thickness of crust (in km)    
    
    ClassDef(WorldModel,2);
};


#endif //WORLDMODEL_H
